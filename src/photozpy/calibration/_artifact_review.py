"""File identity and no-overwrite rejection operations for cosmic-ray review."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import tempfile

from astropy.io import fits

from ._calibration_io import check_destinations, is_corrected


@dataclass(frozen=True)
class RejectionPlan:
    """Explicit paths for one original exposure and its known derived files."""

    raw: Path
    after: Path
    before: Path
    destination: Path

    @property
    def record(self) -> Path:
        """Return the adjacent audit record, outside the FITS input selection."""
        return self.destination.with_name(f'{self.destination.name}.rejection.json')

    @property
    def derived(self) -> tuple[Path, ...]:
        """Return only known products; never search by a filename prefix."""
        return (
            self.after, self.before,
            self.after.with_name(f'{self.after.stem}_cr_mask.png'),
        )


def plan_rejections(pairs: Mapping[Path, Path]) -> dict[Path, RejectionPlan]:
    """Match unchanged filenames plus observation metadata before opening DS9.

    Infer each root from root/analysis/image.fits. Never search recursively.
    Other layouts, renamed files, or ambiguous metadata require resolution
    before review; guessing a raw-file location is unsafe.
    """
    plans: dict[Path, RejectionPlan] = {}
    for after, before in pairs.items():
        if after.parent.name != 'analysis':
            raise ValueError('Review requires working images in root/analysis/.')
        root = after.parent.parent
        rejected = root / 'rejected'
        if rejected.is_symlink():
            raise ValueError('The rejected directory must not be a symbolic link.')
        if rejected.exists() and not rejected.is_dir():
            raise NotADirectoryError(rejected)
        raw = root / after.name
        for path in (raw, after, before):
            if path.is_symlink() or not path.is_file():
                raise ValueError(f'Expected a regular, retained FITS file: {path}')
        if raw.samefile(after) or raw.samefile(before):
            raise ValueError('Raw originals must be separate from working images and backups.')
        original, cleaned, backup = (fits.getheader(path) for path in (raw, after, before))
        for keyword in ('NAXIS1', 'NAXIS2', 'DATE-OBS', 'EXPTIME'):
            if original.get(keyword) is None or not (
                original[keyword] == cleaned.get(keyword) == backup.get(keyword)
            ):
                raise ValueError(f'{after.name}: cannot match raw exposure using {keyword}.')
        if cleaned.get('IMTYPE') != 'Light' or not is_corrected(cleaned, 'CRCORR'):
            raise ValueError(f'{after.name}: review requires a cleaned Light exposure.')
        if cleaned.get('ALIGN') in ('Reference', 'Registered'):
            raise ValueError('Review cosmic-ray outputs before registration.')
        if any(is_corrected(original, key) for key in ('BIASCORR', 'DARKCORR', 'FLATCORR', 'CRCORR')):
            raise ValueError(f'{raw}: expected an uncalibrated original exposure.')
        plans[after] = RejectionPlan(raw, after, before, rejected / raw.name)

    products = [path for plan in plans.values() for path in (plan.raw, *plan.derived)]
    if len(products) != len(set(products)):
        raise ValueError('Exposure paths overlap; rejection is ambiguous.')
    destinations = [path for plan in plans.values() for path in (plan.destination, plan.record)]
    if set(products).intersection(destinations):
        raise ValueError('Rejected destinations overlap exposure inputs or products.')
    check_destinations(destinations, overwrite=False)
    return plans


def retain_original(plan: RejectionPlan, copy: Callable[[Path, Path], None]) -> None:
    """Archive raw bytes and an audit record before removing the root original.

    Existing archives and records are never replaced. If publication fails,
    remove only this call's new artifacts and leave the root original intact.
    Derived-file cleanup is the caller's responsibility after this succeeds.
    """
    copied = recorded = False
    temporary: Path | None = None
    try:
        copy(plan.raw, plan.destination)
        copied = True
        with tempfile.NamedTemporaryFile(
            mode='w', encoding='utf-8', dir=plan.destination.parent,
            prefix='.', suffix='.json', delete=False,
        ) as output:
            temporary = Path(output.name)
            json.dump({
                'stage': 'post_cosmic_ray_review',
                'reason': 'Original exposure rejected during visual quality review.',
                'timestamp_utc': datetime.now(timezone.utc).isoformat(),
                'original': str(plan.raw), 'rejected': str(plan.destination),
                'cleaned': str(plan.after), 'before': str(plan.before),
            }, output, indent=2)
            output.write('\n')
        os.link(temporary, plan.record)
        recorded = True
        plan.raw.unlink()
    except BaseException:
        # An interrupt can arrive immediately after unlink succeeded. Never
        # remove the retained original once its root path has disappeared.
        if plan.raw.exists():
            if recorded:
                plan.record.unlink()
            if copied:
                plan.destination.unlink()
        raise
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
