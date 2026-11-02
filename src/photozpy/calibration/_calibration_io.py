"""Internal collection and output bookkeeping for calibration workflows.

Numerical calibration algorithms belong to their caller, not this module.
"""

from __future__ import annotations

import os
from pathlib import Path
import tempfile
from typing import TypeAlias

import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy.nddata import CCDData
from ccdproc import ImageFileCollection

from ..collections import ImageCollectionGroup
from ..collections.image_collection_helpers import (
    NoMatchingImagesError,
    filter_image_collection,
    refresh_image_collection,
    refresh_image_collection_group,
)


ImageInput: TypeAlias = ImageFileCollection | ImageCollectionGroup


def children(collection: ImageInput) -> tuple[ImageFileCollection, ...]:
    """Return child collections without changing their directory boundaries."""
    if isinstance(collection, ImageCollectionGroup):
        return collection.collections
    if isinstance(collection, ImageFileCollection):
        return (collection,)
    raise TypeError('Expected an ImageFileCollection or ImageCollectionGroup.')


def refresh(collection: ImageInput) -> None:
    """Refresh in place without broadening explicit file selections."""
    for child in children(collection):
        if child.ext != 0:
            raise ValueError('Calibration currently requires collection ext=0.')

    if isinstance(collection, ImageCollectionGroup):
        refresh_image_collection_group(collection)
    else:
        refresh_image_collection(collection)


def paths(collection: ImageFileCollection) -> list[Path]:
    """Resolve selected filenames in their original order; skip hidden files."""
    return [
        (Path(name) if Path(name).is_absolute()
         else Path(collection.location or '.') / name).resolve()
        for name in collection.files if not Path(name).name.startswith('.')
    ]


def select(collection: ImageFileCollection, **criteria: object) -> list[Path]:
    """Use the refactored multi-value filter; represent no matches as []."""
    try:
        return paths(filter_image_collection(collection, **criteria))
    except NoMatchingImagesError:
        return []


def is_corrected(header: fits.Header, keyword: str) -> bool:
    """Recognize existing Yes/No or FITS boolean processing flags."""
    value = header.get(keyword, False)
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in ('yes', 'true', 't', '1'):
            return True
        if normalized in ('no', 'false', 'f', '0', ''):
            return False
        raise ValueError(f'Unrecognized {keyword} flag: {value!r}.')
    if value in (True, False, 0, 1):
        return bool(value)
    raise ValueError(f'Unrecognized {keyword} flag: {value!r}.')


def exposure(header: fits.Header, path: Path) -> float:
    """Require a finite, positive EXPTIME in seconds, without guessing units."""
    try:
        value = float(header['EXPTIME'])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f'{path.name}: missing or invalid EXPTIME.') from exc

    if not np.isfinite(value) or value <= 0:
        raise ValueError(f'{path.name}: EXPTIME must be positive and finite.')
    return value


def validate_compatible(image_paths: list[Path]) -> None:
    """Check 2D geometry and compatible declared units before writing a batch."""
    expected_shape: tuple[int, int] | None = None
    expected_unit: u.UnitBase | None = None

    for path in image_paths:
        header = fits.getheader(path)
        if header.get('NAXIS') != 2:
            raise ValueError(f'{path.name}: expected a 2D primary image.')
        shape = (header['NAXIS2'], header['NAXIS1'])
        if 'BUNIT' not in header:
            raise ValueError(f'{path.name}: missing BUNIT; correct headers first.')
        # Match CCDData.read's handling of the historical uppercase ADU card
        # written by HeaderCorrection; the physical unit remains unchanged.
        unit_text = str(header['BUNIT']).strip()
        unit = u.adu if unit_text.lower() == 'adu' else u.Unit(unit_text)

        if expected_shape is not None and shape != expected_shape:
            raise ValueError(f'{path.name}: calibration image shapes differ.')
        if expected_unit is not None and unit != expected_unit:
            raise ValueError(f'{path.name}: calibration image units differ.')

        expected_shape, expected_unit = shape, unit


def check_destinations(destinations: list[Path], overwrite: bool) -> None:
    """Reject path collisions or forbidden overwrites before any file is written."""
    if len(destinations) != len(set(destinations)):
        raise ValueError('Multiple images map to the same output path.')
    for path in destinations:
        if path.exists() and not overwrite:
            raise FileExistsError(path)


def write_ccd(
    ccd: CCDData, destination: Path, *, overwrite: bool,
    cmask: np.ndarray | None = None,
) -> None:
    """Publish the caller's CCDData output atomically, without changing its data.

    CCDData controls the output HDUs, including any retained mask/uncertainty.
    This is not a copy of unrelated extensions in the original input FITS file.
    Existing checksum cards are recomputed for the newly serialized CCDData.
    Optional CMASK records cosmic-ray detections independently of CCDData.mask;
    standard CCDData readers do not use it to exclude image values.
    """
    if cmask is not None and cmask.shape != ccd.shape:
        raise ValueError('CMASK must match the output image shape.')
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(
        prefix=f'.{destination.name}.', suffix='.fits', dir=destination.parent,
    )
    os.close(descriptor)
    temporary = Path(name)

    try:
        checksum = 'CHECKSUM' in ccd.meta or 'DATASUM' in ccd.meta
        with ccd.to_hdu() as hdul:
            if cmask is not None:
                hdul.append(fits.ImageHDU(np.asarray(cmask, dtype=np.uint8), name='CMASK'))
            hdul.writeto(temporary, overwrite=True, checksum=checksum)
        if overwrite:
            os.replace(temporary, destination)
        else:
            os.link(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def rebuild(
    source: ImageInput,
    replacements: dict[Path, Path] | None = None,
    additions: dict[int, list[Path]] | None = None,
) -> ImageInput:
    """Return selections containing replacements and newly generated masters.

    Include only known input members and explicit outputs; do not rescan unrelated
    files. Files outside a child's original location retain absolute paths.
    """
    replacements = {} if replacements is None else replacements
    additions = {} if additions is None else additions
    result_children = []

    for child in children(source):
        selected = [replacements.get(path, path) for path in paths(child)]
        selected.extend(additions.get(id(child), []))
        selected = list(dict.fromkeys(selected))
        if not selected:
            continue

        parents = {path.parent for path in selected}
        location = next(iter(parents)) if len(parents) == 1 else Path(child.location or '.').resolve()
        filenames = [
            path.name if path.parent == location else str(path) for path in selected
        ]
        result_children.append(ImageFileCollection(
            location=location, filenames=filenames, keywords='*', ext=child.ext,
        ))

    if not result_children:
        raise ValueError('No images remain in the output collection.')
    if isinstance(source, ImageCollectionGroup):
        return ImageCollectionGroup.from_collections(result_children)
    return result_children[0]
