"""Adapt the legacy zero-point recipe to independent source result slots.

The numerical main line follows MagnitudeInfo.calculate_zero_points,
remove_outlier and calibrate_inst_mags from the pre-refactoring implementation:
reference minus instrumental, one z-score pass, arithmetic mean, then addition
to target instrumental magnitudes with unchanged instrumental errors.
"""

from __future__ import annotations

from dataclasses import dataclass
import logging
from pathlib import Path
from typing import TYPE_CHECKING

import astropy.units as u
import numpy as np
import pandas as pd
from scipy import stats

from .filter_magnitude import FilterMagnitude, MagnitudeSystem, UPPER_LIMIT_ERROR

if TYPE_CHECKING:
    from matplotlib.figure import Figure

    from ..telescope import Telescope
    from .sources import Sources


logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class _ZeroPointStatistics:
    """Keep one inspected calculation and the inputs that produced it."""

    telescope: Telescope
    magnitude_system: MagnitudeSystem
    summary: pd.DataFrame
    details: pd.DataFrame
    signature: tuple[object, ...]


def standard_signature(sources: Sources) -> tuple[object, ...]:
    """Snapshot standard values, including edits through shared QTables.

    Normalize NaN to None so unchanged missing values compare equal. No hashes,
    persistent cache or serialization are needed for this in-memory result.
    """
    signature: list[object] = []
    for single in sources.by_role("standard"):
        info = single.source_info
        measurements: list[object] = []
        for kind in ("instrumental", "reference"):
            measurement = single.get_magnitudes(kind=kind)
            if measurement is None:
                measurements.append(None)
                continue
            values = tuple(
                tuple(None if np.isnan(value) else float(value) for value in column)
                for column in (
                    measurement.magnitudes.to_value(u.mag),
                    measurement.errors.to_value(u.mag),
                    measurement.significance,
                )
            )
            measurements.append(
                (
                    id(measurement.telescope),
                    measurement.magnitude_system,
                    measurement.filter_names,
                    values,
                )
            )
        signature.append(
            (
                info.observation_field.field_name,
                info.source_name,
                id(info.observation_field.telescope),
                tuple(measurements),
            )
        )
    return tuple(signature)


def calculate_zero_points(
    sources: Sources, sigma: float, min_significance: float
) -> _ZeroPointStatistics:
    """Calculate per-star offsets and the legacy unweighted clipped mean.

    The population standard deviation and single-pass z-scores match scipy's
    defaults in the old implementation. Constant and singleton samples are
    retained explicitly: their undefined z-scores formerly rejected every star.
    Scatter is descriptive, with no inferred uncertainty of the mean.
    """
    if not np.isfinite(sigma) or sigma <= 0:
        raise ValueError("sigma must be finite and positive.")
    if not np.isfinite(min_significance) or min_significance < 0:
        raise ValueError("min_significance must be finite and non-negative.")
    standards = sources.by_role("standard")
    if not standards:
        raise ValueError("At least one standard-role source is required.")
    telescope = standards[0].source_info.observation_field.telescope
    systems: set[MagnitudeSystem] = set()
    rows: list[dict[str, object]] = []

    for single in standards:
        info = single.source_info
        if info.observation_field.telescope is not telescope:
            raise ValueError("Standards must share the same Telescope instance.")
        instrumental = single.get_magnitudes(kind="instrumental")
        reference = single.get_magnitudes(kind="reference")
        if reference is not None:
            systems.add(reference.magnitude_system)
        for index, filter_name in enumerate(telescope.filters()):
            m_inst = (
                float(instrumental.magnitudes[index].to_value(u.mag))
                if instrumental is not None
                else np.nan
            )
            m_ref = (
                float(reference.magnitudes[index].to_value(u.mag))
                if reference is not None
                else np.nan
            )
            significance = (
                float(instrumental.significance[index])
                if instrumental is not None
                else np.nan
            )
            upper_limit = bool(
                instrumental is not None
                and instrumental.errors[index] == UPPER_LIMIT_ERROR
            )
            zero_point = np.nan
            if instrumental is None:
                reason = "missing_instrumental"
            elif reference is None:
                reason = "missing_reference"
            elif not np.isfinite(m_inst) or not np.isfinite(m_ref):
                reason = "missing_magnitude"
            elif upper_limit:
                reason = "upper_limit"
            else:
                zero_point = m_ref - m_inst
                if not np.isfinite(significance):
                    reason = "missing_significance"
                elif significance < min_significance:
                    reason = "low_significance"
                else:
                    reason = "retained"
            rows.append(
                dict(
                    source_id=info.source_name,
                    field=info.observation_field.field_name,
                    filter=filter_name,
                    reference_magnitude=m_ref,
                    instrumental_magnitude=m_inst,
                    significance=significance,
                    zero_point=zero_point,
                    used=reason == "retained",
                    reason=reason,
                )
            )

    if len(systems) != 1:
        raise ValueError("Standards require one common AB or Vega reference system.")
    magnitude_system = next(iter(systems))
    details = pd.DataFrame(rows)
    summaries: list[dict[str, object]] = []
    for filter_name in telescope.filters():
        candidates = details.index[(details["filter"] == filter_name) & details["used"]]
        zero = details.loc[candidates, "zero_point"].to_numpy(dtype=float)
        # Preserve remove_outlier's single-pass scipy.stats.zscore recipe.
        keep = np.ones(len(zero), dtype=bool)
        if len(zero) > 1 and not np.all(zero == zero[0]):
            zscores = stats.zscore(zero)
            keep = np.abs(zscores) <= sigma
        rejected = candidates[~keep]
        details.loc[rejected, "used"] = False
        details.loc[rejected, "reason"] = "sigma_clipped"
        retained = zero[keep]
        summaries.append(
            dict(
                filter=filter_name,
                zero_point=float(np.mean(retained)) if len(retained) else np.nan,
                scatter=float(np.std(retained)) if len(retained) else np.nan,
                n_candidates=len(zero),
                n_used=len(retained),
                n_clipped=int(np.count_nonzero(~keep)),
                magnitude_system=magnitude_system,
            )
        )
    summary = pd.DataFrame(summaries)
    if not summary["n_used"].any():
        raise ValueError("No usable standard-star zero points remain after selection.")
    if len(standards) <= 8:
        logger.warning(
            "Only %d standard stars are available; inspect per-filter counts and scatter.",
            len(standards),
        )
    return _ZeroPointStatistics(
        telescope, magnitude_system, summary, details, standard_signature(sources)
    )


def plot_zero_point_statistics(
    result: _ZeroPointStatistics,
    reference_offset: float,
    save_path: str | Path | None,
) -> Figure:
    """Display all finite offsets and distinguish retained/rejected samples.

    Count unavailable offsets separately; missing measurements cannot be plotted
    as zero points. Matplotlib's autoscaling includes both retained and rejected
    offsets, even when a rejected point is far outside the visual guide lines.
    """
    import matplotlib.pyplot as plt

    if not np.isfinite(reference_offset) or reference_offset < 0:
        raise ValueError("reference_offset must be finite and non-negative.")
    details = result.details
    filters = [
        name
        for name in result.summary["filter"]
        if np.isfinite(details.loc[details["filter"] == name, "zero_point"]).any()
    ]
    ncols = min(2, len(filters))
    nrows = (len(filters) + ncols - 1) // ncols
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(6 * ncols, 4.5 * nrows), squeeze=False
    )
    for ax, filter_name in zip(axes.flat, filters):
        stars = details[details["filter"] == filter_name]
        retained = stars[stars["used"]]
        excluded = stars[~stars["used"] & np.isfinite(stars["zero_point"])]
        n_unavailable = len(stars) - len(retained) - len(excluded)
        # Keep the legacy vertical scatter layout, with one panel per filter.
        ax.scatter(
            np.zeros(len(retained)),
            retained["zero_point"],
            label="retained offsets",
            color="skyblue",
        )
        if len(excluded):
            ax.scatter(
                np.zeros(len(excluded)),
                excluded["zero_point"],
                label="rejected offsets",
                color="darkorange",
                marker="x",
                s=55,
                linewidths=1.5,
                zorder=4,
            )
        row = result.summary[result.summary["filter"] == filter_name].iloc[0]
        if np.isfinite(row["zero_point"]):
            mean = row["zero_point"]
            ax.scatter([0], [mean], color="lime", label="clipped average", zorder=3)
            ax.axhline(
                mean - reference_offset,
                color="red",
                label=f"average +/- {reference_offset:g} mag (guide)",
            )
            ax.axhline(mean + reference_offset, color="red")
        counts = f"{len(retained)} retained, {len(excluded)} rejected"
        if n_unavailable:
            counts += f", {n_unavailable} unavailable"
        ax.set_title(f"{filter_name}\n{counts}")
        ax.set_ylabel("Zero point [mag]")
        ax.set_xticks([0], [filter_name])
        ax.set_xlim(-0.1, 1)
        ax.legend()
    for ax in list(axes.flat)[len(filters):]:
        ax.set_visible(False)
    fig.tight_layout()
    if save_path is not None:
        fig.savefig(save_path, dpi=300)
    plt.show()
    return fig


def calibrate_targets(sources: Sources, result: _ZeroPointStatistics) -> pd.DataFrame:
    """Stage all science results, then replace only their calibrated slots.

    No error propagation is added: each instrumental error and significance is
    copied unchanged. Finite upper-limit magnitudes receive the same additive
    offset and keep their -99 error flag. No exposure-time normalization occurs.
    """
    targets = sources.by_role("science")
    if not targets:
        raise ValueError("At least one science-role source is required.")
    zero_points = result.summary.set_index("filter")["zero_point"]
    pending: list[tuple[Sources, FilterMagnitude | None]] = []
    rows: list[dict[str, object]] = []
    for single in targets:
        info = single.source_info
        if info.observation_field.telescope is not result.telescope:
            raise ValueError(
                "Targets and standards must share the same Telescope instance."
            )
        instrumental = single.get_magnitudes(kind="instrumental")
        filters = result.telescope.filters()
        calibrated_values = np.full(len(filters), np.nan)
        for index, filter_name in enumerate(filters):
            m_inst = (
                float(instrumental.magnitudes[index].to_value(u.mag))
                if instrumental is not None
                else np.nan
            )
            error = (
                float(instrumental.errors[index].to_value(u.mag))
                if instrumental is not None
                else np.nan
            )
            significance = (
                float(instrumental.significance[index])
                if instrumental is not None
                else np.nan
            )
            zero_point = float(zero_points[filter_name])
            upper_limit = error == UPPER_LIMIT_ERROR.to_value(u.mag)
            status = "missing_instrumental"
            if np.isfinite(m_inst):
                if not np.isfinite(zero_point):
                    raise ValueError(
                        f"No usable zero point for {info.source_name!r}, {filter_name!r}."
                    )
                calibrated_values[index] = m_inst + zero_point
                status = "upper_limit" if upper_limit else "calibrated"
            rows.append(
                dict(
                    source_id=info.source_name,
                    field=info.observation_field.field_name,
                    filter=filter_name,
                    instrumental_magnitude=m_inst,
                    zero_point=zero_point,
                    calibrated_magnitude=calibrated_values[index],
                    error=error,
                    significance=significance,
                    is_upper_limit=upper_limit,
                    magnitude_system=result.magnitude_system,
                    status=status,
                )
            )
        calibrated = None
        if np.isfinite(calibrated_values).any():
            calibrated = FilterMagnitude(
                telescope=result.telescope,
                filter_names=filters,
                magnitudes=calibrated_values * u.mag,
                errors=instrumental.errors.copy(),
                significance=instrumental.significance.copy(),
                magnitude_system=result.magnitude_system,
            )
        pending.append((single, calibrated))
    table = pd.DataFrame(rows)
    for single, calibrated in pending:
        single.set_magnitudes(kind="calibrated", magnitudes=calibrated)
    return table
