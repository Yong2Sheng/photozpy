"""Regress legacy zero-point arithmetic and shared source calibration results."""

from dataclasses import replace

import astropy.units as u
from astropy.coordinates import SkyCoord
import matplotlib
import numpy as np
import pandas as pd
import pytest
from scipy import stats

from photozpy.sources import FilterMagnitude, ObservationField, SourceInfo, Sources
from photozpy.telescope import SARA_CT, SARA_RM


matplotlib.use("Agg")


def measurement(
    magnitudes,
    *,
    system="instrumental",
    errors=None,
    significance=None,
    filters=("SDSS g", "SDSS r"),
    telescope=SARA_CT,
):
    """Create controlled filter-aligned inputs without image processing."""
    if errors is None:
        errors = np.full(len(magnitudes), 0.03)
    if significance is None:
        significance = np.full(len(magnitudes), 20.0)
    return FilterMagnitude(
        telescope=telescope,
        filter_names=filters,
        magnitudes=magnitudes,
        errors=errors,
        significance=significance,
        magnitude_system=system,
    )


def make_sources(offsets, *, same_field=False, system="AB"):
    """Give every standard an independently owned instrumental/reference pair."""
    coordinate = SkyCoord(150, -5, unit="deg")
    field = ObservationField(
        field_name="standards", telescope=SARA_CT, nominal_pointing=coordinate
    )
    standards = Sources.from_list(
        SourceInfo(
            source_name=f"standard-{index}",
            source_role="standard",
            source_coordinate=coordinate,
            observation_field=field,
        )
        for index in range(len(offsets))
    )
    for single, offset in zip(standards, offsets):
        single.set_magnitudes(kind="instrumental", magnitudes=measurement([-10, -9]))
        single.set_magnitudes(
            kind="reference",
            magnitudes=measurement(
                [offset - 10, offset - 9], system=system, errors=[0.5, 0.8]
            ),
        )
    targets = Sources(
        SourceInfo(
            source_name="target",
            source_role="science",
            source_coordinate=coordinate,
            observation_field=field
            if same_field
            else replace(field, field_name="science"),
        )
    )
    targets.set_magnitudes(kind="instrumental", magnitudes=measurement([-8, -7]))
    return standards, targets


def test_legacy_single_pass_zscore_mean_and_population_scatter():
    offsets = np.array([25 + index / 100 for index in range(11)] + [27.0])
    standards, targets = make_sources(offsets)
    originals = [
        (
            single.get_magnitudes(kind="instrumental"),
            single.get_magnitudes(kind="reference"),
        )
        for single in standards
    ]
    summary = standards.calculate_zero_points()
    keep = np.abs(stats.zscore(offsets)) <= 3
    retained = offsets[keep]
    g = summary.set_index("filter").loc["SDSS_g'"]
    assert g.zero_point == pytest.approx(np.mean(retained))
    assert g.scatter == pytest.approx(np.std(retained))
    assert g.n_candidates == 12
    assert g.n_used == 11
    assert g.n_clipped == 1
    assert g.magnitude_system == "AB"
    details = standards.zero_point_details
    clipped = details[(details["filter"] == "SDSS_g'") & ~details["used"]]
    assert clipped.source_id.tolist() == ["standard-11"]
    assert clipped.reason.tolist() == ["sigma_clipped"]
    # Snapshot tables are independent of the authoritative calculation.
    summary.loc[:, "zero_point"] = 0
    details.loc[:, "used"] = False
    assert standards.zero_points.set_index("filter").loc["SDSS_g'", "zero_point"] > 25
    assert standards.zero_point_details["used"].any()
    for single, (instrumental, reference) in zip(standards, originals):
        assert single.get_magnitudes(kind="instrumental") is instrumental
        assert single.get_magnitudes(kind="reference") is reference
        assert single.get_magnitudes(kind="calibrated") is None
    assert targets.get_magnitudes(kind="calibrated") is None


@pytest.mark.parametrize("offsets", [[25.0], [25.0, 25.0, 25.0]])
def test_singleton_and_constant_samples_remain_usable(offsets):
    standards, _ = make_sources(offsets)
    row = standards.calculate_zero_points().set_index("filter").loc["SDSS_g'"]
    assert row.zero_point == 25
    assert row.scatter == 0
    assert row.n_used == len(offsets)


def test_selection_is_per_filter_and_does_not_use_catalog_errors():
    standards, _ = make_sources([25.0, 25.2, 30, 40, 50, 60])
    standards[1].set_magnitudes(
        kind="instrumental", magnitudes=measurement([-10, -9], significance=[2.9, 3])
    )
    standards[2].set_magnitudes(
        kind="instrumental", magnitudes=measurement([-10, -9], errors=[-99, 0.03])
    )
    standards[3].set_magnitudes(
        kind="instrumental", magnitudes=measurement([np.nan, -9])
    )
    standards[4].set_magnitudes(kind="reference", magnitudes=None)
    standards[5].set_magnitudes(
        kind="instrumental",
        magnitudes=measurement([-10, -9], significance=[np.nan, 20]),
    )
    summary = standards.calculate_zero_points().set_index("filter")
    assert summary.loc["SDSS_g'", "zero_point"] == 25
    assert summary.loc["SDSS_g'", "n_used"] == 1
    assert summary.loc["SDSS_r'", "zero_point"] == pytest.approx(
        np.mean([25, 25.2, 30, 40, 60])
    )
    reasons = standards.zero_point_details.query('filter == "SDSS_g\'"').reason.tolist()
    assert reasons == [
        "retained",
        "low_significance",
        "upper_limit",
        "missing_magnitude",
        "missing_reference",
        "missing_significance",
    ]
    assert np.isnan(summary.loc["SDSS_u'", "zero_point"])


@pytest.mark.parametrize("system", ["AB", "Vega"])
@pytest.mark.parametrize("same_field", [False, True])
def test_unified_calibration_preserves_errors_references_limits_and_shared_views(
    system, same_field
):
    standards, targets = make_sources(
        [25, 25.1, 24.9], same_field=same_field, system=system
    )
    instrumental = measurement([-8, -7], errors=[0.04, -99], significance=[30, 1.5])
    reference = measurement([17, 18], system=system)
    targets.set_magnitudes(kind="instrumental", magnitudes=instrumental)
    targets.set_magnitudes(kind="reference", magnitudes=reference)
    combined = Sources.combine(standards, targets)
    # Default calculation also works without an explicit statistics call.
    table = combined.calibrate_targets(standard_sources=standards)
    calibrated = targets.get_magnitudes(kind="calibrated")
    assert calibrated.magnitude_system == system
    np.testing.assert_allclose(
        calibrated.magnitudes.value, [np.nan, 17, 18, np.nan, np.nan], equal_nan=True
    )
    np.testing.assert_array_equal(calibrated.errors.value, instrumental.errors.value)
    np.testing.assert_array_equal(calibrated.significance, instrumental.significance)
    assert targets.get_magnitudes(kind="reference") is reference
    assert targets.get_magnitudes(kind="instrumental") is instrumental
    assert combined["target"].get_magnitudes(kind="calibrated") is calibrated
    assert all(single.get_magnitudes(kind="calibrated") is None for single in standards)
    assert table.source_id.unique().tolist() == ["target"]
    g = table.set_index("filter").loc["SDSS_g'"]
    r = table.set_index("filter").loc["SDSS_r'"]
    assert g.calibrated_magnitude == 17
    assert g.error == 0.04  # No catalog/scatter/zero-point error contribution.
    assert r.calibrated_magnitude == 18
    assert r.is_upper_limit and r.error == -99 and r.status == "upper_limit"
    repeated = targets.calibrate_targets(standard_sources=standards)
    pd.testing.assert_frame_equal(table, repeated)


def test_custom_clipping_result_is_reused_for_calibration():
    standards, targets = make_sources([25, 25.01, 25.02, 26])
    row = standards.calculate_zero_points(sigma=1).set_index("filter").loc["SDSS_g'"]
    assert row.n_used == 3
    table = targets.calibrate_targets(standard_sources=standards)
    assert table.set_index("filter").loc[
        "SDSS_g'", "calibrated_magnitude"
    ] == pytest.approx(17.01)


@pytest.mark.parametrize(
    "kind,column",
    [
        ("instrumental", "magnitude"),
        ("instrumental", "significance"),
        ("reference", "magnitude"),
    ],
)
def test_mutable_standard_measurements_require_explicit_recalculation(kind, column):
    standards, targets = make_sources([25, 25.1])
    standards.calculate_zero_points()
    original = targets.get_magnitudes(kind="calibrated")
    qtable = standards[0].get_magnitudes(kind=kind).qtable
    index = list(qtable["filter"]).index("SDSS_g'")
    qtable[column][index] += 1 if column == "significance" else 1 * u.mag
    with pytest.raises(RuntimeError, match="measurements changed"):
        targets.calibrate_targets(standard_sources=standards)
    with pytest.raises(RuntimeError, match="measurements changed"):
        standards.zero_points
    assert targets.get_magnitudes(kind="calibrated") is original
    standards.calculate_zero_points()
    targets.calibrate_targets(standard_sources=standards)


def test_missing_zero_point_fails_before_any_target_is_updated():
    standards, first = make_sources([25])
    standards[0].get_magnitudes(kind="reference").qtable["magnitude"][2] = (
        np.nan * u.mag
    )
    second = Sources(replace(first.source_info, source_name="second"))
    first.set_magnitudes(
        kind="instrumental", magnitudes=measurement([-8], filters=["SDSS g"])
    )
    old = measurement([20], system="AB", filters=["SDSS g"])
    first.set_magnitudes(kind="calibrated", magnitudes=old)
    second.set_magnitudes(kind="instrumental", magnitudes=measurement([-8, -7]))
    with pytest.raises(ValueError, match="No usable zero point"):
        Sources.combine(first, second).calibrate_targets(standard_sources=standards)
    assert first.get_magnitudes(kind="calibrated") is old
    assert second.get_magnitudes(kind="calibrated") is None


def test_missing_instrumental_does_not_invent_calibrated_values():
    standards, targets = make_sources([25])
    targets.set_magnitudes(kind="instrumental", magnitudes=None)
    table = targets.calibrate_targets(standard_sources=standards)
    assert table.calibrated_magnitude.isna().all()
    assert (table.status == "missing_instrumental").all()
    assert targets.get_magnitudes(kind="calibrated") is None


def test_invalid_standard_sets_and_telescope_mismatch():
    standards, targets = make_sources([25, 25])
    with pytest.raises(RuntimeError, match="calculate_zero_points"):
        standards.zero_point_details
    with pytest.raises(ValueError, match="standard-role"):
        targets.calculate_zero_points()
    with pytest.raises(ValueError, match="science-role"):
        standards.calibrate_targets(standard_sources=standards)
    with pytest.raises(TypeError, match="Sources"):
        targets.calibrate_targets(standard_sources=[])
    standards[1].set_magnitudes(
        kind="reference", magnitudes=measurement([15, 16], system="Vega")
    )
    with pytest.raises(ValueError, match="common AB or Vega"):
        standards.calculate_zero_points()
    standards[1].set_magnitudes(kind="reference", magnitudes=None)
    standards.calculate_zero_points()
    other = Sources(
        replace(
            targets.source_info,
            observation_field=replace(
                targets.source_info.observation_field, telescope=SARA_RM
            ),
        )
    )
    with pytest.raises(ValueError, match="same Telescope"):
        other.calibrate_targets(standard_sources=standards)


@pytest.mark.parametrize(
    "kwargs",
    [
        dict(sigma=0),
        dict(sigma=np.nan),
        dict(min_significance=-1),
        dict(min_significance=np.inf),
    ],
)
def test_invalid_thresholds(kwargs):
    standards, _ = make_sources([25])
    with pytest.raises(ValueError):
        standards.calculate_zero_points(**kwargs)


def test_all_ineligible_or_all_clipped_are_explicit():
    standards, _ = make_sources([25, 27])
    with pytest.raises(ValueError, match="No usable"):
        standards.calculate_zero_points(sigma=0.1)
    with pytest.raises(ValueError, match="No usable"):
        standards.calculate_zero_points(min_significance=100)


@pytest.mark.parametrize(
    "filter_names", [("SDSS g",), ("SDSS g", "SDSS r"), ("SDSS g", "SDSS r", "SDSS i")]
)
def test_plot_handles_one_two_and_three_filters_without_writing(
    tmp_path, monkeypatch, filter_names
):
    import matplotlib.pyplot as plt

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(plt, "show", lambda: None)
    standards, _ = make_sources([25, 25.1])
    for single in standards:
        single.set_magnitudes(
            kind="instrumental",
            magnitudes=measurement([-10] * len(filter_names), filters=filter_names),
        )
        single.set_magnitudes(
            kind="reference",
            magnitudes=measurement(
                [15] * len(filter_names), filters=filter_names, system="AB"
            ),
        )
    before = standards.calculate_zero_points()
    fig = standards.plot_zero_point_statistics()
    axes = [ax for ax in fig.axes if ax.get_visible()]
    assert len(axes) == len(filter_names)
    for ax in axes:
        np.testing.assert_allclose(
            [line.get_ydata()[0] for line in ax.lines], [24.9, 25.1]
        )
    assert not list(tmp_path.iterdir())
    pd.testing.assert_frame_equal(before, standards.zero_points)
    plt.close(fig)
    saved = tmp_path / "zero_points.png"
    fig = standards.plot_zero_point_statistics(save_path=saved)
    assert saved.stat().st_size > 0
    plt.close(fig)


def test_plot_includes_rejected_extremes_and_distinguishes_missing_offsets(monkeypatch):
    import matplotlib.pyplot as plt

    monkeypatch.setattr(plt, "show", lambda: None)
    offsets = [25 + index / 100 for index in range(11)] + [40, 30, 29]
    standards, _ = make_sources(offsets)
    # Keep a finite low-significance offset visible alongside the clipped outlier.
    standards[-2].set_magnitudes(
        kind="instrumental",
        magnitudes=measurement([-10, -9], significance=[2, 2]),
    )
    # A missing reference has no zero point and must not be plotted as a number.
    standards[-1].set_magnitudes(kind="reference", magnitudes=None)
    summary = standards.calculate_zero_points()
    details = standards.zero_point_details
    assert summary.set_index("filter").loc["SDSS_g'", "n_clipped"] == 1

    fig = standards.plot_zero_point_statistics()
    try:
        for ax in (ax for ax in fig.axes if ax.get_visible()):
            assert "11 retained, 2 rejected, 1 unavailable" in ax.get_title()
            rejected = next(
                collection for collection in ax.collections
                if collection.get_label() == "rejected offsets"
            )
            np.testing.assert_array_equal(rejected.get_offsets()[:, 1], [40, 30])
            lower, upper = ax.get_ylim()
            assert lower < 25 < 40 < upper
        pd.testing.assert_frame_equal(summary, standards.zero_points)
        pd.testing.assert_frame_equal(details, standards.zero_point_details)
    finally:
        plt.close(fig)
