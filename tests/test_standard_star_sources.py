"""Validate catalog-to-source association without external catalogs or writes."""

from types import SimpleNamespace

from astropy.coordinates import SkyCoord
from astropy.io import fits
import numpy as np
import pandas as pd
import pytest

from photozpy.sources import ObservationField
from photozpy.standard_stars import StandardStarCatalog
from photozpy.telescope import SARA_CT


@pytest.fixture
def catalog():
    instance = StandardStarCatalog.__new__(StandardStarCatalog)
    instance.metadata = SimpleNamespace(ra_unit="deg", dec_unit="deg")
    instance._matched_stars = None
    return instance


@pytest.fixture
def field():
    return ObservationField(
        field_name="DLS_field", telescope=SARA_CT,
        nominal_pointing=SkyCoord(163, -5, unit="deg"),
    )


@pytest.fixture
def rows():
    return pd.DataFrame({
        "id_name": [b"star-B", b"star-A"],
        "ra": [163.01, 163.02], "dec": [-5.01, -5.02],
        "g": [14.2, 16.3], "r": [13.8, np.nan],
        "g_err": [0.02, 0.03],
    }, index=[9, 4])


def convert(catalog, rows, field, **options):
    arguments = dict(
        observation_field=field,
        magnitude_columns={"SDSS r": "r", "SDSS g": "g"},
        error_columns={"SDSS_g'": "g_err"},
        magnitude_system="AB",
    )
    arguments.update(options)
    return catalog.to_sources(rows, **arguments)


def test_catalog_values_follow_identity_and_canonical_filter_order(catalog, rows, field):
    original = rows.copy(deep=True)
    sources = convert(catalog, rows, field)
    assert sources.source_names == ("star-B", "star-A")
    for index, single in enumerate(sources):
        info = single.source_info
        assert info.observation_field is field
        assert info.source_role == "standard"
        assert info.source_coordinate.icrs.ra.deg == pytest.approx(rows.iloc[index].ra)
        reference = single.get_magnitudes(kind="reference")
        assert reference.telescope is field.telescope
        assert reference.magnitude_system == "AB"
        assert reference.magnitudes.value[1] == rows.iloc[index].g
        assert reference.errors.value[1] == rows.iloc[index].g_err
        assert np.isnan(reference.errors.value[2])
        assert np.all(np.isnan(reference.significance))
        assert single.get_magnitudes(kind="instrumental") is None
        assert single.get_magnitudes(kind="calibrated") is None
    assert np.isnan(sources[1].get_magnitudes(kind="reference").magnitudes.value[2])
    pd.testing.assert_frame_equal(rows, original)


@pytest.mark.parametrize("options,exception", [
    ({"magnitude_system": "instrumental"}, ValueError),
    ({"magnitude_columns": {}}, ValueError),
    ({"magnitude_columns": {"SDSS g": "g", "SDSS_g'": "r"}}, ValueError),
    ({"error_columns": {"SDSS z": "g_err"}}, ValueError),
    ({"error_columns": {"SDSS g": "absent"}}, KeyError),
])
def test_invalid_column_contracts_fail(catalog, rows, field, options, exception):
    with pytest.raises(exception):
        convert(catalog, rows, field, **options)


def test_duplicate_identity_and_invalid_coordinates_fail(catalog, rows, field):
    rows.loc[4, "id_name"] = b"star-B"
    with pytest.raises(ValueError, match="Duplicate source"):
        convert(catalog, rows, field)
    rows.loc[4, "id_name"] = b"star-A"
    rows.loc[4, "ra"] = np.nan
    with pytest.raises(ValueError, match="Non-finite"):
        convert(catalog, rows, field)


def test_empty_catalog_selection_and_missing_errors(catalog, rows, field):
    assert len(convert(catalog, rows.iloc[:0], field)) == 0
    sources = convert(catalog, rows, field, error_columns=None)
    assert np.all(np.isnan(sources[0].get_magnitudes(kind="reference").errors))


def test_query_forwards_old_cuts_and_rejects_wrong_field(catalog, rows, field, monkeypatch):
    captured = {}

    def query(**kwargs):
        captured.update(kwargs)
        return rows

    monkeypatch.setattr(catalog, "find_standard_stars", query)
    monkeypatch.setattr(fits, "getheader", lambda *args, **kwargs: {"FIELD": "DLS_field"})
    sources = catalog.find_standard_sources(
        "standard.fits", observation_field=field,
        magnitude_columns={"SDSS g": "g"}, error_columns={"SDSS g": "g_err"},
        magnitude_system="AB", edge_margin_pix=100, mag_limit=23,
    )
    assert len(sources) == 2
    assert catalog.matched_stars.id_name.tolist() == list(sources.source_names)
    assert catalog.matched_stars.index.tolist() == rows.index.tolist()
    assert captured == dict(
        fits_path="standard.fits", hdu_index=0, edge_margin_pix=100,
        add_pixel_columns=True, mag_limit=23, mag_columns=["g"],
    )
    captured.clear()
    monkeypatch.setattr(fits, "getheader", lambda *args, **kwargs: {"FIELD": "other"})
    with pytest.raises(ValueError, match="FIELD"):
        catalog.find_standard_sources(
            "other.fits", observation_field=field,
            magnitude_columns={"SDSS g": "g"}, magnitude_system="AB",
        )
    assert not captured


def test_automatic_columns_use_telescope_aliases_and_errors(catalog, rows, field):
    rows = rows.rename(columns={"g": "SDSS_gp", "r": "SDSS_rp", "g_err": "SDSS_gp_err"})
    sources = catalog.to_sources(rows, observation_field=field, magnitude_system="AB")
    reference = sources[0].get_magnitudes(kind="reference")
    assert reference.magnitudes.value[1] == 14.2
    assert reference.errors.value[1] == 0.02
    assert reference.magnitudes.value[2] == 13.8
    assert np.isnan(reference.errors.value[2])


def test_automatic_columns_reject_ambiguity_and_do_not_substitute_passbands(field):
    resolve = StandardStarCatalog.resolve_magnitude_columns
    with pytest.raises(ValueError, match="Multiple catalog columns"):
        resolve(("SDSS_gp", "SDSS_g'"), field.telescope)
    with pytest.raises(ValueError, match="No catalog"):
        resolve(("PanSTARRS_zs", "PanSTARRS_Y"), field.telescope)
    magnitudes, errors = resolve(("SDSS_gp", "PanSTARRS_zs"), field.telescope)
    assert magnitudes == {"SDSS_g'": "SDSS_gp"}
    assert errors == {}


def test_automatic_query_uses_metadata_columns(catalog, field, monkeypatch):
    catalog.metadata.table_columns = ("id_name", "ra", "dec", "SDSS_gp", "SDSS_gp_err")
    rows = pd.DataFrame({"id_name": [b"star"], "ra": [163.0], "dec": [-5.0],
                         "SDSS_gp": [14.2], "SDSS_gp_err": [0.02]})
    captured = {}

    def query(**kwargs):
        captured.update(kwargs)
        return rows

    monkeypatch.setattr(catalog, "find_standard_stars", query)
    monkeypatch.setattr(fits, "getheader", lambda *args, **kwargs: {"FIELD": "DLS_field"})
    sources = catalog.find_standard_sources(
        "standard.fits", observation_field=field, magnitude_system="AB", mag_limit=23,
    )
    assert len(sources) == 1
    assert captured["mag_columns"] == ["SDSS_gp"]
    assert sources[0].get_magnitudes(kind="reference").errors.value[1] == 0.02


def test_ds9_preview_uses_catalog_crosses_without_mutating_sources(catalog, rows, field, monkeypatch):
    sources = convert(catalog, rows, field)
    calls = []

    class Viewer:
        def ensure_running(self):
            calls.append("ready")

        def load_fits(self, path, **kwargs):
            calls.append((path, kwargs))

        def verify_image_loaded(self):
            calls.append("verified")

        def clear_regions(self):
            calls.append("clear")

        def add_regions(self, text):
            calls.append(text)

    monkeypatch.setattr(fits, "getheader", lambda *args, **kwargs: {"FIELD": "DLS_field"})
    viewer = Viewer()
    assert catalog.preview_standard_stars("standard.fits", sources, adapter=viewer) is viewer
    assert calls[:4] == ["ready", ("standard.fits", {"hdu": 0}), "verified", "clear"]
    assert "icrs\npoint(163.0100000000,-5.0100000000)" in calls[-1]
    assert calls[-1].count("point=cross") == 2
    assert "text={star-B}" in calls[-1]
    assert sources.source_names == ("star-B", "star-A")


def test_matched_table_requires_a_completed_query(catalog):
    with pytest.raises(RuntimeError, match="run a catalog query"):
        _ = catalog.matched_stars


@pytest.mark.parametrize("mag_limit,expected_names", [
    (None, ["star-B", "star-A"]),
    (15.0, ["star-B"]),
    (10.0, []),
])
def test_matched_table_applies_all_cuts_and_preserves_legacy_result(
    catalog, rows, monkeypatch, mag_limit, expected_names,
):
    monkeypatch.setattr(catalog, "_filter_star_in_fov", lambda **kwargs: rows)
    monkeypatch.setattr(catalog, "_filter_stars_near_edge", lambda **kwargs: rows)
    result = catalog.find_standard_stars(
        "standard.fits", mag_limit=mag_limit, mag_columns=["g"],
    )
    table = catalog.matched_stars
    assert table.id_name.tolist() == expected_names
    assert table.columns.tolist() == result.columns.tolist()
    assert table.index.tolist() == result.index.tolist()
    np.testing.assert_array_equal(table.g.to_numpy(), result.g.to_numpy())
    # The existing return value keeps original byte IDs and pre-limit state.
    assert all(isinstance(name, bytes) for name in result.id_name)
    assert len(catalog.final_stars) == 2
    # Neither the return value nor the display table can mutate the snapshot.
    if len(result):
        result.iloc[0, result.columns.get_loc("g")] = -100.0
        table.iloc[0, table.columns.get_loc("g")] = -200.0
        assert catalog.matched_stars.g.iloc[0] == 14.2


def test_failed_query_preserves_previous_table(catalog, rows, monkeypatch):
    monkeypatch.setattr(catalog, "_filter_star_in_fov", lambda **kwargs: rows)
    monkeypatch.setattr(catalog, "_filter_stars_near_edge", lambda **kwargs: rows)
    catalog.find_standard_stars("standard.fits")
    previous = catalog.matched_stars
    with pytest.raises(ValueError, match="magnitude names"):
        catalog.find_standard_stars("standard.fits", mag_limit=15)
    pd.testing.assert_frame_equal(catalog.matched_stars, previous)
