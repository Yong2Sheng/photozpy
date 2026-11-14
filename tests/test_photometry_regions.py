"""Scientific geometry and DS9 edit round trips for per-image region sessions."""

import json
import shutil
from dataclasses import replace
from unittest.mock import Mock

from astropy.io import fits
from astropy.wcs import WCS
from ccdproc import ImageFileCollection
import numpy as np
import pytest
from regions import Regions

from photozpy.sources import ObservationField, SourceInfo, Sources
from photozpy.telescope import SARA_CT
from photozpy.photometry import PhotometryRegions
from photozpy.photometry import regions as module
from photozpy.interactive_ds9.region_reviewer import RegionReviewer
from photozpy.interactive_ds9.adapter import DS9Adapter
from test_interactive_ds9_image_reviewer import fake_tkinter as fake_tkinter


@pytest.fixture(autouse=True)
def isolate_review_debug_logs(tmp_path, monkeypatch):
    """Keep the existing reviewer's diagnostic output in test directories."""
    monkeypatch.chdir(tmp_path)


def inputs(tmp_path, *, mask=False, cmask=False):
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [65, 65]
    wcs.wcs.crval = [163, -5]
    wcs.wcs.cdelt = [-0.0001, 0.0001]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    y, x = np.indices((128, 128))
    data = 100 + 2000 * np.exp(-((x - 61.3) ** 2 + (y - 59.6) ** 2) / (2 * 1.3**2))
    data += 1700 * np.exp(-((x - 92.2) ** 2 + (y - 91.7) ** 2) / (2 * 1.3**2))
    header = wcs.to_header()
    header.update(
        IMTYPE="Master Light",
        FIELD="DLS_field",
        FILTER="SDSS_g'",
        FWHM=3.0,
        BUNIT="adu",
    )
    flags = np.zeros(data.shape, dtype=np.uint8)
    if mask:
        flags[55:65, 61:65] = 1
    path = tmp_path / "Master_DLS_g.fits"
    hdus = [fits.PrimaryHDU(data, header), fits.ImageHDU(flags, name="MASK")]
    if cmask:
        hdus.append(fits.ImageHDU(np.ones(data.shape, dtype=np.uint8), name="CMASK"))
    fits.HDUList(hdus).writeto(path)
    field = ObservationField(
        field_name="DLS_field",
        telescope=SARA_CT,
        nominal_pointing=wcs.pixel_to_world(64, 64),
    )
    sources = Sources(
        *[
            SourceInfo(
                source_name=name,
                source_role="standard",
                observation_field=field,
                source_coordinate=wcs.pixel_to_world(x, y),
            )
            for name, x, y in [("catalog-A", 60, 60), ("catalog-B", 92, 92)]
        ]
    )
    return path, sources


def combined(session):
    return (
        session.aperture_path.read_text()
        + session.background_path.read_text()
        + session.initial_regions()
    )


def test_centroids_radii_pixel_origin_ids_and_image_unchanged(tmp_path):
    path, sources = inputs(tmp_path, cmask=True)
    original = path.read_bytes()
    session = PhotometryRegions.generate(path, sources, root=tmp_path)
    first = session.records[0]
    assert first["x"] == pytest.approx(61.3, abs=0.25)
    assert first["y"] == pytest.approx(59.6, abs=0.25)
    assert first["centroid_method"] in ("quadratic", "com_fallback")
    assert first["generated_radii"] == [9, 15, 24]
    parsed = Regions.read(session.aperture_path, format="ds9")
    assert parsed[0].center.x == pytest.approx(first["x"])
    assert parsed[0].meta["tag"] == ["id:catalog-A", "aperture"]
    session.save_review(combined(session))
    assert session.records[0]["reviewed"]
    assert session.records[0]["initial_x"] == pytest.approx(60)
    assert path.read_bytes() == original
    assert sources[0].get_magnitudes(kind="instrumental") is None


def test_custom_scaling_and_fallback_use_initial_coordinate(tmp_path, monkeypatch):
    path, sources = inputs(tmp_path)
    monkeypatch.setattr(
        module, "centroid_sources", Mock(side_effect=ValueError("fit failed"))
    )
    session = PhotometryRegions.generate(
        path,
        sources,
        root=tmp_path,
        aperture_scaling=2,
        inner_annulus_scaling=4,
        outer_annulus_scaling=7,
    )
    row = session.records[0]
    assert row["centroid_method"] == "initial_fallback"
    assert row["centroid_reason"] == "fit failed"
    assert row["x"] == pytest.approx(60)
    assert row["generated_radii"] == [6, 12, 21]


@pytest.mark.parametrize(
    "options",
    [
        {"aperture_scaling": 5},
        {"outer_annulus_scaling": 4},
        {"inner_annulus_scaling": np.nan},
        {"box_size": 20},
        {"max_centroid_shift": 0},
    ],
)
def test_invalid_controls_do_not_generate_files(tmp_path, options):
    path, sources = inputs(tmp_path)
    with pytest.raises(ValueError):
        PhotometryRegions.generate(path, sources, root=tmp_path, **options)
    assert not (tmp_path / "photometry").exists()


def test_edit_move_resize_and_delete_pairs_by_id_not_order(tmp_path):
    path, sources = inputs(tmp_path)
    session = PhotometryRegions.generate(path, sources, root=tmp_path)
    # Only B remains, with a manual center and non-default radii.
    text = (
        "image\nannulus(94,95,14,22) # tag={background} tag={id:catalog-B}\n"
        "circle(94,95,7) # tag={id:catalog-B} tag={aperture}\n"
    )
    session.save_review(text)
    assert not session.records[0]["accepted"]
    second = session.records[1]
    assert (second["x"], second["y"], second["aperture_radius"]) == (93, 94, 7)
    assert (second["inner_radius"], second["outer_radius"]) == (14, 22)
    assert len(Regions.read(session.aperture_path, format="ds9")) == 1
    assert "point(" not in session.aperture_path.read_text()
    with pytest.raises(FileExistsError):
        PhotometryRegions.generate(path, sources, root=tmp_path)


@pytest.mark.parametrize(
    "text",
    [
        "image\ncircle(61,61,9) # tag={id:catalog-A} tag={aperture}\n",
        "image\ncircle(61,61,9) # tag={id:unknown} tag={aperture}\n",
        "image\ncircle(61,61,9)\n",
    ],
)
def test_invalid_review_does_not_overwrite_files(tmp_path, text):
    path, sources = inputs(tmp_path)
    session = PhotometryRegions.generate(path, sources, root=tmp_path)
    files = (
        session.aperture_path,
        session.background_path,
        tmp_path / session.manifest_path,
    )
    originals = [file.read_bytes() for file in files]
    with pytest.raises(ValueError):
        session.save_review(text)
    assert [file.read_bytes() for file in files] == originals


def test_resume_after_moving_root_and_record_copies(tmp_path):
    original_root = tmp_path / "original"
    original_root.mkdir()
    path, sources = inputs(original_root)
    session = PhotometryRegions.generate(path, sources, root=original_root)
    state = json.loads((original_root / session.manifest_path).read_text())
    assert not any(str(original_root) in str(value) for value in state.values())
    moved = tmp_path / "moved"
    shutil.move(original_root, moved)
    resumed = PhotometryRegions.load(moved, session.manifest_path)
    resumed.save_review(combined(resumed))
    assert resumed.image_path == moved / path.name
    records = resumed.records
    records[0]["x"] = -100
    assert resumed.records[0]["x"] != -100


def test_mask_excludes_values_but_cmask_does_not(tmp_path):
    path, sources = inputs(tmp_path)
    session = PhotometryRegions.generate(path, sources, root=tmp_path)
    original_x = session.records[0]["x"]
    with fits.open(path, mode="update") as hdul:
        hdul[0].data[58, 60] = 1e10
        hdul["MASK"].data[58, 60] = 1
        hdul.append(
            fits.ImageHDU(np.ones(hdul[0].data.shape, dtype=np.uint8), name="CMASK")
        )
    masked = PhotometryRegions.generate(path, sources, root=tmp_path, overwrite=True)
    assert masked.records[0]["x"] == pytest.approx(original_x, abs=0.4)
    assert masked.records[0]["aperture_valid_fraction"] < 1


@pytest.mark.usefixtures("fake_tkinter")
def test_reviewer_loads_overlays_and_only_advances_after_valid_save(tmp_path):
    path, sources = inputs(tmp_path)
    session = PhotometryRegions.generate(path, sources, root=tmp_path)
    adapter = Mock(spec=DS9Adapter)
    reviewer = RegionReviewer([session], adapter=adapter, isolated=False)
    reviewer.show_current()
    adapter.load_region.assert_any_call(session.aperture_path)
    assert "text=" not in adapter.add_regions.call_args_list[0].args[0]
    adapter.get.return_value = (
        "image\ncircle(61,61,9) # tag={id:catalog-A} tag={aperture}\n"
    )
    reviewer.mark_good()
    assert reviewer.index == 0
    adapter.get.return_value = combined(session)
    reviewer.mark_good()
    assert reviewer.index == 1
    assert all(row["reviewed"] for row in session.records)
    reviewer.go_back()
    assert reviewer.index == 0
    reviewer.quit_review()


def test_offset_centers_saved_independently_with_warning_and_single_label(tmp_path):
    path, sources = inputs(tmp_path)
    session = PhotometryRegions.generate(path, sources, root=tmp_path)
    text = (
        "image\ncircle(61,61,9) # tag={id:catalog-A} tag={aperture}\n"
        "annulus(64,65,15,24) # tag={id:catalog-A} tag={background}\n"
    )
    messages = session.save_review(text)
    row = session.records[0]
    assert (row["x"], row["y"]) == (60, 60)
    assert (row["background_x"], row["background_y"]) == (63, 64)
    assert row["center_offset"] == 5
    assert len(messages) == 1 and "5.00 pixels" in messages[0]
    aperture = Regions.read(session.aperture_path, format="ds9")[0]
    background = Regions.read(session.background_path, format="ds9")[0]
    assert (background.center.x, background.center.y) == (63, 64)
    assert aperture.meta["text"] == "catalog-A"
    assert "text" not in background.meta
    assert background.meta["tag"] == ["id:catalog-A", "background"]
    # Subsequent saves must not snap the annulus back to the aperture center.
    session.save_review(combined(session))
    assert session.records[0]["background_x"] == 63


def test_legacy_background_labels_are_hidden_without_modifying_file(tmp_path):
    path, sources = inputs(tmp_path)
    session = PhotometryRegions.generate(path, sources, root=tmp_path)
    text = session.background_path.read_text().replace(
        "tag={background}", "tag={background} text={old label}"
    )
    session.background_path.write_text(text)
    display = session.background_display_regions()
    assert "old label" not in display
    assert "id:catalog-A" in display
    assert session.background_path.read_text() == text


@pytest.mark.usefixtures("fake_tkinter")
def test_offset_warning_is_red_and_persists_after_next(tmp_path):
    path, sources = inputs(tmp_path)
    session = PhotometryRegions.generate(path, sources, root=tmp_path)
    adapter = Mock(spec=DS9Adapter)
    reviewer = RegionReviewer([session], adapter=adapter, isolated=False)
    reviewer.show_current()
    adapter.get.return_value = (
        "image\ncircle(61,61,9) # tag={id:catalog-A} tag={aperture}\n"
        "annulus(64,65,15,24) # tag={id:catalog-A} tag={background}\n"
    )
    reviewer.mark_good()
    assert reviewer.index == 1
    assert "Warning" in reviewer.error_label.options["text"]
    assert "5.00 pixels" in reviewer.error_label.options["text"]
    assert reviewer.error_label.options["fg"] == "firebrick"
    assert reviewer.error_message == ""
    reviewer.quit_review()


def test_region_review_uses_existing_isolated_worker(tmp_path, monkeypatch):
    from photozpy.interactive_ds9 import _review_process

    path, sources = inputs(tmp_path)
    session = PhotometryRegions.generate(path, sources, root=tmp_path)
    worker = Mock()
    monkeypatch.setattr(_review_process, "run_isolated", worker)
    reviewer = RegionReviewer([session], isolated=True)
    reviewer.run()
    assert worker.call_args.kwargs["region_sessions"] == [
        {"root": str(tmp_path), "manifest": session.manifest_path.as_posix()},
    ]


def test_fwhm_is_read_from_each_image_and_large_centroid_shift_falls_back(
    tmp_path, monkeypatch
):
    path, sources = inputs(tmp_path)
    first = PhotometryRegions.generate(path, sources, root=tmp_path)
    second_path = tmp_path / "Master_DLS_r.fits"
    with fits.open(path) as hdul:
        hdul[0].header["FWHM"] = 4.0
        hdul[0].header["FILTER"] = "SDSS_r'"
        hdul.writeto(second_path)
    monkeypatch.setattr(
        module,
        "centroid_sources",
        lambda *args, **kwargs: (np.array([120.0]), np.array([120.0])),
    )
    second = PhotometryRegions.generate(second_path, sources, root=tmp_path)
    assert first.records[0]["generated_radii"] == [9, 15, 24]
    assert second.records[0]["generated_radii"] == [12, 20, 32]
    assert second.records[0]["x"] == pytest.approx(60)
    assert second.records[0]["centroid_method"] == "initial_fallback"


def test_background_failure_is_recorded_and_all_deleted_regions_are_saved(
    tmp_path, monkeypatch
):
    path, sources = inputs(tmp_path)
    monkeypatch.setattr(
        module, "Background2D", Mock(side_effect=ValueError("insufficient sky"))
    )
    session = PhotometryRegions.generate(path, sources, root=tmp_path)
    assert all(row["centroid_method"] == "initial_fallback" for row in session.records)
    assert "insufficient sky" in session.records[0]["centroid_reason"]
    session.save_review(session.initial_regions())
    assert not any(row["accepted"] for row in session.records)
    assert "circle(" not in session.aperture_path.read_text()
    assert not Regions.read(session.background_path, format="ds9")


def test_worker_dispatches_region_sessions_without_notebook_tk(tmp_path, monkeypatch):
    from photozpy.interactive_ds9 import _review_process, adapter, region_reviewer

    path, sources = inputs(tmp_path)
    session = PhotometryRegions.generate(path, sources, root=tmp_path)
    request, state = tmp_path / "request.json", tmp_path / "state.json"
    request.write_text(
        json.dumps(
            {
                "region_sessions": [
                    {
                        "root": str(tmp_path),
                        "manifest": session.manifest_path.as_posix(),
                    }
                ],
                "adapter": None,
            }
        )
    )
    constructor = Mock()
    monkeypatch.setattr(adapter, "DS9Adapter", Mock())
    monkeypatch.setattr(region_reviewer, "RegionReviewer", constructor)
    monkeypatch.setattr(
        _review_process.sys, "argv", ["worker", str(request), str(state)]
    )
    _review_process.main()
    assert constructor.call_args.kwargs["isolated"] is False
    assert constructor.call_args.args[0][0].image_path == path
    constructor.return_value.run.assert_called_once()


def second_filter(path, filter_name="SDSS_r'", field_name="DLS_field"):
    destination = path.parent / f"Master_{field_name}_{filter_name}.fits"
    with fits.open(path) as hdul:
        hdul[0].header["FILTER"] = filter_name
        hdul[0].header["FIELD"] = field_name
        hdul[0].header["FWHM"] = 4.0
        hdul.writeto(destination)
    return destination


def test_batch_standards_share_directory_and_preserve_image_filter_metadata(tmp_path):
    path, sources = inputs(tmp_path)
    red = second_filter(path)
    originals = {image: image.read_bytes() for image in (path, red)}
    collection = ImageFileCollection(
        tmp_path, filenames=[red.name, path.name], keywords="*"
    )
    batch = PhotometryRegions.generate_regions(collection, sources, root=tmp_path)
    assert len(batch) == 2
    green, red_session = batch.sessions
    assert green.aperture_path.parent == red_session.aperture_path.parent
    assert green.aperture_path.parent.name == "DLS_field_standards"
    assert green.aperture_path.name == "SDSS_g'_aperture.reg"
    assert red_session.background_path.name == "SDSS_r'_background.reg"
    assert green.records[0]["generated_radii"] == [9, 15, 24]
    assert red_session.records[0]["generated_radii"] == [12, 20, 32]
    assert len(batch.records) == 4
    for session in batch.sessions:
        for file in (session.aperture_path, session.background_path):
            text = file.read_text()
            assert '# FIELD_NAME="DLS_field"' in text
            assert "# FILTER_NAME=" in text
            assert session.image_path.name in text
            assert "catalog-A" in text and "catalog-B" in text
            assert str(tmp_path) not in text
    assert {image: image.read_bytes() for image in originals} == originals


def test_batch_science_targets_get_individual_files_and_combined_records(tmp_path):
    path, standards = inputs(tmp_path)
    second_filter(path)
    targets = Sources(
        *[
            replace(info, source_name=f"target {index}", source_role="science")
            for index, info in enumerate(standards.source_infos)
        ]
    )
    batch = PhotometryRegions.generate_regions(
        ImageFileCollection(tmp_path), targets, root=tmp_path
    )
    assert len(batch) == 4
    assert len(batch.records) == 4
    assert [row["source_id"] for row in batch.records] == [
        "target 0",
        "target 0",
        "target 1",
        "target 1",
    ]
    assert [row["filter_name"] for row in batch.records] == [
        "SDSS_g'",
        "SDSS_r'",
        "SDSS_g'",
        "SDSS_r'",
    ]
    for session in batch.sessions:
        assert len(session.records) == 1
        name = session.records[0]["source_id"]
        assert session.aperture_path.parent.name == name.replace(" ", "_")
        assert session.aperture_path.name.startswith("DLS_field_SDSS_")
        assert session._state["centroid_controls"]["box_size"] == 25
    assert all(target.get_magnitudes(kind="instrumental") is None for target in targets)


def test_batch_preflight_rejects_duplicate_outputs_and_missing_fields_before_writes(
    tmp_path,
):
    path, sources = inputs(tmp_path)
    duplicate = tmp_path / "duplicate.fits"
    shutil.copy2(path, duplicate)
    with pytest.raises(ValueError, match="Multiple images"):
        PhotometryRegions.generate_regions(
            ImageFileCollection(tmp_path), sources, root=tmp_path
        )
    assert not (tmp_path / "photometry").exists()
    missing_info = replace(
        sources.source_infos[0],
        source_name="missing",
        observation_field=replace(
            sources.source_infos[0].observation_field, field_name="absent"
        ),
    )
    with pytest.raises(ValueError, match="No Master Light"):
        PhotometryRegions.generate_regions(
            ImageFileCollection(tmp_path, filenames=[path.name]),
            Sources(*sources.source_infos, missing_info),
            root=tmp_path,
        )
    assert not (tmp_path / "photometry").exists()


def test_science_same_name_in_different_fields_has_distinct_filter_files(tmp_path):
    path, standards = inputs(tmp_path)
    second_filter(path, filter_name="SDSS_g'", field_name="other_field")
    first = replace(
        standards.source_infos[0], source_name="target", source_role="science"
    )
    second = replace(
        first,
        observation_field=replace(first.observation_field, field_name="other_field"),
    )
    batch = PhotometryRegions.generate_regions(
        ImageFileCollection(tmp_path),
        Sources(first, second),
        root=tmp_path,
    )
    assert len(batch) == 2
    assert (
        batch.sessions[0].aperture_path.parent == batch.sessions[1].aperture_path.parent
    )
    assert batch.sessions[0].aperture_path != batch.sessions[1].aperture_path


def test_batch_review_uses_one_queue_and_records_are_copies(tmp_path, monkeypatch):
    from photozpy.interactive_ds9 import region_reviewer

    path, sources = inputs(tmp_path)
    second_filter(path)
    batch = PhotometryRegions.generate_regions(
        ImageFileCollection(tmp_path), sources, root=tmp_path
    )
    constructor = Mock()
    monkeypatch.setattr(region_reviewer, "RegionReviewer", constructor)
    batch.review(isolated=True)
    constructor.assert_called_once_with(batch.sessions, isolated=True)
    constructor.return_value.run.assert_called_once()
    records = batch.records
    records[0]["x"] = -100
    assert batch.records[0]["x"] != -100


@pytest.mark.usefixtures("fake_tkinter")
def test_region_navigation_preserves_view_within_target_and_resets_between_targets(
    tmp_path,
):
    path, standards = inputs(tmp_path)
    second_filter(path)
    targets = Sources(
        *[
            replace(info, source_name=f"target {index}", source_role="science")
            for index, info in enumerate(standards.source_infos)
        ]
    )
    batch = PhotometryRegions.generate_regions(
        ImageFileCollection(tmp_path), targets, root=tmp_path
    )
    adapter = Mock(spec=DS9Adapter)
    reviewer = RegionReviewer(batch.sessions, adapter=adapter, isolated=False)
    reviewer.show_current()
    assert adapter.load_fits.call_args.kwargs == {"preserve_view": False}
    adapter.get.return_value = combined(batch.sessions[0])
    reviewer.mark_good()
    assert adapter.load_fits.call_args.args[0] == batch.sessions[1].image_path
    assert adapter.load_fits.call_args.kwargs == {"preserve_view": True}
    reviewer.go_back()
    assert reviewer.index == 0
    assert adapter.load_fits.call_args.kwargs == {"preserve_view": True}
    reviewer.mark_good()
    adapter.get.return_value = combined(batch.sessions[1])
    reviewer.mark_good()
    assert reviewer.index == 2
    assert adapter.load_fits.call_args.kwargs == {"preserve_view": False}
    # A second target's next filter inherits its view, including at completion.
    adapter.get.return_value = combined(batch.sessions[2])
    reviewer.mark_good()
    assert adapter.load_fits.call_args.kwargs == {"preserve_view": True}
    adapter.get.return_value = combined(batch.sessions[3])
    reviewer.mark_good()
    reviewer.go_back()
    assert reviewer.index == 3
    assert adapter.load_fits.call_args.kwargs == {"preserve_view": True}
    reviewer.quit_review()


def test_load_regions_preserves_manual_geometry_without_generation_or_writes(
    tmp_path, monkeypatch
):
    path, sources = inputs(tmp_path)
    second_filter(path)
    collection = ImageFileCollection(tmp_path)
    batch = PhotometryRegions.generate_regions(collection, sources, root=tmp_path)
    batch.sessions[0].save_review(
        "image\ncircle(61,61,7) # tag={id:catalog-A} tag={aperture}\n"
        "annulus(64,65,14,22) # tag={id:catalog-A} tag={background}\n"
    )
    files = list((tmp_path / "photometry").rglob("*"))
    snapshots = {file: file.read_bytes() for file in files if file.is_file()}
    monkeypatch.setattr(
        PhotometryRegions,
        "generate",
        Mock(side_effect=AssertionError("Must not generate")),
    )
    monkeypatch.setattr(
        module,
        "centroid_sources",
        Mock(side_effect=AssertionError("Must not centroid")),
    )
    restored = PhotometryRegions.load_regions(collection, sources, root=tmp_path)
    assert len(restored) == 2
    assert restored.sessions[0].records == batch.sessions[0].records
    row = restored.sessions[0].records[0]
    assert row["aperture_radius"] == 7
    assert row["background_x"] == 63 and row["background_y"] == 64
    assert not restored.sessions[0].records[1]["accepted"]
    assert {file: file.read_bytes() for file in snapshots} == snapshots


def test_load_regions_reports_missing_files_without_recreating_them(tmp_path):
    path, sources = inputs(tmp_path)
    collection = ImageFileCollection(tmp_path)
    with pytest.raises(FileNotFoundError, match="nothing was regenerated"):
        PhotometryRegions.load_regions(collection, sources, root=tmp_path)
    assert not (tmp_path / "photometry").exists()


def test_load_regions_rejects_changed_source_identities(tmp_path):
    path, sources = inputs(tmp_path)
    PhotometryRegions.generate_regions(
        ImageFileCollection(tmp_path), sources, root=tmp_path
    )
    different = Sources(replace(sources.source_infos[0], source_name="different-star"))
    with pytest.raises(ValueError, match="source identities"):
        PhotometryRegions.load_regions(
            ImageFileCollection(tmp_path), different, root=tmp_path
        )
