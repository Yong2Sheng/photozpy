"""
Unit tests for photozpy.collections.image_collection_helpers.

These tests generate minimal FITS files at runtime instead of storing binary
FITS fixtures in the repository. The FITS data arrays are intentionally tiny:
the helpers tested here operate on ImageFileCollection metadata, file
selection, and FITS headers rather than image science data.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits
from ccdproc import ImageFileCollection

from photozpy.collections.image_collection_group import ImageCollectionGroup
from photozpy.collections.image_collection_helpers import (
    NoMatchingImagesError,
    filter_image_collection,
    filter_image_collection_group,
    get_group_header_values,
    get_header_values,
    refresh_image_collection,
    refresh_image_collection_group,
)


@pytest.fixture
def write_test_fits() -> Callable[..., Path]:
    """Return a helper that writes a minimal FITS image with selected headers."""

    def _write_test_fits(
        path: Path,
        *,
        imagetyp: str,
        filter_name: str | None = None,
        object_name: str | None = None,
        exptime: float | None = None,
    ) -> Path:
        header = fits.Header()
        header["IMAGETYP"] = imagetyp

        if filter_name is not None:
            header["FILTER"] = filter_name

        if object_name is not None:
            header["OBJECT"] = object_name

        if exptime is not None:
            header["EXPTIME"] = exptime

        fits.PrimaryHDU(
            data=np.zeros((2, 2), dtype=np.float32),
            header=header,
        ).writeto(path)

        return path

    return _write_test_fits


@pytest.fixture
def image_collection_directory(
    tmp_path: Path,
    write_test_fits: Callable[..., Path],
) -> Path:
    """Create one representative directory of calibration and science FITS."""
    directory = tmp_path / "night_1"
    directory.mkdir()

    write_test_fits(
        directory / "bias_001.fits",
        imagetyp="BIAS",
        exptime=0.0,
    )
    write_test_fits(
        directory / "flat_b_001.fits",
        imagetyp="FLAT",
        filter_name="B",
        exptime=5.0,
    )
    write_test_fits(
        directory / "light_b_001.fits",
        imagetyp="LIGHT",
        filter_name="B",
        object_name="TARGET_A",
        exptime=30.0,
    )
    write_test_fits(
        directory / "light_b_002.fits",
        imagetyp="LIGHT",
        filter_name="B",
        object_name="TARGET_A",
        exptime=30.0,
    )
    write_test_fits(
        directory / "light_no_filter_001.fits",
        imagetyp="LIGHT",
        object_name="TARGET_A",
        exptime=30.0,
    )
    write_test_fits(
        directory / "light_v_001.fits",
        imagetyp="LIGHT",
        filter_name="V",
        object_name="TARGET_A",
        exptime=60.0,
    )

    return directory


@pytest.fixture
def image_collection(
    image_collection_directory: Path,
) -> ImageFileCollection:
    """Create one directory-backed collection for helper tests."""
    return ImageFileCollection(
        location=image_collection_directory,
        keywords=["IMAGETYP", "FILTER", "OBJECT", "EXPTIME"],
        glob_include="*.fits",
    )


@pytest.fixture
def image_collection_group(
    tmp_path: Path,
    write_test_fits: Callable[..., Path],
) -> ImageCollectionGroup:
    """Create a two-directory group with distinct FITS-header contents."""
    night_1 = tmp_path / "group_night_1"
    night_2 = tmp_path / "group_night_2"
    night_1.mkdir()
    night_2.mkdir()

    write_test_fits(
        night_1 / "bias_001.fits",
        imagetyp="BIAS",
        exptime=0.0,
    )
    write_test_fits(
        night_1 / "light_b_001.fits",
        imagetyp="LIGHT",
        filter_name="B",
        object_name="TARGET_A",
        exptime=30.0,
    )

    write_test_fits(
        night_2 / "flat_v_001.fits",
        imagetyp="FLAT",
        filter_name="V",
        exptime=5.0,
    )
    write_test_fits(
        night_2 / "light_v_001.fits",
        imagetyp="LIGHT",
        filter_name="V",
        object_name="TARGET_B",
        exptime=60.0,
    )

    return ImageCollectionGroup(
        directories=[night_1, night_2],
        keywords=["IMAGETYP", "FILTER", "OBJECT", "EXPTIME"],
        glob_include="*.fits",
    )


def test_refresh_image_collection_discovers_new_fits_file(
    image_collection: ImageFileCollection,
    image_collection_directory: Path,
    write_test_fits: Callable[..., Path],
) -> None:
    """Refreshing a primary directory-backed collection discovers new FITS."""
    original_files = set(image_collection.files)

    write_test_fits(
        image_collection_directory / "master_bias.fits",
        imagetyp="MASTER BIAS",
        exptime=0.0,
    )

    returned_collection = refresh_image_collection(image_collection)

    assert returned_collection is image_collection
    assert set(image_collection.files) == original_files | {"master_bias.fits"}


def test_refresh_filtered_snapshot_does_not_discover_new_matching_file(
    image_collection: ImageFileCollection,
    image_collection_directory: Path,
    write_test_fits: Callable[..., Path],
) -> None:
    """Refreshing a filtered snapshot remains limited to explicit filenames."""
    selected = filter_image_collection(
        image_collection,
        imagetyp="LIGHT",
        filter="B",
    )

    original_snapshot_files = list(selected.files)

    write_test_fits(
        image_collection_directory / "light_b_later.fits",
        imagetyp="LIGHT",
        filter_name="B",
        object_name="TARGET_A",
        exptime=30.0,
    )

    returned_collection = refresh_image_collection(selected)

    assert returned_collection is selected
    assert selected.files == original_snapshot_files
    assert "light_b_later.fits" not in selected.files


def test_filter_image_collection_uses_or_within_one_header(
    image_collection: ImageFileCollection,
) -> None:
    """Multiple accepted values for one header are combined with OR."""
    selected = filter_image_collection(
        image_collection,
        imagetyp=["BIAS", "FLAT"],
    )

    expected_names = {"bias_001.fits", "flat_b_001.fits"}
    expected_files = [
        filename
        for filename in image_collection.files
        if filename in expected_names
    ]

    assert selected.files == expected_files
    assert selected.location == image_collection.location
    assert selected.keywords == image_collection.keywords
    assert selected.ext == image_collection.ext


def test_filter_image_collection_uses_and_between_headers(
    image_collection: ImageFileCollection,
) -> None:
    """Criteria for different headers are combined with AND."""
    selected = filter_image_collection(
        image_collection,
        imagetyp=["LIGHT", "FLAT"],
        filter="B",
    )

    expected_names = {
        "flat_b_001.fits",
        "light_b_001.fits",
        "light_b_002.fits",
    }
    expected_files = [
        filename
        for filename in image_collection.files
        if filename in expected_names
    ]

    assert selected.files == expected_files


def test_filter_image_collection_preserves_source_file_order(
    image_collection: ImageFileCollection,
) -> None:
    """Selection order follows the source collection, not criterion order."""
    selected = filter_image_collection(
        image_collection,
        imagetyp=["LIGHT", "BIAS"],
        filter=["V", "B"],
    )

    expected_names = {
        "light_b_001.fits",
        "light_b_002.fits",
        "light_v_001.fits",
    }
    expected_files = [
        filename
        for filename in image_collection.files
        if filename in expected_names
    ]

    assert selected.files == expected_files


def test_filter_image_collection_treats_string_as_scalar(
    image_collection: ImageFileCollection,
) -> None:
    """A string criterion is one accepted value, not an iterable of letters."""
    selected = filter_image_collection(
        image_collection,
        imagetyp="LIGHT",
    )

    expected_names = {
        "light_b_001.fits",
        "light_b_002.fits",
        "light_no_filter_001.fits",
        "light_v_001.fits",
    }
    expected_files = [
        filename
        for filename in image_collection.files
        if filename in expected_names
    ]

    assert selected.files == expected_files


def test_filter_image_collection_accepts_tuple_criteria(
    image_collection: ImageFileCollection,
) -> None:
    """Tuple criteria are treated as multiple accepted values."""
    selected = filter_image_collection(
        image_collection,
        filter=("B", "V"),
    )

    expected_names = {
        "flat_b_001.fits",
        "light_b_001.fits",
        "light_b_002.fits",
        "light_v_001.fits",
    }

    assert set(selected.files) == expected_names


def test_filter_image_collection_rejects_missing_criteria(
    image_collection: ImageFileCollection,
) -> None:
    """Filtering requires at least one FITS-header criterion."""
    with pytest.raises(
        ValueError,
        match="requires at least one header criterion",
    ):
        filter_image_collection(image_collection)


def test_filter_image_collection_rejects_empty_allowed_values(
    image_collection: ImageFileCollection,
) -> None:
    """A criterion with no accepted values is invalid input."""
    with pytest.raises(
        ValueError,
        match="must contain at least one accepted value",
    ):
        filter_image_collection(
            image_collection,
            filter=[],
        )


def test_filter_image_collection_raises_when_no_fits_match(
    image_collection: ImageFileCollection,
) -> None:
    """A valid criterion with no matches raises the project-specific error."""
    with pytest.raises(
        NoMatchingImagesError,
        match="Cannot create an empty ImageFileCollection snapshot",
    ):
        filter_image_collection(
            image_collection,
            imagetyp="MASTER DARK",
        )


def test_filter_image_collection_rejects_invalid_collection() -> None:
    """The source object must be an ImageFileCollection."""
    with pytest.raises(TypeError, match="must be an ImageFileCollection"):
        filter_image_collection(
            object(),
            imagetyp="LIGHT",
        )


def test_filter_image_collection_rejects_invalid_header_name(
    image_collection: ImageFileCollection,
) -> None:
    """Header criterion names must be non-empty strings."""
    with pytest.raises(TypeError, match="non-empty string"):
        filter_image_collection(
            image_collection,
            **{"": "LIGHT"},
        )


def test_filter_image_collection_group_preserves_directory_boundaries(
    image_collection_group: ImageCollectionGroup,
) -> None:
    """Group filtering returns one snapshot selection per source directory."""
    selections = filter_image_collection_group(
        image_collection_group,
        imagetyp="LIGHT",
    )

    assert isinstance(selections, tuple)
    assert [directory for directory, _ in selections] == list(
        image_collection_group.directories
    )

    expected_files_by_directory = {
        image_collection_group.directories[0]: ["light_b_001.fits"],
        image_collection_group.directories[1]: ["light_v_001.fits"],
    }

    assert [
        (directory, selected_collection.files)
        for directory, selected_collection in selections
    ] == [
        (directory, expected_files_by_directory[directory])
        for directory in image_collection_group.directories
    ]


def test_filter_image_collection_group_propagates_no_matching_images_error(
    image_collection_group: ImageCollectionGroup,
) -> None:
    """Group filtering fails when any child collection has no matching files."""
    with pytest.raises(NoMatchingImagesError):
        filter_image_collection_group(
            image_collection_group,
            imagetyp="MASTER DARK",
        )


def test_filter_image_collection_group_rejects_invalid_group() -> None:
    """Group helpers require an ImageCollectionGroup."""
    with pytest.raises(TypeError, match="must be an ImageCollectionGroup"):
        filter_image_collection_group(
            object(),
            imagetyp="LIGHT",
        )


def test_get_header_values_reads_all_files_in_collection_order(
    image_collection: ImageFileCollection,
) -> None:
    """Header values follow the collection's current file order."""
    image_type_by_filename = {
        "bias_001.fits": "BIAS",
        "flat_b_001.fits": "FLAT",
        "light_b_001.fits": "LIGHT",
        "light_b_002.fits": "LIGHT",
        "light_no_filter_001.fits": "LIGHT",
        "light_v_001.fits": "LIGHT",
    }

    expected_values = [
        image_type_by_filename[filename]
        for filename in image_collection.files
    ]

    assert get_header_values(image_collection, "IMAGETYP") == expected_values


def test_get_header_values_unique_preserves_first_occurrence_order(
    image_collection: ImageFileCollection,
) -> None:
    """Unique extraction removes later duplicates without reordering."""
    all_values = get_header_values(image_collection, "IMAGETYP")

    expected_unique_values: list[object] = []
    for value in all_values:
        if value not in expected_unique_values:
            expected_unique_values.append(value)

    assert get_header_values(
        image_collection,
        "IMAGETYP",
        unique=True,
    ) == expected_unique_values


def test_get_header_values_handles_missing_header_policies(
    image_collection: ImageFileCollection,
) -> None:
    """Missing FITS headers obey raise, ignore, and none policies."""
    light_collection = filter_image_collection(
        image_collection,
        imagetyp="LIGHT",
    )
    filter_by_filename = {
        "light_b_001.fits": "B",
        "light_b_002.fits": "B",
        "light_no_filter_001.fits": None,
        "light_v_001.fits": "V",
    }

    expected_none_values = [
        filter_by_filename[filename]
        for filename in light_collection.files
    ]
    expected_ignored_values = [
        value
        for value in expected_none_values
        if value is not None
    ]

    with pytest.raises(KeyError, match="FILTER"):
        get_header_values(light_collection, "FILTER")

    assert get_header_values(
        light_collection,
        "FILTER",
        missing="none",
    ) == expected_none_values

    assert get_header_values(
        light_collection,
        "FILTER",
        missing="ignore",
    ) == expected_ignored_values


def test_get_header_values_rejects_invalid_missing_policy(
    image_collection: ImageFileCollection,
) -> None:
    """The missing-header policy is restricted to known values."""
    with pytest.raises(ValueError, match="'missing' must be one of"):
        get_header_values(
            image_collection,
            "FILTER",
            missing="skip",
        )


def test_get_header_values_rejects_invalid_header_name(
    image_collection: ImageFileCollection,
) -> None:
    """The requested FITS-header keyword must be a non-empty string."""
    with pytest.raises(TypeError, match="non-empty string"):
        get_header_values(
            image_collection,
            "",
        )


def test_get_header_values_rejects_invalid_collection() -> None:
    """The source object must be an ImageFileCollection."""
    with pytest.raises(TypeError, match="must be an ImageFileCollection"):
        get_header_values(
            object(),
            "FILTER",
        )


def test_get_group_header_values_preserves_directory_boundaries(
    image_collection_group: ImageCollectionGroup,
) -> None:
    """Group header extraction returns ordered directory-value pairs."""
    values_by_directory = get_group_header_values(
        image_collection_group,
        "IMAGETYP",
    )

    assert [directory for directory, _ in values_by_directory] == list(
        image_collection_group.directories
    )

    expected_values_by_directory = [
        get_header_values(collection, "IMAGETYP")
        for collection in image_collection_group
    ]

    assert [
        values for _, values in values_by_directory
    ] == expected_values_by_directory


def test_get_group_header_values_forwards_unique_option(
    image_collection_group: ImageCollectionGroup,
) -> None:
    """The unique option is applied independently within each directory."""
    values_by_directory = get_group_header_values(
        image_collection_group,
        "IMAGETYP",
        unique=True,
    )

    expected = tuple(
        (
            directory,
            get_header_values(
                collection,
                "IMAGETYP",
                unique=True,
            ),
        )
        for directory, collection in image_collection_group.items()
    )

    assert values_by_directory == expected


def test_get_group_header_values_rejects_invalid_group() -> None:
    """Group header extraction requires an ImageCollectionGroup."""
    with pytest.raises(TypeError, match="must be an ImageCollectionGroup"):
        get_group_header_values(
            object(),
            "FILTER",
        )


def test_refresh_image_collection_group_refreshes_only_selected_child(
    image_collection_group: ImageCollectionGroup,
    write_test_fits: Callable[..., Path],
) -> None:
    """Refreshing one group index does not refresh unselected children."""
    first_collection = image_collection_group[0]
    second_collection = image_collection_group[1]
    first_original_files = set(first_collection.files)
    second_original_files = set(second_collection.files)

    first_directory = Path(first_collection.location)
    second_directory = Path(second_collection.location)

    write_test_fits(
        first_directory / "master_bias.fits",
        imagetyp="MASTER BIAS",
        exptime=0.0,
    )
    write_test_fits(
        second_directory / "master_dark.fits",
        imagetyp="MASTER DARK",
        exptime=0.0,
    )

    returned_group = refresh_image_collection_group(
        image_collection_group,
        indices=1,
    )

    assert returned_group is image_collection_group
    assert set(first_collection.files) == first_original_files
    assert "master_bias.fits" not in first_collection.files
    assert set(second_collection.files) == second_original_files | {
        "master_dark.fits",
    }


def test_refresh_image_collection_group_supports_negative_indices(
    image_collection_group: ImageCollectionGroup,
    write_test_fits: Callable[..., Path],
) -> None:
    """Negative indices select child collections with Python index semantics."""
    last_collection = image_collection_group[-1]
    last_directory = Path(last_collection.location)

    write_test_fits(
        last_directory / "master_flat.fits",
        imagetyp="MASTER FLAT",
        filter_name="V",
        exptime=5.0,
    )

    refresh_image_collection_group(
        image_collection_group,
        indices=-1,
    )

    assert "master_flat.fits" in last_collection.files


def test_refresh_image_collection_group_supports_slice_indices(
    image_collection_group: ImageCollectionGroup,
    write_test_fits: Callable[..., Path],
) -> None:
    """Slice indices refresh only the selected child collections."""
    first_collection = image_collection_group[0]
    second_collection = image_collection_group[1]

    write_test_fits(
        Path(first_collection.location) / "first_later.fits",
        imagetyp="LIGHT",
    )
    write_test_fits(
        Path(second_collection.location) / "second_later.fits",
        imagetyp="LIGHT",
    )

    refresh_image_collection_group(
        image_collection_group,
        indices=slice(0, 1),
    )

    assert "first_later.fits" in first_collection.files
    assert "second_later.fits" not in second_collection.files


def test_refresh_image_collection_group_refreshes_duplicate_indices_once(
    image_collection_group: ImageCollectionGroup,
    write_test_fits: Callable[..., Path],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Duplicate indices are normalized before refresh calls are made."""
    first_collection = image_collection_group[0]
    write_test_fits(
        Path(first_collection.location) / "first_later.fits",
        imagetyp="LIGHT",
    )

    original_refresh = first_collection.refresh
    call_count = 0

    def counted_refresh() -> None:
        nonlocal call_count
        call_count += 1
        original_refresh()

    monkeypatch.setattr(first_collection, "refresh", counted_refresh)

    refresh_image_collection_group(
        image_collection_group,
        indices=[0, 0, -2],
    )

    assert call_count == 1
    assert "first_later.fits" in first_collection.files


@pytest.mark.parametrize(
    "indices",
    [
        "0",
        b"0",
        True,
        [False],
        [1.5],
    ],
)
def test_refresh_image_collection_group_rejects_invalid_indices(
    image_collection_group: ImageCollectionGroup,
    indices: object,
) -> None:
    """Invalid group index specifications raise TypeError."""
    with pytest.raises(TypeError):
        refresh_image_collection_group(
            image_collection_group,
            indices=indices,
        )


def test_refresh_image_collection_group_rejects_out_of_range_index(
    image_collection_group: ImageCollectionGroup,
) -> None:
    """Out-of-range group indices raise IndexError."""
    with pytest.raises(IndexError, match="out of range"):
        refresh_image_collection_group(
            image_collection_group,
            indices=99,
        )


def test_refresh_image_collection_group_rejects_invalid_group() -> None:
    """Group refresh requires an ImageCollectionGroup."""
    with pytest.raises(TypeError, match="must be an ImageCollectionGroup"):
        refresh_image_collection_group(object())


def test_filter_image_collection_treats_bytes_as_scalar(
    image_collection: ImageFileCollection,
) -> None:
    """Bytes criteria are scalar values, not iterables of integers."""
    with pytest.raises(NoMatchingImagesError):
        filter_image_collection(
            image_collection,
            imagetyp=b"LIGHT",
        )


def test_filter_image_collection_accepts_generator_criteria(
    image_collection: ImageFileCollection,
) -> None:
    """Non-string iterables are expanded into accepted values."""
    accepted_filters = (filter_name for filter_name in ["B", "V"])

    selected = filter_image_collection(
        image_collection,
        filter=accepted_filters,
    )

    expected_names = {
        "flat_b_001.fits",
        "light_b_001.fits",
        "light_b_002.fits",
        "light_v_001.fits",
    }

    assert set(selected.files) == expected_names


def test_filter_image_collection_accepts_scalar_numeric_criteria(
    image_collection: ImageFileCollection,
) -> None:
    """Non-iterable scalar criteria become one accepted value."""
    selected = filter_image_collection(
        image_collection,
        exptime=0.0,
    )

    assert selected.files == ["bias_001.fits"]


def test_refresh_image_collection_group_refreshes_all_when_indices_is_none(
    image_collection_group: ImageCollectionGroup,
    write_test_fits: Callable[..., Path],
) -> None:
    """indices=None refreshes every child collection."""
    first_collection = image_collection_group[0]
    second_collection = image_collection_group[1]

    write_test_fits(
        Path(first_collection.location) / "first_later.fits",
        imagetyp="LIGHT",
    )
    write_test_fits(
        Path(second_collection.location) / "second_later.fits",
        imagetyp="LIGHT",
    )

    refresh_image_collection_group(image_collection_group)

    assert "first_later.fits" in first_collection.files
    assert "second_later.fits" in second_collection.files
