from pathlib import Path

import pytest
from astropy.io import fits
from ccdproc import ImageFileCollection

from photozpy.collections import ImageCollectionGroup


def write_minimal_fits(path: Path) -> None:
    """Write one minimal valid FITS file."""
    fits.PrimaryHDU().writeto(path)


@pytest.fixture
def image_directories(tmp_path: Path) -> tuple[Path, Path]:
    """Create two directories containing minimal FITS files."""
    directory_a = tmp_path / "source_a"
    directory_b = tmp_path / "source_b"

    directory_a.mkdir()
    directory_b.mkdir()

    write_minimal_fits(directory_a / "a.fits")
    write_minimal_fits(directory_b / "b.fits")

    return directory_a, directory_b


@pytest.fixture
def image_collections(
    image_directories: tuple[Path, Path],
) -> tuple[ImageFileCollection, ImageFileCollection]:
    """Create two independently configured ImageFileCollection objects."""
    directory_a, directory_b = image_directories

    collection_a = ImageFileCollection(
        location=directory_a,
        glob_include="a.fits",
    )

    collection_b = ImageFileCollection(
        location=directory_b,
        glob_include="b.fits",
    )

    return collection_a, collection_b


def test_init_creates_one_collection_for_single_path(
    image_directories: tuple[Path, Path],
) -> None:
    directory_a, _ = image_directories

    group = ImageCollectionGroup(directory_a)

    assert len(group) == 1
    assert isinstance(group[0], ImageFileCollection)
    assert Path(group[0].location).resolve() == directory_a.resolve()


def test_init_creates_one_collection_for_single_string_path(
    image_directories: tuple[Path, Path],
) -> None:
    directory_a, _ = image_directories

    group = ImageCollectionGroup(str(directory_a))

    assert len(group) == 1
    assert Path(group[0].location).resolve() == directory_a.resolve()


def test_init_creates_collections_in_directory_order(
    image_directories: tuple[Path, Path],
) -> None:
    directory_a, directory_b = image_directories

    group = ImageCollectionGroup(
        [directory_a, directory_b],
    )

    assert len(group) == 2
    assert Path(group[0].location).resolve() == directory_a.resolve()
    assert Path(group[1].location).resolve() == directory_b.resolve()


def test_init_accepts_generator_of_directories(
    image_directories: tuple[Path, Path],
) -> None:
    directory_a, directory_b = image_directories

    directories = (
        directory
        for directory in (directory_a, directory_b)
    )

    group = ImageCollectionGroup(directories)

    assert len(group) == 2
    assert group._directories == (
        directory_a.resolve(),
        directory_b.resolve(),
    )


def test_init_forwards_image_file_collection_kwargs(
    tmp_path: Path,
) -> None:
    directory = tmp_path / "source"
    directory.mkdir()

    write_minimal_fits(directory / "raw.fits")
    write_minimal_fits(directory / "reduced.fits")

    group = ImageCollectionGroup(
        directory,
        glob_include="raw.fits",
    )

    assert group[0].files == ["raw.fits"]


def test_init_rejects_non_iterable_directories() -> None:
    with pytest.raises(
        TypeError,
        match="must be a str, pathlib.Path, or an iterable",
    ):
        ImageCollectionGroup(42)


@pytest.mark.parametrize(
    "directories",
    [
        [],
        (),
        (directory for directory in ()),
    ],
)
def test_init_rejects_empty_directory_iterable(
    directories,
) -> None:
    with pytest.raises(
        ValueError,
        match="must contain at least one directory",
    ):
        ImageCollectionGroup(directories)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"location": "other_directory"},
        {"filenames": ["a.fits"]},
    ],
)
def test_init_rejects_reserved_image_file_collection_kwargs(
    image_directories: tuple[Path, Path],
    kwargs: dict,
) -> None:
    directory_a, _ = image_directories

    with pytest.raises(
        TypeError,
        match="must not be passed to ImageCollectionGroup",
    ):
        ImageCollectionGroup(
            directory_a,
            **kwargs,
        )


def test_iter_returns_child_collections_in_order(
    image_collections: tuple[ImageFileCollection, ImageFileCollection],
) -> None:
    collection_a, collection_b = image_collections

    group = ImageCollectionGroup.from_collections(
        [collection_a, collection_b],
    )

    assert tuple(group) == (
        collection_a,
        collection_b,
    )


def test_len_returns_number_of_child_collections(
    image_collections: tuple[ImageFileCollection, ImageFileCollection],
) -> None:
    group = ImageCollectionGroup.from_collections(
        image_collections,
    )

    assert len(group) == 2


def test_from_collections_accepts_single_collection(
    image_collections: tuple[ImageFileCollection, ImageFileCollection],
) -> None:
    collection_a, _ = image_collections

    group = ImageCollectionGroup.from_collections(collection_a)

    assert len(group) == 1
    assert group[0] is collection_a


def test_from_collections_preserves_child_identity(
    image_collections: tuple[ImageFileCollection, ImageFileCollection],
) -> None:
    collection_a, collection_b = image_collections

    group = ImageCollectionGroup.from_collections(
        [collection_a, collection_b],
    )

    assert group[0] is collection_a
    assert group[1] is collection_b


def test_from_collections_accepts_generator_of_collections(
    image_collections: tuple[ImageFileCollection, ImageFileCollection],
) -> None:
    collection_a, collection_b = image_collections

    collections = (
        collection
        for collection in (collection_a, collection_b)
    )

    group = ImageCollectionGroup.from_collections(collections)

    assert len(group) == 2
    assert group[0] is collection_a
    assert group[1] is collection_b


def test_from_collections_allows_same_location_with_different_configs(
    tmp_path: Path,
) -> None:
    directory = tmp_path / "source"
    directory.mkdir()

    write_minimal_fits(directory / "raw.fits")
    write_minimal_fits(directory / "reduced.fits")

    raw_collection = ImageFileCollection(
        location=directory,
        glob_include="raw.fits",
    )

    reduced_collection = ImageFileCollection(
        location=directory,
        glob_include="reduced.fits",
    )

    group = ImageCollectionGroup.from_collections(
        [raw_collection, reduced_collection],
    )

    assert len(group) == 2
    assert group[0] is raw_collection
    assert group[1] is reduced_collection

    assert tuple(
        directory_
        for directory_, _ in group.items()
    ) == (
        directory.resolve(),
        directory.resolve(),
    )


@pytest.mark.parametrize(
    "collections",
    [
        [],
        (),
        (collection for collection in ()),
    ],
)
def test_from_collections_rejects_empty_iterable(
    collections,
) -> None:
    with pytest.raises(
        ValueError,
        match="must contain at least one ImageFileCollection",
    ):
        ImageCollectionGroup.from_collections(collections)


@pytest.mark.parametrize(
    "collections",
    [
        [None],
        ["not an image collection"],
        [ImageFileCollection, None],
    ],
)
def test_from_collections_rejects_non_collection_elements(
    collections,
) -> None:
    with pytest.raises(
        TypeError,
        match="All items in 'collections' must be ImageFileCollection",
    ):
        ImageCollectionGroup.from_collections(collections)


def test_getitem_returns_collection_for_integer_index(
    image_collections: tuple[ImageFileCollection, ImageFileCollection],
) -> None:
    collection_a, collection_b = image_collections

    group = ImageCollectionGroup.from_collections(
        [collection_a, collection_b],
    )

    assert group[0] is collection_a
    assert group[1] is collection_b


def test_getitem_returns_collection_for_negative_index(
    image_collections: tuple[ImageFileCollection, ImageFileCollection],
) -> None:
    collection_a, collection_b = image_collections

    group = ImageCollectionGroup.from_collections(
        [collection_a, collection_b],
    )

    assert group[-1] is collection_b
    assert group[-2] is collection_a


def test_getitem_returns_group_for_nonempty_slice(
    image_collections: tuple[ImageFileCollection, ImageFileCollection],
) -> None:
    collection_a, collection_b = image_collections

    group = ImageCollectionGroup.from_collections(
        [collection_a, collection_b],
    )

    subgroup = group[1:2]

    assert isinstance(subgroup, ImageCollectionGroup)
    assert len(subgroup) == 1
    assert subgroup[0] is collection_b
    assert subgroup[0] is not collection_a


def test_getitem_reverse_slice_preserves_child_identity(
    image_collections: tuple[ImageFileCollection, ImageFileCollection],
) -> None:
    collection_a, collection_b = image_collections

    group = ImageCollectionGroup.from_collections(
        [collection_a, collection_b],
    )

    reverse_group = group[::-1]

    assert isinstance(reverse_group, ImageCollectionGroup)
    assert reverse_group[0] is collection_b
    assert reverse_group[1] is collection_a


def test_getitem_raises_value_error_for_empty_slice(
    image_collections: tuple[ImageFileCollection, ImageFileCollection],
) -> None:
    group = ImageCollectionGroup.from_collections(
        image_collections,
    )

    with pytest.raises(
        ValueError,
        match="must contain at least one ImageFileCollection",
    ):
        group[0:0]


def test_items_returns_directory_collection_pairs_in_order(
    image_collections: tuple[ImageFileCollection, ImageFileCollection],
) -> None:
    collection_a, collection_b = image_collections

    group = ImageCollectionGroup.from_collections(
        [collection_a, collection_b],
    )

    items = tuple(group.items())

    assert items[0][0] == Path(collection_a.location).resolve()
    assert items[1][0] == Path(collection_b.location).resolve()

    assert items[0][1] is collection_a
    assert items[1][1] is collection_b


def test_from_collections_rejects_collection_without_location(
    tmp_path: Path,
) -> None:
    fits_path = tmp_path / "image.fits"
    write_minimal_fits(fits_path)

    collection = ImageFileCollection(
        filenames=[fits_path],
    )

    assert not collection.location

    with pytest.raises(
        ValueError,
        match="must have a non-empty location",
    ):
        ImageCollectionGroup.from_collections(collection)


def test_from_collections_rejects_non_iterable_input() -> None:
    with pytest.raises(
        TypeError,
        match="must be an ImageFileCollection or an iterable",
    ):
        ImageCollectionGroup.from_collections(42)
