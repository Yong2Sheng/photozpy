"""
Helper functions for working with ccdproc.ImageFileCollection objects.

This module contains operations performed on image collections, while
ImageCollectionGroup remains a lightweight container for one or more
directory-backed ImageFileCollection objects.

Design principles
-----------------
1. A primary ImageFileCollection is normally directory-backed:

       ImageFileCollection(location=directory, filenames=None, ...)

   Calling ``refresh()`` on such a collection is delegated to ccdproc. It is
   expected to update the collection according to its own directory-discovery
   configuration and rebuild its header summary.

2. A filtered collection is a snapshot selection:

       filter_image_collection(primary_collection, imagetyp="LIGHT")

   The result represents the files selected at that moment. It should not be
   treated as a live directory view that automatically discovers later-created
   FITS products.

3. ImageCollectionGroup is intended to contain primary directory-backed
   collections. Group-level filtering therefore returns directory-labelled
   snapshot selections instead of another ImageCollectionGroup.

4. This module does not modify FITS files, delete files from disk, or manually
   mutate ImageFileCollection internal state.
"""


from __future__ import annotations

from collections.abc import Iterable, Sequence
from itertools import product
from operator import index as operator_index
from pathlib import Path
from typing import Literal, TypeAlias

from ccdproc import ImageFileCollection

from .image_collection_group import ImageCollectionGroup

MissingHeaderPolicy: TypeAlias = Literal["raise", "ignore", "none"]


class NoMatchingImagesError(ValueError):
    """
    Raised when FITS-header filtering selects no images.
    """


GroupCollectionPairs: TypeAlias = tuple[
    tuple[Path, ImageFileCollection],
    ...,
]


__all__ = [
    "NoMatchingImagesError",
    "filter_image_collection",
    "filter_image_collection_group",
    "get_group_header_values",
    "get_header_values",
    "refresh_image_collection",
    "refresh_image_collection_group",
]


def refresh_image_collection(
    collection: ImageFileCollection,
) -> ImageFileCollection:
    """
    Refresh one ImageFileCollection in place.

    For a directory-backed collection, this delegates file discovery and
    header-summary rebuilding to ``ccdproc.ImageFileCollection.refresh()``.

    For an explicit-filename collection, the same ccdproc refresh operation is
    used, but file discovery remains limited to its explicit filename list.

    Parameters
    ----------
    collection
        The collection to refresh.

    Returns
    -------
    ccdproc.ImageFileCollection
        The same object passed in, after in-place refresh.

    Raises
    ------
    TypeError
        If ``collection`` is not an ImageFileCollection.

    Notes
    -----
    This function intentionally does not reconstruct a new ImageFileCollection.
    Reconstructing one from partial public attributes can silently lose original
    construction settings or change collection semantics.
    """
    _validate_image_collection(collection, argument_name="collection")

    collection.refresh()

    return collection


def refresh_image_collection_group(
    group: ImageCollectionGroup,
    indices: int | slice | Iterable[int] | None = None,
) -> ImageCollectionGroup:
    """
    Refresh selected child collections in an ImageCollectionGroup in place.

    Parameters
    ----------
    group
        The ImageCollectionGroup whose child collections should be refreshed.

    indices
        Child collection index or indices to refresh.

        - ``None``: refresh all child collections.
        - ``int``: refresh one child collection.
        - ``slice``: refresh the selected slice of child collections.
        - iterable of int: refresh the specified child collections.

        Negative indices are supported. Duplicate indices are refreshed only
        once, preserving their first-occurrence order.

    Returns
    -------
    ImageCollectionGroup
        The same group object passed in, after in-place refresh.

    Raises
    ------
    TypeError
        If ``group`` is not an ImageCollectionGroup, or if an index is not an
        integer.

    IndexError
        If an index is outside the valid range.

    Examples
    --------
    Refresh all collections:

    >>> refresh_image_collection_group(group)

    Refresh only the first collection:

    >>> refresh_image_collection_group(group, indices=0)

    Refresh the first and third collections:

    >>> refresh_image_collection_group(group, indices=[0, 2])
    """
    _validate_image_collection_group(group)

    selected_indices = _normalize_group_indices(
        indices=indices,
        group_size=len(group),
    )

    for collection_index in selected_indices:
        refresh_image_collection(group[collection_index])

    return group


def filter_image_collection(
    collection: ImageFileCollection,
    **header_criteria: object,
) -> ImageFileCollection:
    """
    Create a snapshot selection from one ImageFileCollection.

    Values for the same header are combined with logical OR. Constraints for
    different headers are combined with logical AND.

    The source collection is expected to use filenames in the same form
    returned by ``collection.files_filtered()``. Photozpy primary collections
    satisfy this because they are directory-backed collections created with
    ``filenames=None``.

    For example:

    >>> selected = filter_image_collection(
    ...     collection,
    ...     imagetyp=["LIGHT", "BIAS"],
    ...     filter=["B", "V"],
    ... )

    corresponds conceptually to:

    .. code-block:: text

        (IMAGETYP == "LIGHT" OR IMAGETYP == "BIAS")
        AND
        (FILTER == "B" OR FILTER == "V")

    Parameters
    ----------
    collection
        Source collection to filter.

    **header_criteria
        FITS-header selection criteria. Each value may be either:

        - a scalar value, such as ``imagetyp="LIGHT"``;
        - an iterable of accepted values, such as
          ``filter=["B", "V", "R"]``.

        Strings and bytes are always treated as scalar values, not as
        iterables of characters.

    Returns
    -------
    ccdproc.ImageFileCollection
        A newly created explicit-filename snapshot selection.

    Raises
    ------
    TypeError
        If ``collection`` is not an ImageFileCollection.

    ValueError
        If no header criteria are supplied.

    NoMatchingImagesError
        If no FITS files match the requested header criteria.

    Notes
    -----
    The returned collection is intentionally an explicit-filename collection.
    It preserves the selected filenames, their order, the source location,
    source keywords, and source extension.

    It should not be expected to discover FITS files created later in the
    directory. Refresh the primary directory-backed collection and filter it
    again when a new selection is needed.
    """
    _validate_image_collection(collection, argument_name="collection")

    if not header_criteria:
        raise ValueError(
            "filter_image_collection() requires at least one header criterion."
        )

    criteria_combinations = _expand_header_criteria(header_criteria)

    matched_filenames: set[str] = set()

    for criteria in criteria_combinations:
        matched_filenames.update(
            str(filename)
            for filename in collection.files_filtered(**criteria)
        )

    # ``files_filtered`` is called once for each header-value combination.
    # Reconstruct the result in the original collection order so that the
    # output does not depend on combination order or unordered set iteration.
    ordered_selected_filenames = [
        filename
        for filename in collection.files
        if str(filename) in matched_filenames
    ]

    return _new_snapshot_collection(
        source_collection=collection,
        filenames=ordered_selected_filenames,
    )


def filter_image_collection_group(
    group: ImageCollectionGroup,
    **header_criteria: object,
) -> GroupCollectionPairs:
    """
    Filter each child collection in a group independently.

    Parameters
    ----------
    group
        Source ImageCollectionGroup.

    **header_criteria
        FITS-header selection criteria forwarded to
        :func:`filter_image_collection`.

    Returns
    -------
    tuple[tuple[pathlib.Path, ccdproc.ImageFileCollection], ...]
        Ordered ``(directory, selected_collection)`` pairs.

        The result is deliberately not an ImageCollectionGroup. Each selected
        collection is a snapshot selection rather than a primary directory-
        backed collection that should be refreshed to discover new FITS files.

    Examples
    --------
    >>> selections = filter_image_collection_group(
    ...     group,
    ...     imagetyp="LIGHT",
    ... )

    >>> for directory, light_collection in selections:
    ...     print(directory, len(light_collection.files))
    """
    _validate_image_collection_group(group)

    return tuple(
        (
            directory,
            filter_image_collection(
                collection,
                **header_criteria,
            ),
        )
        for directory, collection in group.items()
    )


def get_header_values(
    collection: ImageFileCollection,
    header: str,
    *,
    unique: bool = False,
    missing: MissingHeaderPolicy = "raise",
) -> list[object]:
    """
    Extract one FITS-header value from every file in a collection.

    Parameters
    ----------
    collection
        Source collection.

    header
        FITS-header keyword to read.

    unique
        If True, remove repeated values while preserving first-occurrence
        order. Default is False.

    missing
        Policy for files whose headers do not contain ``header``.

        - ``"raise"``: raise KeyError.
        - ``"ignore"``: omit that file from the returned values.
        - ``"none"``: append None for that file.

    Returns
    -------
    list
        Header values in the collection's current file order.

    Raises
    ------
    TypeError
        If ``collection`` is invalid, ``header`` is not a non-empty string,
        or ``missing`` is invalid.

    KeyError
        If a file is missing the requested header and ``missing="raise"``.

    Notes
    -----
    Headers are read through ``collection.headers()`` so that lookup follows
    the collection's configured FITS extension and current file order. This
    function does not rely on ``collection.summary``, whose available keyword
    columns depend on how the collection was constructed.
    """
    _validate_image_collection(collection, argument_name="collection")
    _validate_header_name(header)
    _validate_missing_policy(missing)

    values: list[object] = []

    for fits_header, filename in collection.headers(return_fname=True):
        if header in fits_header:
            values.append(fits_header[header])
            continue

        if missing == "raise":
            raise KeyError(
                f"Header {header!r} is missing from FITS file {filename!r}."
            )

        if missing == "none":
            values.append(None)

    if unique:
        return _unique_preserving_order(values)

    return values


def get_group_header_values(
    group: ImageCollectionGroup,
    header: str,
    *,
    unique: bool = False,
    missing: MissingHeaderPolicy = "raise",
) -> tuple[tuple[Path, list[object]], ...]:
    """
    Extract header values from every child collection while preserving directory
    boundaries.

    Parameters
    ----------
    group
        Source ImageCollectionGroup.

    header
        FITS-header keyword to read.

    unique
        Passed to :func:`get_header_values`.

    missing
        Passed to :func:`get_header_values`.

    Returns
    -------
    tuple[tuple[pathlib.Path, list], ...]
        Ordered ``(directory, header_values)`` pairs.

    Notes
    -----
    A tuple of pairs is used instead of ``dict[Path, list]`` because
    ImageCollectionGroup permits repeated directory entries. A dictionary
    would silently discard earlier entries for duplicate paths.

    Examples
    --------
    >>> values_by_directory = get_group_header_values(
    ...     group,
    ...     header="FILTER",
    ...     unique=True,
    ... )

    >>> for directory, filters in values_by_directory:
    ...     print(directory, filters)
    """
    _validate_image_collection_group(group)

    return tuple(
        (
            directory,
            get_header_values(
                collection,
                header,
                unique=unique,
                missing=missing,
            ),
        )
        for directory, collection in group.items()
    )


def _new_snapshot_collection(
    source_collection: ImageFileCollection,
    filenames: Sequence[str],
) -> ImageFileCollection:
    """
    Construct a non-empty explicit-filename snapshot collection.

    Snapshot collections represent a concrete selection of existing FITS files.
    An empty selection is considered an invalid pipeline state and is rejected.
    """
    selected_filenames = list(filenames)

    if not selected_filenames:
        raise NoMatchingImagesError(
            "Cannot create an empty ImageFileCollection snapshot. "
            "No FITS files matched the requested selection criteria."
        )

    location = source_collection.location

    collection_kwargs: dict[str, object] = {
        "filenames": selected_filenames,
        "keywords": source_collection.keywords,
        "ext": source_collection.ext,
    }

    if location:
        collection_kwargs["location"] = location

    return ImageFileCollection(**collection_kwargs)


def _expand_header_criteria(
    header_criteria: dict[str, object],
) -> tuple[dict[str, object], ...]:
    """
    Expand multi-value header criteria into scalar-value combinations.

    Example
    -------
    Input:

    >>> {
    ...     "imagetyp": ["LIGHT", "BIAS"],
    ...     "filter": ["B", "V"],
    ... }

    Output:

    >>> (
    ...     {"imagetyp": "LIGHT", "filter": "B"},
    ...     {"imagetyp": "LIGHT", "filter": "V"},
    ...     {"imagetyp": "BIAS", "filter": "B"},
    ...     {"imagetyp": "BIAS", "filter": "V"},
    ... )
    """
    normalized_values: list[tuple[object, ...]] = []
    headers: list[str] = []

    for header, value in header_criteria.items():
        _validate_header_name(header)

        normalized_value = _normalize_criterion_values(value)

        if not normalized_value:
            raise ValueError(
                f"Header criterion '{header}' must contain at least one "
                "accepted value."
            )

        headers.append(header)
        normalized_values.append(normalized_value)

    return tuple(
        dict(zip(headers, values, strict=True))
        for values in product(*normalized_values)
    )


def _normalize_criterion_values(
    value: object,
) -> tuple[object, ...]:
    """
    A scalar value, list, tuple, set, or other iterable is normalized into a
    tuple of accepted values. (无论你传入的是单个值、list、tuple、set，
    最后都得到一个 tuple。)
    Scalar values become one-element tuples. Lists, tuples, sets, and other
    non-string iterables become tuples containing their elements.

    Strings and bytes remain scalar values. Other iterables represent multiple
    accepted values.
    """
    if isinstance(value, (str, bytes)):
        return (value,)

    if isinstance(value, Iterable):
        return tuple(value)

    return (value,)


def _normalize_group_indices(
    indices: int | slice | Iterable[int] | None,
    group_size: int,
) -> tuple[int, ...]:
    """
    Normalize a group index specification into unique valid integer indices.
    """
    if indices is None:
        return tuple(range(group_size))

    if isinstance(indices, slice):
        return tuple(range(group_size)[indices])

    if isinstance(indices, (str, bytes)):
        raise TypeError(
            "'indices' must be an int, slice, iterable of int, or None."
        )

    if _is_index_like(indices):
        raw_indices = (indices,)
    else:
        try:
            raw_indices = tuple(indices)
        except TypeError as exc:
            raise TypeError(
                "'indices' must be an int, slice, iterable of int, or None."
            ) from exc

    normalized_indices: list[int] = []

    for raw_index in raw_indices:
        normalized_index = _normalize_single_group_index(
            raw_index,
            group_size=group_size,
        )

        if normalized_index not in normalized_indices:
            normalized_indices.append(normalized_index)

    return tuple(normalized_indices)


def _normalize_single_group_index(
    raw_index: object,
    *,
    group_size: int,
) -> int:
    """
    Convert one index-like object into a validated non-negative index.
    """
    if isinstance(raw_index, bool):
        raise TypeError("Boolean values are not valid collection indices.")

    try:
        normalized_index = operator_index(raw_index)
    except TypeError as exc:
        raise TypeError(
            f"Collection index must be an integer, got {raw_index!r}."
        ) from exc

    if normalized_index < 0:
        normalized_index += group_size

    if not 0 <= normalized_index < group_size:
        raise IndexError(
            f"Collection index {raw_index!r} is out of range "
            f"for a group of length {group_size}."
        )

    return normalized_index


def _is_index_like(value: object) -> bool:
    """
    Return True when ``value`` can be interpreted as a Python integer index.
    """
    if isinstance(value, bool):
        return False

    try:
        operator_index(value)
    except TypeError:
        return False

    return True


def _unique_preserving_order(
    values: Iterable[object],
) -> list[object]:
    """
    Remove duplicate values without relying on hashability.

    FITS header values are normally hashable scalars, but this implementation
    intentionally avoids assuming that every possible header value is hashable.
    """
    unique_values: list[object] = []

    for value in values:
        if value not in unique_values:
            unique_values.append(value)

    return unique_values


def _validate_image_collection(
    collection: object,
    *,
    argument_name: str,
) -> None:
    """
    Validate one ImageFileCollection argument.
    """
    if not isinstance(collection, ImageFileCollection):
        raise TypeError(
            f"'{argument_name}' must be an ImageFileCollection, "
            f"got {type(collection).__name__}."
        )


def _validate_image_collection_group(
    group: object,
) -> None:
    """
    Validate one ImageCollectionGroup argument.
    """
    if not isinstance(group, ImageCollectionGroup):
        raise TypeError(
            "'group' must be an ImageCollectionGroup, "
            f"got {type(group).__name__}."
        )


def _validate_header_name(
    header: object,
) -> None:
    """
    Validate one FITS-header keyword.
    """
    if not isinstance(header, str) or not header.strip():
        raise TypeError(
            "'header' must be a non-empty string."
        )


def _validate_missing_policy(
    missing: object,
) -> None:
    """
    Validate the missing-header policy.
    """
    valid_policies = {"raise", "ignore", "none"}

    if missing not in valid_policies:
        raise ValueError(
            "'missing' must be one of "
            f"{sorted(valid_policies)!r}, got {missing!r}."
        )
