from pathlib import Path
from collections.abc import Iterable

from ccdproc import ImageFileCollection


class ImageCollectionGroup:
    """
    A group of one or more `ccdproc.ImageFileCollection` objects.

    This class provides a lightweight container for working with image
    collections from multiple directories. Each child collection represents
    one directory and remains an independent `ImageFileCollection`.

    The group can be created either from directories that share the same
    `ImageFileCollection` configuration, or from existing, independently
    configured `ImageFileCollection` objects.

    Parameters
    ----------
    directories : str, pathlib.Path, or iterable of str or pathlib.Path
        One directory or an iterable of directories used to create child
        `ImageFileCollection` objects. Relative paths and ``~`` are converted
        to normalized absolute paths before child collections are created.

    **image_file_collection_kwargs
        Keyword arguments forwarded unchanged to each child
        `ccdproc.ImageFileCollection`, except ``location`` and ``filenames``.

        ``location`` is controlled by ``directories``. ``filenames`` is not
        accepted because one shared filename list does not have an unambiguous
        meaning across multiple directories.

    Raises
    ------
    TypeError
        If ``directories`` is neither a path-like object nor an iterable of
        path-like objects.

        Also raised if ``location`` or ``filenames`` is passed through
        ``image_file_collection_kwargs``.

    ValueError
        If ``directories`` is empty.

    Notes
    -----
    The group stores child collections in the same order as their corresponding
    directories. Therefore, for every valid index ``i``:

    ``self._directories[i]`` corresponds to ``self._collections[i]``.

    The group does not merge child collections into one large
    `ImageFileCollection`. Each child retains its own configuration, summary
    table, filtering behavior, and refresh state.

    See Also
    --------
    ccdproc.ImageFileCollection
        The single-directory collection class used by each child collection.

    ImageCollectionGroup.from_collections
        Create a group from existing `ImageFileCollection` objects with
        potentially different configurations.
    """
    
    def __init__(
        self,
        directories: str | Path | Iterable[str | Path],
        **image_file_collection_kwargs,
    ) -> None:
        """
        Create child `ImageFileCollection` objects from one or more directories.

        A separate `ImageFileCollection` is created for every input directory.
        All child collections receive the same keyword arguments.

        Parameters
        ----------
        directories : str, pathlib.Path, or iterable of str or pathlib.Path
            One directory or an iterable of directories. Each directory is
            normalized with ``Path(...).expanduser().resolve()`` before being
            passed to `ImageFileCollection`.

        **image_file_collection_kwargs
            Keyword arguments forwarded to every child
            `ccdproc.ImageFileCollection`.

            ``location`` and ``filenames`` are reserved and cannot be passed.
            ``location`` is determined by each item in ``directories``.
            ``filenames`` is intentionally unsupported because different
            directories may require different filename selections.

        Raises
        ------
        TypeError
            If ``directories`` is not a path-like object or an iterable of
            path-like objects.

            Also raised if ``location`` or ``filenames`` is included in
            ``image_file_collection_kwargs``.

        ValueError
            If ``directories`` is empty.
        """

        # These keyword arguments cannot be forwarded to ImageFileCollection,
        # because ImageCollectionGroup controls child collection locations.
        forbidden = {"location", "filenames"}
        invalid = forbidden & image_file_collection_kwargs.keys()
        if invalid:
            names = ", ".join(sorted(invalid))
            raise TypeError(
                f"{names} must not be passed to ImageCollectionGroup. "
                "Use 'directories' to define child collection locations."
            )

        # check directories type
        if isinstance(directories, (str, Path)):
            directories = (directories,)
        elif not isinstance(directories, Iterable):
            raise TypeError(
                "'directories' must be a str, pathlib.Path, "
                "or an iterable of str/pathlib.Path objects."
            )
        
        self._directories = tuple(
            Path(directory).expanduser().resolve()
            for directory in directories
        )
        
        # Reject an empty directory iterable.
        if not self._directories:
            raise ValueError(
                "'directories' must contain at least one directory."
            )

        # Create one child ImageFileCollection for each directory.
        self._collections = tuple(
            ImageFileCollection(
                location=directory,
                **image_file_collection_kwargs,
            )
            for directory in self._directories
        )

    @classmethod
    def from_collections(
        cls, 
        collections: ImageFileCollection | Iterable[ImageFileCollection],
    ) -> "ImageCollectionGroup":
        """
        Create an `ImageCollectionGroup` from existing image collections.

        This constructor preserves the original `ImageFileCollection` objects
        rather than reconstructing them from their locations. It is intended
        for cases where child collections use different configurations, such
        as different ``glob_include``, ``glob_exclude``, ``ext``, ``keywords``,
        or filename-selection rules.

        Parameters
        ----------
        collections : ccdproc.ImageFileCollection or iterable of ccdproc.ImageFileCollection
            One existing `ImageFileCollection` or an iterable of existing
            `ImageFileCollection` objects.

            Every collection must have a non-``None`` ``location`` attribute.
            This requirement ensures that the group can maintain a directory
            entry corresponding to every child collection.

        Returns
        -------
        ImageCollectionGroup
            A group containing the original `ImageFileCollection` objects.

            The returned group preserves object identity. For example, if
            ``collection_a`` is passed as the first item, then
            ``group[0] is collection_a`` is ``True``.

        Raises
        ------
        TypeError
            If ``collections`` is neither an `ImageFileCollection` nor an
            iterable of `ImageFileCollection` objects.

            Also raised if any iterable element is not an
            `ImageFileCollection`.

        ValueError
            If ``collections`` is empty.

            Also raised if any child collection has ``location is None``.

        Notes
        -----
        This method bypasses the normal ``__init__`` construction path so that
        existing child collections are retained exactly as provided.

        In particular, this method does not recreate child collections and
        does not attempt to infer or reproduce their original construction
        keyword arguments.
        """
        
        # initialize collections tuple
        if isinstance(collections, ImageFileCollection):
            collections = (collections,)
        elif not isinstance(collections, Iterable):
            raise TypeError(
                "'collections' must be an ImageFileCollection or an "
                "iterable of ImageFileCollection objects."
            )
        else:
            collections = tuple(collections)

        # make sure the input collection is not an empty list
        if not collections:
            raise ValueError(
                "'collections' must contain at least one ImageFileCollection."
            )

        # First check element types.
        if not all(
            isinstance(collection, ImageFileCollection)
            for collection in collections
        ):
            raise TypeError(
                "All items in 'collections' must be ImageFileCollection objects."
            )
        
        # Then check whether each valid collection has one location.
        for collection in collections:
            if collection.location is None:
                raise ValueError(
                    "Each ImageFileCollection passed to "
                    "ImageCollectionGroup.from_collections() "
                    "must have a non-None location."
                )

        directories = tuple(
            Path(collection.location).expanduser().resolve()
            for collection in collections
        )

        image_group_collection = cls.__new__(cls)
        image_group_collection._collections = collections
        image_group_collection._directories = directories

        return image_group_collection

    def __iter__(self):
        """
        Iterate over child `ImageFileCollection` objects.

        Yields
        ------
        ccdproc.ImageFileCollection
            Child collections in the same order as the input directories or
            input collections.

        Examples
        --------
        >>> for collection in group:
        ...     print(collection.location)
        """
        return iter(self._collections)

    def __len__(self):
        """
        Return the number of child image collections in the group.

        Returns
        -------
        int
            Number of stored `ImageFileCollection` objects.

        Examples
        --------
        >>> len(group)
        3
        """
        return len(self._collections)

    def __getitem__(self, index):
        """
        Parameters
        ----------
        index : int or slice
            Integer index or slice applied to the child collection tuple.
        
        Returns
        -------
        ccdproc.ImageFileCollection or ImageCollectionGroup
            A single child collection for an integer index, or a new
            `ImageCollectionGroup` for a non-empty slice.
        
        Raises
        ------
        ValueError
            If a slice selects no child collections. Empty groups are not valid
            `ImageCollectionGroup` instances.
        """
        result = self._collections[index]
    
        if isinstance(index, slice):
            return type(self).from_collections(result)
    
        return result

    def items(self):
        """
        Iterate over directory and child-collection pairs.

        Returns
        -------
        zip
            An iterator yielding ``(directory, collection)`` pairs, where
            ``directory`` is a normalized `pathlib.Path` and ``collection`` is
            the corresponding `ccdproc.ImageFileCollection`.

        Notes
        -----
        The returned object is a one-pass iterator. Convert it to a tuple or
        list if the pairs need to be reused:

        >>> pairs = tuple(group.items())

        Repeated directories are allowed. Therefore, the directory value should
        be treated as location metadata rather than a unique identifier for a
        child collection.

        Examples
        --------
        >>> for directory, collection in group.items():
        ...     print(directory, collection.location)
        """
        return zip(self._directories, self._collections)