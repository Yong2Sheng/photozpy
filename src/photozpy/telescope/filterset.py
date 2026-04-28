from __future__ import annotations

from collections.abc import Mapping
from types import MappingProxyType
from typing import TypeAlias

from tabulate import tabulate


AliasInput: TypeAlias = None | str | tuple[str, ...] | list[str]
AliasTuple: TypeAlias = tuple[str, ...]

EMPTY_VALUES = (None, "", (), [])

class UnknownFilterAliasError(ValueError):
    """Raised when a filter alias cannot be resolved."""

class FilterSet:
    """
    Define standard filter names and aliases.

    Lookup directions:

        standard filter name -> aliases
        alias -> standard filter name

    This object is intended to be immutable after construction.
    Add new aliases by editing the corresponding filter-set module.
    """

    def __init__(
        self,
        filterset_name: str,
        standard_to_alias: Mapping[str, AliasInput],
    ) -> None:

        self._filterset_name = filterset_name

        raw_standard_to_alias = standard_to_alias
        processed_standard_to_alias: dict[str, AliasTuple] = {}
        alias_to_standard: dict[str, str] = {}

        for standard_name, alias in raw_standard_to_alias.items():
            # first deal with empty values
            if alias in EMPTY_VALUES:
                # No extra aliases. The standard name itself is still added to
                # alias_to_standard below, so resolve("SDSS_g") will work.
                alias_tuple: AliasTuple = ()
                
            # second deal with string inputs
            elif isinstance(alias, str):
                if standard_name == alias:
                    raise ValueError(
                        f"Alias {alias!r} is the same as standard name {standard_name!r}. "
                        "Do not include the standard name in aliases; it is added automatically."
                    )
                else:
                    alias_tuple: AliasTuple  = (alias,)
                    
            # third deal with tuple and list inputs  
            elif isinstance(alias, (tuple, list)):
                if standard_name in alias:
                    raise ValueError(
                        f"Aliases for {standard_name!r} contain the standard name itself. "
                        "Do not include the standard name in aliases; it is added automatically."
                    )
                else:
                    alias_tuple: AliasTuple = tuple(alias)
                    
            # raise an error if the input type is not expected
            else:
                raise ValueError(f"Unsupported alias input type: {type(alias)}!")

            # standard_to_alias stores only user-declared aliases.
            # The standard name itself is intentionally not inserted here,
            # so this table stays close to the filter-set module definition.
            processed_standard_to_alias[standard_name] = alias_tuple
            
            # Build the reverse lookup table:
            #
            #     alias -> standard filter name
            #
            # Each alias must appear only once in the whole filter set.
            # This is intentionally strict. Even if the duplicated alias points to
            # the same standard filter name, it is still treated as an error because
            # the alias table should be clean and explicit.
            #
            # Example of a real conflict:
            #
            #     "SDSS_g": ("g", "sdss_g")
            #     "SDSS_r": ("g", "sdss_r")
            #
            # Here alias "g" would point to both "SDSS_g" and "SDSS_r".
            #
            # Example of redundant data:
            #
            #     "SDSS_g": ("g", "g", "sdss_g")
            #
            # This is not ambiguous, but it is still unnecessary duplication.
            # Since aliases are maintained manually in the filter-set modules, we
            # reject both conflicts and redundant duplicates early.
            #
            # Without this check, assigning to the same dict key again would silently
            # overwrite the previous value.
            #
            # alias_to_standard is a derived lookup table.
            # Here we add the standard name itself so resolve("SDSS_g")
            # works even after the FITS header has already been standardized.
            alias_to_standard[standard_name] = standard_name
            for alias in alias_tuple:
                if alias in alias_to_standard:
                    raise ValueError(
                        f"Duplicated alias {alias!r}. "
                        f"It was already mapped to {alias_to_standard[alias]!r}, "
                        f"but now maps to {standard_name!r}."
                    )
                alias_to_standard[alias] = standard_name
        
        self._standard_to_alias = MappingProxyType(processed_standard_to_alias)
        self._alias_to_standard = MappingProxyType(alias_to_standard)
        
        return

    @property
    def filterset_name(self) -> str:
        return self._filterset_name

    @property
    def standard_names(self) -> tuple[str, ...]:
        return tuple(self._standard_to_alias.keys())

    @property
    def alias(self) -> tuple[str, ...]:
        return tuple(self._alias_to_standard.keys())

    @property
    def standard_to_alias(self) -> Mapping[str, AliasTuple]:
        return self._standard_to_alias

    @property
    def alias_to_standard(self) -> Mapping[str, str]:
        return self._alias_to_standard

    def alias_of(
        self,
        standard_name: str,
    ) -> AliasTuple:
        try:
            return self._standard_to_alias[standard_name]
        except KeyError:
            raise KeyError(
                f"Unknown standard filter name {standard_name!r} "
                f"in filter set {self.standard_names!r}"
            ) from None

    def resolve(
        self,
        alias: str,
    ) -> str:
        try:
            return self._alias_to_standard[alias]
        except KeyError:
            raise UnknownFilterAliasError(
                f"Unknown filter alias {alias!r} in filter set {self.filterset_name!r}.\n\n"
                f"Available filters and aliases:\n"
                f"{self.to_table()}"
            ) from None

    def to_table(self) -> str:
        rows: list[list[object]] = [
            [key, value]
            for key, value in self._standard_to_alias.items()
        ]
        return tabulate(rows, headers=["Standard Filter Names", "Alias"], tablefmt="fancy_grid")

    def print_summary(self) -> None:
        print(self.to_table())