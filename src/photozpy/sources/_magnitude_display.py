"""Notebook copy controls for vertical photometry spreadsheets."""

from __future__ import annotations

from collections.abc import Sequence
import csv
from html import escape
from io import StringIO
import math
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .sources import MagnitudeKind, Sources


def _tsv(rows: list[list[str]]) -> str:
    buffer = StringIO()
    writer = csv.writer(buffer, delimiter="\t", lineterminator="\n")
    for row in rows:
        if row == [""]:
            buffer.write("\n")
        else:
            writer.writerow(row)
    return buffer.getvalue()


def _copy_box(label: str, value: str, rows: int) -> str:
    # A real text field isolates selection from all surrounding columns.
    # Selection and Cmd/Ctrl+C also work when notebook scripts are disabled.
    return (
        '<div style="margin:8px 0">'
        f'<label style="display:block;font-weight:600">{escape(label)}'
        f'<textarea readonly aria-label="{escape(label, quote=True)}" '
        f'rows="{max(1, rows)}" wrap="off" spellcheck="false" '
        'style="display:block;box-sizing:border-box;width:100%;min-width:140px;'
        'font:14px/1.6 monospace;padding:8px;white-space:pre;resize:vertical"'
        f'>{escape(value)}</textarea></label>'
        '<button type="button" style="margin-top:4px" onclick="'
        'const t=this.previousElementSibling.querySelector(\'textarea\');'
        't.focus();t.select();'
        'if(navigator.clipboard &amp;&amp; window.isSecureContext){'
        'navigator.clipboard.writeText(t.value).catch(()=&gt;document.execCommand(\'copy\'));'
        '}else{document.execCommand(\'copy\');}">复制 / Copy</button></div>'
    )


def _number(result, filter_name: str, column: str) -> float:
    if result is None:
        return math.nan
    row = result.qtable[result.filter_names.index(filter_name)]
    value = row[column]
    return float(value.value) if column != "significance" else float(value)


def _text(value: float, precision: int = 4) -> str:
    return f"{value:.{precision}f}" if math.isfinite(value) else ""


def google_sheets_html(
    sources: Sources,
    kind: MagnitudeKind,
    filters: Sequence[str] | None,
    columns: Sequence[str] | None,
) -> str:
    """Build escaped HTML; all copy payloads contain data cells only."""
    panels = []
    for single in sources:
        info = single.source_info
        telescope = info.observation_field.telescope
        griz = tuple(f"SDSS_{b}'" for b in "griz")
        filter_names = (
            (griz if all(f in telescope.filters() for f in griz) else telescope.filters())
            if filters is None else tuple(telescope.resolve_filter_alias(f) for f in filters)
        )
        if len(set(filter_names)) != len(filter_names):
            raise ValueError("filters must resolve to distinct telescope filters.")
        coord = info.source_coordinate.icrs
        ra = coord.ra.to_string(unit="hourangle", sep=":", precision=3, pad=True)
        dec = coord.dec.to_string(
            unit="deg", sep=":", precision=3, pad=True, alwayssign=True,
        )
        result = single.get_magnitudes(kind=kind)
        instrumental = single.get_magnitudes(kind="instrumental")
        sheet_layout = columns is None and kind == "calibrated"
        headings = (
            ("obs mag", "STD correction", "corrected mag", "error", "Detection σ")
            if sheet_layout else tuple(columns or ("magnitude", "error", "significance"))
        )
        data = []
        for filter_name in filter_names:
            magnitude = _number(result, filter_name, "magnitude")
            error = _number(result, filter_name, "error")
            significance = _number(result, filter_name, "significance")
            if sheet_layout:
                obs = _number(instrumental, filter_name, "magnitude")
                data.append([
                    _text(obs), _text(magnitude - obs, 6), _text(magnitude),
                    _text(error), _text(significance, 2),
                ])
            else:
                values = dict(
                    source_name=info.source_name, ra=ra, dec=dec,
                    field_name=info.observation_field.field_name, filter=filter_name,
                    magnitude_system=result.magnitude_system if result else "",
                    magnitude=_text(magnitude), error=_text(error),
                    significance=_text(significance, 2),
                )
                data.append([values[name] for name in headings])
        row_order = " → ".join(filter_names)
        system = result.magnitude_system if result is not None else "unavailable"
        panels.append(
            '<section style="border:1px solid #aaa;border-radius:6px;'
            'padding:14px;margin:12px 0;max-width:1000px">'
            f'<h3 style="margin:0 0 6px">{escape(info.source_name)}</h3>'
            + _copy_box("Source name", info.source_name, 1)
            + _copy_box(
                "RA / Dec (ICRS)",
                f"{ra} {dec}", 1,
            )
            + f'<p>Field: {escape(info.observation_field.field_name)} · '
            f'System: {escape(system)}</p>'
            f'<p>行顺序 / Row order: {escape(row_order)}</p>'
            + _copy_box("整块复制 / Paste block: " + " | ".join(headings), _tsv(data), len(data))
            + '<details><summary>按列复制 / Copy individual columns</summary>'
            '<div style="display:flex;gap:12px;overflow-x:auto">'
            + "".join(
                _copy_box(name, _tsv([[row[i]] for row in data]), len(data))
                for i, name in enumerate(headings)
            )
            + '</div></details></section>'
        )
    return (
        '<div><p>点击复制按钮，或点击框内后按 Cmd/Ctrl+A、Cmd/Ctrl+C。'
        '整块数据粘贴到 obs mag 的第一行；坐标框粘贴到一个单元格。'
        '复制内容不含标题和 filter 列。</p>'
        + ("".join(panels) if panels else '<p>No sources to display.</p>')
        + '</div>'
    )


def show_google_sheets(sources, kind, filters, columns) -> None:
    from IPython.display import HTML, display

    display(HTML(google_sheets_html(sources, kind, filters, columns)))
