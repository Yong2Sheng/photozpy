"""Check the actual cell payloads, rather than the appearance of TSV text."""

import csv
from html.parser import HTMLParser
from io import StringIO

from astropy.coordinates import SkyCoord
import pytest

from photozpy.sources import FilterMagnitude, ObservationField, SourceInfo, Sources
from photozpy.telescope import SARA_RM


class CopyBoxes(HTMLParser):
    def __init__(self, html):
        super().__init__()
        self.boxes = []
        self.current = None
        self.feed(html)

    def handle_starttag(self, tag, attrs):
        if tag == "textarea":
            self.current = [dict(attrs)["aria-label"], ""]

    def handle_data(self, data):
        if self.current is not None:
            self.current[1] += data

    def handle_endtag(self, tag):
        if tag == "textarea":
            self.boxes.append(self.current)
            self.current = None


@pytest.mark.parametrize("declination, dec_text", [
    (-2, "-02:00:00.000"), (2, "+02:00:00.000"),
])
def test_google_sheets_pastes_into_existing_four_by_five_block(
    monkeypatch, declination, dec_text,
):
    coordinate = SkyCoord(150, declination, unit="deg")
    field = ObservationField(
        field_name="field", telescope=SARA_RM, nominal_pointing=coordinate,
    )
    sources = Sources(SourceInfo(
        source_name='target <script>alert("x")</script>', source_role="science",
        source_coordinate=coordinate, observation_field=field,
    ))
    for kind, values, system in (
        ("instrumental", (-28.2, -27.9, -27.4, -26.3), "instrumental"),
        ("calibrated", (19.7, 19.2, 18.5, 18.8), "AB"),
    ):
        sources.set_magnitudes(kind=kind, magnitudes=FilterMagnitude(
            telescope=SARA_RM, filter_names=tuple(f"SDSS_{b}" for b in "griz"),
            magnitudes=values, errors=(.08, .09, .07, -99),
            significance=(12, 13, 14, 2.6), magnitude_system=system,
        ))
    displayed = []
    monkeypatch.setattr("IPython.display.display", displayed.append)
    sources.show_magnitudes(format="google_sheets")
    html = displayed[0].data
    assert '<script>alert' not in html
    boxes = CopyBoxes(html).boxes
    assert boxes.pop(0) == ["Source name", sources.source_info.source_name]
    assert boxes[0][1] == f"10:00:00.000 {dec_text}"
    assert "\t" not in boxes[0][1] and "\n" not in boxes[0][1]
    assert list(csv.reader(StringIO(boxes[0][1]), delimiter="\t")) == [
        [f"10:00:00.000 {dec_text}"],
    ]
    rows = list(csv.reader(StringIO(boxes[1][1]), delimiter="\t"))
    assert len(rows) == 4 and all(len(row) == 5 for row in rows)
    assert rows[0] == ["-28.2000", "47.900000", "19.7000", "0.0800", "12.00"]
    assert rows[3][-2:] == ["-99.0000", "2.60"]
    # The single-column payload has no coordinates, source names, or headers.
    assert boxes[2][1] == "-28.2000\n-27.9000\n-27.4000\n-26.3000\n"
    assert sources.get_magnitudes(kind="calibrated").magnitudes[1].value == 19.7

    sources.show_magnitudes(format="google_sheets", columns=("ra", "dec"))
    coordinate_rows = list(csv.reader(
        StringIO(CopyBoxes(displayed[-1].data).boxes[2][1]), delimiter="\t",
    ))
    assert coordinate_rows == [["10:00:00.000", dec_text]] * 4


def test_missing_results_keep_four_rows_and_explicit_column_order(monkeypatch):
    coordinate = SkyCoord(150, 2, unit="deg")
    field = ObservationField(
        field_name="field", telescope=SARA_RM, nominal_pointing=coordinate,
    )
    sources = Sources(SourceInfo(
        source_name="missing", source_role="science",
        source_coordinate=coordinate, observation_field=field,
    ))
    displayed = []
    monkeypatch.setattr("IPython.display.display", displayed.append)
    sources.show_magnitudes(
        format="google_sheets", columns=("magnitude", "error"),
    )
    boxes = CopyBoxes(displayed[0].data).boxes
    assert boxes.pop(0) == ["Source name", "missing"]
    assert boxes[1][1] == "\t\n" * 4
    assert boxes[2][1] == "\n" * 4
