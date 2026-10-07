#!/usr/bin/env python3
"""Minimal .xlsx reader: no openpyxl in any env on this machine.

An .xlsx is a zip of XML. This pulls out the shared-string table and walks
each sheet's rows, returning cells as plain strings keyed by column letter.
Enough to read the crosslink supplement; not a general xlsx implementation.

Handles what this file actually uses: shared strings (t="s"), inline strings
(t="inlineStr"), plain numbers, and empty cells. Formulas are read as their
cached value.
"""
import re, zipfile
from xml.etree import ElementTree as ET

NS = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main',
      'r': 'http://schemas.openxmlformats.org/officeDocument/2006/relationships',
      'pr': 'http://schemas.openxmlformats.org/package/2006/relationships'}


def _text(el):
    """All text under an element, joining the <t> runs of a rich string."""
    return ''.join(t.text or '' for t in el.iter(f"{{{NS['m']}}}t"))


def sheet_map(z):
    """Sheet name -> worksheet part name, resolved through the rels file."""
    wb = ET.fromstring(z.read('xl/workbook.xml'))
    rels = ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
    target = {r.get('Id'): r.get('Target') for r in rels}
    out = {}
    for s in wb.iter(f"{{{NS['m']}}}sheet"):
        t = target[s.get(f"{{{NS['r']}}}id")].lstrip('/')
        out[s.get('name')] = t if t.startswith('xl/') else 'xl/' + t
    return out


def shared_strings(z):
    if 'xl/sharedStrings.xml' not in z.namelist():
        return []
    root = ET.fromstring(z.read('xl/sharedStrings.xml'))
    return [_text(si) for si in root.iter(f"{{{NS['m']}}}si")]


def rows(path, sheet_name):
    """Yield {column letter: value} per row, in sheet order."""
    with zipfile.ZipFile(path) as z:
        part = sheet_map(z)[sheet_name]
        sst = shared_strings(z)
        root = ET.fromstring(z.read(part))
        for row in root.iter(f"{{{NS['m']}}}row"):
            cells = {}
            for c in row.iter(f"{{{NS['m']}}}c"):
                ref = c.get('r') or ''
                col = re.match(r'([A-Z]+)', ref)
                if not col:
                    continue
                typ = c.get('t')
                if typ == 'inlineStr':
                    v = _text(c)
                else:
                    node = c.find(f"{{{NS['m']}}}v")
                    v = node.text if node is not None and node.text else ''
                    if typ == 's' and v != '':
                        v = sst[int(v)]
                v = (v or '').strip()
                if v:
                    cells[col.group(1)] = v
            yield cells


if __name__ == '__main__':
    import sys
    for i, r in enumerate(rows(sys.argv[1], sys.argv[2])):
        print(i, r)
        if i > 25:
            break
