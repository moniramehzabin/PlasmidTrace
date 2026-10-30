"""Export PlasmidTrace evidence tables to a filterable XLSX using the standard library.

No spreadsheet formulas are generated. Source identifiers and evidence strings
remain text; known numeric columns are emitted as numbers.
"""
import csv
import json
from pathlib import Path
import re
import math
import zipfile
from xml.sax.saxutils import escape

NUMERIC = {'length_bp', 'gc_percent', 'coverage_header', 'coverage_ratio', 'rrna_fraction', 'external_graph_neighbors', 'start', 'end'}


def clean(value):
    # XML 1.0 cannot contain most low control characters.
    return re.sub(r'[\x00-\x08\x0b\x0c\x0e-\x1f]', '', str(value))[:32767]


def column(n):
    result = ''
    while n:
        n, remainder = divmod(n-1, 26)
        result = chr(65+remainder) + result
    return result


def sheet_xml(rows):
    count = len(rows[0]) if rows else 1
    parts = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">',
        '<sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/></sheetView></sheetViews>',
        '<cols>']
    widths = []
    for j in range(count):
        header = rows[0][j] if rows else ''
        width = 42 if any(x in header for x in ('detail', 'interpretation', 'evidence', 'pathway_name', 'enzyme_name')) else 24
        widths.append(width)
        parts.append(f'<col min="{j+1}" max="{j+1}" width="{width}" customWidth="1"/>')
    parts.append('</cols><sheetData>')
    for i, row in enumerate(rows, 1):
        lines = max((math.ceil(len(clean(value)) / max(1, widths[j]-2))
                     for j, value in enumerate(row)), default=1)
        height = min(409, max(32, lines * 16 + 10))
        parts.append(f'<row r="{i}" ht="{height}" customHeight="1">')
        for j, raw in enumerate(row, 1):
            text = clean(raw); address = f'{column(j)}{i}'
            if i == 1: text = text.replace('_', ' ')
            numeric = i > 1 and rows[0][j-1] in NUMERIC and re.fullmatch(r'-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?', text)
            style = 1 if i==1 else (3 if numeric and rows[0][j-1]=='evalue' else 2)
            if numeric:
                parts.append(f'<c r="{address}" s="{style}"><v>{text}</v></c>')
            else:
                parts.append(f'<c r="{address}" s="{style}" t="inlineStr"><is><t xml:space="preserve">{escape(text)}</t></is></c>')
        parts.append('</row>')
    parts.append('</sheetData>')
    if rows:
        parts.append(f'<autoFilter ref="A1:{column(count)}{len(rows)}"/>')
    parts.append('</worksheet>')
    return ''.join(parts)


def export_results(folder):
    folder = Path(folder)
    sheets = [('Read me', [
        ['Item', 'Explanation'],
        ['Source', 'TSV files in output/ are the authoritative full outputs.'],
        ['Interpretation', 'Evidence for review, not confirmed plasmid identity or gene activity.'],
        ['Missing checks', 'Unassessed checks do not establish absence.'],
        ['Long fields', 'Cells are limited to 32767 characters; consult original TSV for complete evidence.'],
    ])]
    assessment_path = folder / 'assessment.json'
    if assessment_path.is_file():
        assessment = json.loads(assessment_path.read_text(encoding='utf-8'))
        sheets[0][1].append(['Run assessment', assessment.get('label', '')])
        for warning in assessment.get('warnings', []):
            sheets[0][1].append(['Limitation', warning])
    used = {'Read me'}
    for path in sorted((folder / 'output').rglob('*.tsv')):
        with path.open(encoding='utf-8-sig', newline='') as handle:
            rows = list(csv.reader(handle, delimiter='\t'))
        if not rows: continue
        width = max(map(len, rows))
        rows = [r + [''] * (width - len(r)) for r in rows]
        base = re.sub(r'[\[\]:*?/\\]', '_', str(path.relative_to(folder / 'output')).removesuffix('.tsv'))[:27]
        name = base
        n = 2
        while name in used:
            name = base[:26] + '_' + str(n); n += 1
        used.add(name)
        sheets.append((name, rows))
    output = folder / 'PlasmidTrace_results.xlsx'
    with zipfile.ZipFile(output, 'w', zipfile.ZIP_DEFLATED) as z:
        content = ['<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/><Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>']
        for i in range(1, len(sheets)+1):
            content.append(f'<Override PartName="/xl/worksheets/sheet{i}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>')
        z.writestr('[Content_Types].xml', ''.join(content)+'</Types>')
        z.writestr('_rels/.rels', '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>')
        workbook = ['<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>']
        relationships = ['<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">']
        for i, (name, rows) in enumerate(sheets, 1):
            workbook.append(f'<sheet name="{escape(name)}" sheetId="{i}" r:id="rId{i}"/>')
            relationships.append(f'<Relationship Id="rId{i}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{i}.xml"/>')
            z.writestr(f'xl/worksheets/sheet{i}.xml', sheet_xml(rows))
        relationships.append(f'<Relationship Id="rId{len(sheets)+1}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/></Relationships>')
        z.writestr('xl/workbook.xml', ''.join(workbook)+'</sheets></workbook>')
        z.writestr('xl/_rels/workbook.xml.rels', ''.join(relationships))
        z.writestr('xl/styles.xml', '''<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
<fonts count="2"><font><sz val="11"/><name val="Calibri"/></font><font><b/><sz val="11"/><color rgb="FFFFFFFF"/><name val="Calibri"/></font></fonts>
<fills count="3"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill><fill><patternFill patternType="solid"><fgColor rgb="FF1D654F"/><bgColor indexed="64"/></patternFill></fill></fills>
<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>
<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>
<cellXfs count="4"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>
<xf numFmtId="0" fontId="1" fillId="2" borderId="0" xfId="0" applyAlignment="1"><alignment vertical="center" wrapText="1"/></xf>
<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0" applyAlignment="1"><alignment vertical="top" wrapText="1"/></xf>
<xf numFmtId="11" fontId="0" fillId="0" borderId="0" xfId="0" applyAlignment="1"><alignment vertical="top" wrapText="1"/></xf></cellXfs>
<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles></styleSheet>''')
    return output
