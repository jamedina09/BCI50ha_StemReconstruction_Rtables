"""Patch pandoc's default reference.docx: Arial 11 pt like the original document,
black headings, compact table text and a ruled table style."""
import re
import shutil
import sys
import zipfile

src, dst = sys.argv[1], sys.argv[2]

with zipfile.ZipFile(src) as z:
    files = {n: z.read(n) for n in z.namelist()}

s = files["word/styles.xml"].decode("utf8")
ARIAL = '<w:rFonts w:ascii="Arial" w:hAnsi="Arial" w:eastAsia="Arial" w:cs="Arial" />'


def style_block(sid):
    m = re.search(r'<w:style [^>]*w:styleId="%s".*?</w:style>' % sid, s, re.S)
    if not m:
        raise SystemExit("style not found: " + sid)
    return m


def replace_style(sid, fn):
    global s
    m = style_block(sid)
    s = s[: m.start()] + fn(m.group(0)) + s[m.end():]


# Document defaults: Arial 11 pt
s = re.sub(r'<w:rPrDefault>\s*<w:rPr>\s*<w:rFonts[^>]*/>\s*<w:sz w:val="24" />\s*<w:szCs w:val="24" />',
           '<w:rPrDefault><w:rPr>' + ARIAL + '<w:sz w:val="22" /><w:szCs w:val="22" />', s)
s = s.replace('<w:pPrDefault>\n      <w:pPr>\n        <w:spacing w:after="200" />',
              '<w:pPrDefault>\n      <w:pPr>\n        <w:spacing w:after="120" w:line="276" w:lineRule="auto" />')


def headings(block, size):
    block = re.sub(r'<w:rFonts[^>]*/>', ARIAL, block, flags=re.S)
    block = re.sub(r'<w:color [^>]*/>', '<w:color w:val="000000" />', block, flags=re.S)
    block = re.sub(r'<w:sz w:val="\d+" />', '<w:sz w:val="%d" />' % size, block)
    block = re.sub(r'<w:szCs w:val="\d+" />', '<w:szCs w:val="%d" />' % size, block)
    if '<w:b />' not in block and '<w:b/>' not in block:
        block = block.replace('<w:rPr>', '<w:rPr><w:b /><w:bCs />', 1)
    return block


for sid, size in [("Heading1", 32), ("Heading2", 26), ("Heading3", 22), ("Heading4", 22)]:
    replace_style(sid, lambda b, size=size: headings(b, size))
replace_style("Title", lambda b: headings(b, 36).replace('<w:jc w:val="center" />', '<w:jc w:val="left" />'))
replace_style("Subtitle", lambda b: headings(b, 24).replace('<w:rPr><w:b /><w:bCs />', '<w:rPr>', 1))

# Body text: no space before, 6 pt after
replace_style("BodyText", lambda b: b.replace('<w:spacing w:before="180" w:after="180" />',
                                              '<w:spacing w:before="0" w:after="120" />'))
# Normal: Arial 11 pt set on the style too (viewers that ignore docDefaults)
replace_style("Normal", lambda b: b.replace(
    '<w:qFormat />', '<w:qFormat />\n    <w:rPr>' + ARIAL + '<w:sz w:val="22" /><w:szCs w:val="22" /></w:rPr>', 1))
# Tight list items keep the body size, with little space between items
replace_style("Compact", lambda b: b.replace(
    '<w:spacing w:before="36" w:after="36" />', '<w:spacing w:before="0" w:after="40" />'))
# Table text (applied to table cells after conversion): 9 pt, single spacing
s = s.replace('</w:styles>', '''<w:style w:type="paragraph" w:customStyle="1" w:styleId="TableText">
    <w:name w:val="Table Text" />
    <w:basedOn w:val="Normal" />
    <w:qFormat />
    <w:pPr><w:spacing w:before="20" w:after="20" w:line="240" w:lineRule="auto" /></w:pPr>
    <w:rPr><w:sz w:val="18" /><w:szCs w:val="18" /></w:rPr>
  </w:style>
</w:styles>''')
# Column and file names: Courier New 10 pt
replace_style("VerbatimChar", lambda b: re.sub(
    r'<w:rFonts[^>]*/>\s*<w:sz w:val="\d+" />',
    '<w:rFonts w:ascii="Courier New" w:hAnsi="Courier New" w:cs="Courier New" /><w:sz w:val="20" />', b))
replace_style("Hyperlink", lambda b: re.sub(r'<w:color [^>]*/>', '<w:color w:val="1155CC" />', b))
replace_style("Caption", lambda b: b.replace('<w:i />', '<w:sz w:val="20" /><w:szCs w:val="20" />'))

# Ruled tables: top and bottom rules, light inner horizontal lines, shaded bold header
TABLE = '''<w:style w:type="table" w:default="1" w:styleId="Table">
    <w:name w:val="Table" />
    <w:basedOn w:val="TableNormal" />
    <w:qFormat />
    <w:tblPr>
      <w:tblInd w:w="0" w:type="dxa" />
      <w:tblBorders>
        <w:top w:val="single" w:sz="8" w:space="0" w:color="000000" />
        <w:bottom w:val="single" w:sz="8" w:space="0" w:color="000000" />
        <w:insideH w:val="single" w:sz="4" w:space="0" w:color="BFBFBF" />
      </w:tblBorders>
      <w:tblCellMar>
        <w:top w:w="40" w:type="dxa" />
        <w:left w:w="80" w:type="dxa" />
        <w:bottom w:w="40" w:type="dxa" />
        <w:right w:w="80" w:type="dxa" />
      </w:tblCellMar>
    </w:tblPr>
    <w:tblStylePr w:type="firstRow">
      <w:rPr><w:b /><w:bCs /></w:rPr>
      <w:tcPr>
        <w:tcBorders>
          <w:bottom w:val="single" w:sz="8" w:space="0" w:color="000000" />
        </w:tcBorders>
        <w:shd w:val="clear" w:color="auto" w:fill="F2F2F2" />
        <w:vAlign w:val="bottom" />
      </w:tcPr>
    </w:tblStylePr>
  </w:style>'''
replace_style("Table", lambda b: TABLE)

files["word/styles.xml"] = s.encode("utf8")
with zipfile.ZipFile(dst, "w", zipfile.ZIP_DEFLATED) as z:
    for n, data in files.items():
        z.writestr(n, data)
print("reference written:", dst)
