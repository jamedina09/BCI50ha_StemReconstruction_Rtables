"""Finish the pandoc output: table cells take the TableText style (9 pt)
instead of pandoc's Compact (also used for tight lists), and the page margins
are those of the Ecology Data Paper 2.0 manuscript template (0.75 in)."""
import re
import sys
import zipfile

src, dst = sys.argv[1], sys.argv[2]
with zipfile.ZipFile(src) as z:
    files = {n: z.read(n) for n in z.namelist()}

doc = files["word/document.xml"].decode("utf8")
doc = re.sub(r"<w:tbl>.*?</w:tbl>",
             lambda m: m.group(0).replace('<w:pStyle w:val="Compact" />', '<w:pStyle w:val="TableText" />'),
             doc, flags=re.S)
# Page of the template: US Letter, margins 0.75 in (pgSz and pgMar follow
# footnotePr in the section properties).
doc = re.sub(r"<w:pgSz [^>]*/>|<w:pgMar [^>]*/>", "", doc)
if doc.count("</w:sectPr>") != 1:
    raise SystemExit("expected one section, found %d" % doc.count("</w:sectPr>"))
doc = doc.replace(
    "</w:sectPr>",
    '<w:pgSz w:w="12240" w:h="15840" />'
    '<w:pgMar w:top="1080" w:right="1080" w:bottom="1080" w:left="1080" w:header="720" w:footer="720" w:gutter="0" />'
    "</w:sectPr>")
files["word/document.xml"] = doc.encode("utf8")

with zipfile.ZipFile(dst, "w", zipfile.ZIP_DEFLATED) as z:
    for name, data in files.items():
        z.writestr(name, data)
print("written:", dst)
