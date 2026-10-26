#!/bin/sh
# Builds manuscript/bci-50ha_data-paper_manuscript.docx from build/manuscript.md
# with pandoc, in the format of the Ecology Data Paper 2.0 manuscript template
# (Arial 11 pt, line spacing 1.15, margins 0.75 in). Highlighted text marks
# what the authors must still supply. Needs pandoc and python3.
#   sh data_paper_and_repo_publication/ecology_data_paper/build/build-manuscript.sh
set -e
cd "$(dirname "$0")"
mkdir -p ../manuscript
pandoc -o pandoc_reference.docx --print-default-data-file reference.docx
python3 patch_reference.py pandoc_reference.docx reference.docx
pandoc manuscript.md -f markdown -t docx --reference-doc=reference.docx -o manuscript_raw.docx
python3 finish_docx.py manuscript_raw.docx ../manuscript/bci-50ha_data-paper_manuscript.docx
rm -f pandoc_reference.docx reference.docx manuscript_raw.docx
