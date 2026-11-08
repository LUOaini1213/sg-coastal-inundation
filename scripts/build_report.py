"""Render report/report.md from its template (numbers from tools/check_readme_numbers.values)
and, when pandoc and a Chromium browser are available, print it to report/report.pdf.

    python scripts/build_report.py            # report.md + report.pdf (+ docs/report.pdf)
    python scripts/build_report.py --md-only
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import check_readme_numbers as c  # noqa: E402

CSS = """
body{font-family:'Segoe UI',Roboto,Helvetica,Arial,sans-serif;font-size:10.5pt;line-height:1.42;color:#111;
max-width:none;margin:0}
h1{font-size:14pt;margin:1.1em 0 .4em;border-bottom:1px solid #ccc}
header h1.title{font-size:17pt;border:0;margin-bottom:.1em}
p.subtitle{font-size:12pt;color:#444;margin:.1em 0}
p.author,p.date{color:#444;margin:.1em 0}
table{border-collapse:collapse;margin:.6em 0;font-size:9.5pt}
th,td{border-bottom:1px solid #ddd;padding:3px 7px;text-align:left}
th{background:#f3f3f1}
img{display:block;margin:.5em auto}
figure{margin:.4em 0;page-break-inside:avoid}
figcaption{font-size:9pt;color:#555;text-align:center}
@page{size:A4;margin:16mm 15mm}
"""

BROWSERS = [r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
            r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
            "chromium", "google-chrome", "chromium-browser"]


def render_md() -> str:
    V = c.values()
    tpl = (ROOT / "report" / "report.tmpl").read_text(encoding="utf8")
    missing = sorted(set(re.findall(r"\{\{(\w[\w.]*)\}\}", tpl)) - set(V))
    if missing:
        raise SystemExit(f"report template keys without values: {missing}")
    return re.sub(r"\{\{(\w[\w.]*)\}\}", lambda m: V[m.group(1)], tpl)


def main():
    md = render_md()
    out_md = ROOT / "report" / "report.md"
    out_md.write_text(md, encoding="utf8")
    print("wrote", out_md)
    if "--md-only" in sys.argv:
        return
    pandoc = shutil.which("pandoc")
    if not pandoc:
        raise SystemExit("pandoc not found")
    css = ROOT / "report" / "_report.css"
    css.write_text(CSS, encoding="utf8")
    html = ROOT / "report" / "_report.html"
    subprocess.run([pandoc, str(out_md), "-s", "--embed-resources", "--css", str(css), "-o", str(html),
                    "--resource-path", str(ROOT / "report")], check=True, cwd=ROOT / "report")
    browser = next((b for b in BROWSERS if Path(b).exists() or shutil.which(b)), None)
    if browser is None:
        raise SystemExit("no Chromium-based browser for PDF printing")
    pdf = ROOT / "report" / "report.pdf"
    subprocess.run([browser, "--headless=new", "--disable-gpu", "--no-pdf-header-footer",
                    f"--print-to-pdf={pdf}", html.as_uri()], check=True, timeout=180)
    shutil.copy(pdf, ROOT / "docs" / "report.pdf")
    css.unlink()
    html.unlink()
    print("wrote", pdf)


if __name__ == "__main__":
    main()
