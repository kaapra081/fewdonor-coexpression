"""Split the built paper into main_paper.pdf (title, body, references) and supplementary.pdf (cover with contents, appendices, reference list).
python make_supplement.py <build dir with main.pdf, main.aux, main.bbl> <tectonic> <bundle> <out dir>"""
import os, sys, re, subprocess, pypdfium2 as pdfium
B, TEC, BUNDLE, OUTD = sys.argv[1:5]; W = os.path.join(B, "_supp"); os.makedirs(W, exist_ok=True)
src = pdfium.PdfDocument(os.path.join(B, "main.pdf")); n = len(src)
fa = next(i for i in range(n) if src[i].get_textpage().get_text_range().startswith("A Proofs"))
aux = open(os.path.join(B, "main.aux")).read()
secs = re.findall(r"\\contentsline \{section\}\{\\numberline \{([A-Z])\}(.*?)\}\{(\d+)\}", aux)
rows = "\n".join("%s & %s & %s\\\\" % (L, re.sub(r"\\label\{[^}]*\}", "", t), p) for L, t, p in secs)
tex = r"""\documentclass[10pt,letterpaper]{article}
\usepackage[margin=1in]{geometry}
\usepackage[numbers,sort&compress]{natbib}
\usepackage[hidelinks]{hyperref}
\usepackage{xurl}
\usepackage{microtype}
\begin{document}
\thispagestyle{empty}
\begin{center}
{\Large\bfseries Supplementary Material}\\[10pt]
{\large Why pseudobulk gene co-expression fails in large single-cell cohorts,\\ and how to fix it}\\[10pt]
Kaavin Prasanna
\end{center}
\vspace{12pt}
\noindent This file contains Appendices A--LAST of the paper. Section, theorem, table, figure and equation numbers, and the citation numbers [$n$], are those of the main paper; its reference list is repeated at the end of this file. Page numbers continue those of the main paper.

\vspace{8pt}
\noindent Code, frozen analysis plans and result files: \url{https://github.com/kaapra081/fewdonor-coexpression}

\vspace{14pt}
\noindent\begin{tabular}{@{}lp{0.80\linewidth}r@{}}
\multicolumn{3}{@{}l}{\bfseries Contents}\\[2pt]
ROWS
\end{tabular}
\clearpage
\setcounter{page}{PAGE}
\renewcommand{\refname}{References (as in the main paper)}
\input{main.bbl}
\end{document}
""".replace("LAST", secs[-1][0]).replace("ROWS", rows).replace("PAGE", str(n))
open(os.path.join(W, "wrap.tex"), "w").write(tex); subprocess.run(["cp", os.path.join(B, "main.bbl"), W], check=True)
r = subprocess.run([TEC, "-X", "compile", "--bundle", BUNDLE, "wrap.tex"], cwd=W, capture_output=True, text=True); assert r.returncode == 0, r.stderr[-2000:]
wrap = pdfium.PdfDocument(os.path.join(W, "wrap.pdf"))
m = pdfium.PdfDocument.new(); m.import_pages(src, list(range(0, fa))); m.save(os.path.join(OUTD, "main_paper.pdf"))
s = pdfium.PdfDocument.new(); s.import_pages(wrap, [0]); s.import_pages(src, list(range(fa, n))); s.import_pages(wrap, list(range(1, len(wrap)))); s.save(os.path.join(OUTD, "supplementary.pdf"))
print("main_paper", fa, "pages; supplementary", 1 + n - fa + len(wrap) - 1, "pages; appendices", [L for L, _, _ in secs])
