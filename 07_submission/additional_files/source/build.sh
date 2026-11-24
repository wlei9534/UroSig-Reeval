#!/bin/bash
# usage: build.sh <name-without-extension> "<title>" "<subtitle>"
set -e
pandoc "$1.md" -t latex --no-highlight -o body.tex
cat > "$1.tex" <<TEX
\documentclass[10pt,a4paper]{article}
\usepackage[margin=20mm]{geometry}
\usepackage{fontspec}\setmainfont{Liberation Sans}\setmonofont[Scale=0.9]{Liberation Mono}
\usepackage{longtable,booktabs,array,calc}
\usepackage[hidelinks]{hyperref}
\usepackage{xurl}
\providecommand{\tightlist}{\setlength{\itemsep}{2pt}\setlength{\parskip}{0pt}}
\providecommand{\pandocbounded}[1]{#1}
\newcommand{\passthrough}[1]{#1}
\input{header.tex}
\setcounter{secnumdepth}{0}
\begin{document}
{\LARGE\bfseries $2\par}\vspace{4pt}{\large $3\par}\vspace{10pt}
\input{body.tex}
\end{document}
TEX
xelatex -interaction=nonstopmode "$1.tex" > build.log 2>&1 || { grep -A8 "^!" build.log | head -30; exit 1; }
xelatex -interaction=nonstopmode "$1.tex" > build.log 2>&1
grep -c "Overfull" build.log || true
