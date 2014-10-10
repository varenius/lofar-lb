#!/bin/bash
pdflatex main.tex;
#bibtex main;
pdflatex main.tex;
pdflatex main.tex;
# TO SPELLCHECK: aspell --lang=en_GB -c lbchapter.tex 
