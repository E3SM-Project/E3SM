This is the top-level of the SCREAM design documentation. That documentation
is written in LaTeX and is spread over a directory structure mimicking that 
of the actual code within the ../src directory. 

To compile the documentation for a particular process, go to the docs 
directory for that process and issue:

pdflatex main.tex
bibtex main
pdflatex main.tex

If you want to compile the documentation for the whole model, issue the same
command from the top-level docs directory.

Obviously, you need to have a working Latex installation for this to work.