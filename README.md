README
======

This is a package designed to generated secretome predictions by combining 
outputs of various programs (packaged in dependencies). 

There are functions to:
* search for signal peptides in sequences using signalp and then output mature and full sequences
* search for transmembrane domains within mature sequences using tmhmm
* detect signal peptides that target for secretion using targetp
* identify protein sequences that are likely to be extracellular using wolfpsort

A full analysis can be run using the secretome\_pipe.py script (use -h to show
options)

Requies:
* linux 
* python3 
* perl (specifically installed at /usr/bin/perl)
* biopython
