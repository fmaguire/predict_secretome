README
======

[![Build Status](https://travis-ci.org/fmaguire/predict_secretome.svg)](https://travis-ci.org/fmaguire/predict_secretome)

This is a package designed to generated secretome predictions by combining 
outputs of various programs (packaged in dependencies). 

Currently it is hardcoded for eukaryotic fungal sequences.

Specifically the program:
* Uses signalp to identify proteins with signal peptides
* Checks the mature sequences of these proteins to ensure there are no TM domains
  (TM domains in the signal peptide itself will pass this stage) via tmhmm
* Identifies signal peptides targeted for secretion using targetp
* Identifies sequences likely to be extracellular compartment using wolfpsort

By default predictions are combined conservatively i.e. the predicted secretome
is only those sequences that fulfill all of the above criteria but the pipe
can optionally be run permissively where passing any stage will be sufficient for output
in the predicted secretome.

The package can also optionally output predicted transporters based on a minimum number of
TM domains appearing in a sequence (or mature sequence in the case of proteins
with signal peptides). This optional output and TM domains threshold can be specified with '-t NUMBER'

A full analysis can be run using the secretome\_pipe.py script (all options 
can be shown using use -h)

Requies:
* linux (likely will work on MACOSX but is not tested)
* python3.4 or python2.7
* perl 5 (specifically installed at /usr/bin/perl)
* biopython 1.64

Run tests by invoking ```python -m unittest``` from main directory
