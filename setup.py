import os
from setuptools import setup

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "predict_secretome",
    version = "0.1",
    author = "Finlay Maguire",
    author_email = "root@finlaymagui.re",
    description = ("A script generate conservative and permissive secretome "
                   "predictions from a protein fasta")
    license = "MIT",
    keywords = "bioinformatics secretome prediction",
    url = "https://github.com/fmaguire/predict_secretome",
    packages=['eDicer', 'test'],
    long_description=read('README.md'),
)
