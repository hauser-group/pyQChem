## pyQChem - Input/Output-Tools for Q-Chem

[![DOI](https://zenodo.org/badge/16428728.svg)](https://zenodo.org/badge/latestdoi/16428728)
[![Build Status](https://travis-ci.com/hauser-group/pyQChem.svg?branch=master)](https://travis-ci.com/hauser-group/pyQChem) 

PyQChem is a Python module designed for an intuitive manipulation of [Q-Chem](http://www.q-chem.com) input and output files. It was written with special focus on the features of [IPython](http://ipython.org) such as tab completion and easy access to help docstrings via the question mark operator.

Use conda to add the latest release of pqQChem to your environment (recommended),

```
conda install -c awhauser pyqchem 
```

or clone and install the current Github development version via

```
git clone https://github.com/hauser-group/pyQChem.git
cd pyQChem
pip install .
```

Two IPython notebooks (and support files) are provided as an introduction
to basic inputfile handling and outputfile parsing.
Example Python scripts that use pyQChem can be found in the 'demos' directory.

Andreas W. Hauser, Matthew Goldey, Ehud Tsivion and Michael Wormit, January 2015

Ralf Meyer and Thomas Heavey, September 2018
