#!/bin/sh
# Changes into the test directory and uses unittestdiscover to find all
# unittests. This strange construct is needed to force travis-ci to import
# the package from the virtual-env and not the cloned directory
cd test
python -m unittest discover . 'test_*.py'
