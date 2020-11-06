#!/bin/bash

VERSIONS="python3.8 python3.7 python3.6"

for PYVER in $VERSIONS; do
  $PYVER -m pip install -r dev-requirements.txt
  $PYVER setup.py bdist_wheel
done

