#!/bin/bash


jupyter nbconvert --to script post-01-consolidate.ipynb
python3 post-01-consolidate.py
rm post-01-consolidate.py
