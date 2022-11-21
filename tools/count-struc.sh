#!/bin/bash
set -u -e
python3 -c 'import sys, numpy as np; print(len(np.load(sys.argv[1])))' $1