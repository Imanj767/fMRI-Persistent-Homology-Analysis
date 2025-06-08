#!/bin/bash

# fMRI Persistent Homology Analysis App for Brainlife
# This script runs the Python analysis using Singularity container

set -e
set -x

# Run the Python analysis script using Singularity
singularity exec -e docker://brainlife/ga-python:lab328-dipy141-pybrainlife-1.0 python ./app.py
