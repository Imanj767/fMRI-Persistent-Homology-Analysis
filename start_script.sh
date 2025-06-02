#!/bin/bash

# Brainlife App: fMRI Persistent Homology Analysis
# This script processes fMRI data and computes persistent homology features

set -e
set -x

# Display system information
echo "Starting fMRI Persistent Homology Analysis"
echo "Python version: $(python --version)"
echo "Working directory: $(pwd)"
echo "Available files:"
ls -la

# Check if required input files exist
if [ ! -f "func/func.nii.gz" ] && [ ! -f "func.nii.gz" ]; then
    echo "Error: No fMRI functional data found"
    echo "Expected: func/func.nii.gz or func.nii.gz"
    exit 1
fi

# Check for optional mask
if [ -f "mask/mask.nii.gz" ]; then
    echo "Brain mask found: mask/mask.nii.gz"
else
    echo "No brain mask provided (optional)"
fi

# Create output directory
mkdir -p output

# Run the main processing script
echo "Starting persistent homology computation..."
python main.py

# Check if output was generated
if [ -f "output/analysis_summary.json" ]; then
    echo "Analysis completed successfully"
    echo "Output files:"
    ls -la output/
else
    echo "Error: Analysis did not complete successfully"
    exit 1
fi

echo "fMRI Persistent Homology Analysis finished"
