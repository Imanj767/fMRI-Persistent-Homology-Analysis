# fMRI Persistent Homology Analysis

## Overview

This Brainlife app computes persistent homology features from functional MRI (fMRI) time series data. Persistent homology is a topological data analysis (TDA) technique that captures the shape and connectivity patterns in high-dimensional data, making it valuable for analyzing brain functional connectivity.

## Features

- **Topological Data Analysis**: Computes persistent homology up to dimension 2
- **Multiple Distance Metrics**: Supports correlation, Euclidean, and cosine distances
- **ROI Extraction Methods**: PCA, spatial parcellation, or direct voxel selection
- **Sliding Window Analysis**: Temporal analysis of topological features
- **Comprehensive Output**: Persistence diagrams, numerical features, and visualizations

## Inputs

### Required
- **func**: Functional MRI data in NIfTI format (4D: x, y, z, time)

### Optional
- **mask**: Brain mask for region-of-interest analysis

## Outputs

- **persistence_diagrams.npy**: Raw persistence diagrams for each homology dimension
- **persistence_diagrams.png**: Visualization of persistence diagrams
- **topological_features.csv**: Numerical features extracted from global analysis
- **window_features.csv**: Time-resolved topological features from sliding windows
- **feature_evolution.png**: Temporal evolution of topological features
- **analysis_summary.json**: Configuration and processing summary

## Configuration Parameters

- **max_dimension** (default: 2): Maximum homology dimension to compute
- **window_size** (default: 50): Number of timepoints in each sliding window
- **overlap** (default: 0.5): Overlap ratio between consecutive windows
- **distance_metric** (default: "correlation"): Distance metric for persistence computation
  - Options: "correlation", "euclidean", "cosine"
- **roi_extraction** (default: "parcellation"): Method for extracting ROI time series
  - Options: "parcellation", "voxel", "pca"

## Algorithm Details

### 1. Data Preprocessing
- Load fMRI data and apply brain mask if provided
- Extract ROI time series using selected method:
  - **Parcellation**: Divide brain into spatial regions and average within regions
  - **PCA**: Use principal components to reduce dimensionality
  - **Voxel**: Select subset of voxels directly

### 2. Distance Matrix Computation
- **Correlation Distance**: 1 - |correlation coefficient|
- **Euclidean Distance**: Normalized Euclidean distance between time series
- **Cosine Distance**: Angular distance between time series vectors

### 3. Persistent Homology Computation
- Uses Ripser algorithm for efficient persistence computation
- Computes homology in dimensions 0, 1, and optionally 2
- Generates persistence diagrams showing birth and death of topological features

### 4. Feature Extraction
For each homology dimension:
- Number of topological features
- Total persistence (sum of lifetimes)
- Maximum persistence
- Mean and standard deviation of persistence
- Mean birth and death times

### 5. Sliding Window Analysis
- Creates overlapping windows across time
- Computes topological features for each window
- Tracks temporal evolution of brain topology

## Interpretation

### Homology Dimensions
- **Dimension 0**: Connected components (functional modules)
- **Dimension 1**: Loops/cycles in functional connectivity
- **Dimension 2**: Voids/cavities in connectivity structure

### Key Features
- **Number of Features**: Complexity of functional organization
- **Persistence**: Stability/significance of topological structures
- **Temporal Evolution**: Changes in brain topology over time

## Installation and Usage

### For Brainlife Platform
This app is designed to run on the Brainlife platform. Simply:
1. Upload your fMRI data (and optional mask)
2. Configure parameters
3. Submit the job

### Local Development
```bash
# Clone repository
git clone https://github.com/yourusername/app-fmri-persistent-homology
cd app-fmri-persistent-homology

# Install dependencies
pip install -r requirements.txt

# Run analysis
python main.py
```

### Docker Usage
```bash
# Build image
docker build -t fmri-persistent-homology .

# Run container
docker run -v /path/to/data:/app/data fmri-persistent-homology
```

## Dependencies

- Python 3.9+
- NumPy, SciPy, scikit-learn
- NiBabel (neuroimaging data)
- Matplotlib, Seaborn (visualization)
- GUDHI, Ripser (persistent homology)
- Persim (persistence diagram tools)

## References

1. Edelsbrunner, H., & Harer, J. (2010). Computational topology: an introduction. American Mathematical Society.
2. Saggar, M., et al. (2018). Towards a new approach to reveal dynamical organization of the brain using topological data analysis. Nature Communications, 9(1), 1399.
3. Tralie, C., et al. (2018). Ripser.py: A lean persistent homology library for Python. Journal of Open Source Software, 3(29), 925.

## Support

For issues and questions:
- GitHub Issues: [repository-url]/issues
- Brainlife Forum: https://brainlife.io/forum
- Email: [your-email]

## License

MIT License - see LICENSE file for details.

## Citation

If you use this app in your research, please cite:
```
[Your Name] (2024). fMRI Persistent Homology Analysis. 
Brainlife App. doi: [DOI if available]
```
