#!/usr/bin/env python3

"""
Brainlife App: fMRI Persistent Homology Analysis
This app computes persistent homology from fMRI time series data
"""

import json
import numpy as np
import nibabel as nib
import os
import sys
from sklearn.metrics import pairwise_distances
import matplotlib.pyplot as plt
import pandas as pd

# Try to import ripser for persistent homology
try:
    from ripser import ripser
    HAS_RIPSER = True
except ImportError:
    print("Warning: ripser not available, using alternative method")
    HAS_RIPSER = False

# Try to import nilearn components
try:
    from nilearn import image
    from nilearn.regions import RegionExtractor
    from nilearn.connectome import ConnectivityMeasure
    from nilearn.maskers import NiftiLabelsMasker
    HAS_NILEARN = True
except ImportError:
    print("Warning: some nilearn components not available")
    HAS_NILEARN = False

def load_config():
    """Load configuration from config.json"""
    try:
        with open('config.json', 'r') as f:
            config = json.load(f)
        return config
    except FileNotFoundError:
        print("Error: config.json file not found")
        sys.exit(1)
    except json.JSONDecodeError:
        print("Error: Invalid JSON in config.json")
        sys.exit(1)

def load_fmri_data(fmri_path):
    """Load and preprocess fMRI data"""
    print(f"Loading fMRI data from: {fmri_path}")
    
    if not os.path.exists(fmri_path):
        raise FileNotFoundError(f"fMRI file not found: {fmri_path}")
    
    # Load fMRI image
    fmri_img = nib.load(fmri_path)
    print(f"fMRI shape: {fmri_img.shape}")
    
    # Basic validation
    if len(fmri_img.shape) != 4:
        raise ValueError(f"Expected 4D fMRI data, got {len(fmri_img.shape)}D")
    
    return fmri_img

def extract_time_series_simple(fmri_img, n_regions=100):
    """Simple time series extraction without complex dependencies"""
    print(f"Extracting time series using simple method with {n_regions} regions")
    
    # Get fMRI data
    fmri_data = fmri_img.get_fdata()
    
    # Create a simple brain mask (non-zero voxels)
    brain_mask = np.any(fmri_data != 0, axis=3)
    
    # Get brain voxel coordinates
    brain_coords = np.where(brain_mask)
    n_voxels = len(brain_coords[0])
    
    if n_voxels == 0:
        raise ValueError("No brain voxels found")
    
    print(f"Found {n_voxels} brain voxels")
    
    # Sample regions if we have too many voxels
    if n_voxels > n_regions:
        # Randomly sample voxels for regions
        indices = np.random.choice(n_voxels, n_regions, replace=False)
        selected_coords = (
            brain_coords[0][indices],
            brain_coords[1][indices], 
            brain_coords[2][indices]
        )
    else:
        selected_coords = brain_coords
        n_regions = n_voxels
    
    # Extract time series
    time_series = []
    for i in range(len(selected_coords[0])):
        x, y, z = selected_coords[0][i], selected_coords[1][i], selected_coords[2][i]
        ts = fmri_data[x, y, z, :]
        
        # Simple preprocessing: demean and standardize
        ts = ts - np.mean(ts)
        if np.std(ts) > 0:
            ts = ts / np.std(ts)
        
        time_series.append(ts)
    
    time_series = np.array(time_series).T  # Shape: (timepoints, regions)
    
    print(f"Extracted time series shape: {time_series.shape}")
    return time_series

def extract_time_series_advanced(fmri_img, atlas_path=None, n_regions=100):
    """Advanced time series extraction using nilearn (if available)"""
    
    if not HAS_NILEARN:
        return extract_time_series_simple(fmri_img, n_regions)
    
    print("Using advanced nilearn-based extraction")
    
    if atlas_path and os.path.exists(atlas_path):
        # Use provided atlas
        atlas_img = nib.load(atlas_path)
        print(f"Using atlas: {atlas_path}")
        
        # Extract time series using atlas
        masker = NiftiLabelsMasker(
            labels_img=atlas_img,
            standardize=True,
            detrend=True,
            low_pass=0.1,
            high_pass=0.01,
            t_r=2.0
        )
        time_series = masker.fit_transform(fmri_img)
        
    else:
        # Use region extraction for creating parcellation
        print(f"Creating {n_regions} regions using RegionExtractor")
        
        try:
            extractor = RegionExtractor(
                fmri_img,
                n_regions=n_regions,
                clustering_method='hierarchical',
                standardize=True,
                detrend=True,
                low_pass=0.1,
                high_pass=0.01,
                t_r=2.0
            )
            
            time_series = extractor.fit_transform(fmri_img)
        except Exception as e:
            print(f"RegionExtractor failed: {e}, falling back to simple method")
            return extract_time_series_simple(fmri_img, n_regions)
    
    print(f"Extracted time series shape: {time_series.shape}")
    return time_series

def compute_functional_connectivity(time_series, method='correlation'):
    """Compute functional connectivity matrix"""
    print(f"Computing functional connectivity using {method}")
    
    if HAS_NILEARN:
        try:
            connectivity_measure = ConnectivityMeasure(kind=method)
            connectivity_matrix = connectivity_measure.fit_transform([time_series])[0]
        except:
            # Fallback to simple correlation
            connectivity_matrix = np.corrcoef(time_series.T)
    else:
        # Simple correlation calculation
        connectivity_matrix = np.corrcoef(time_series.T)
    
    print(f"Connectivity matrix shape: {connectivity_matrix.shape}")
    return connectivity_matrix

def simple_persistent_homology(distance_matrix, max_edge_length=1.0):
    """Simple implementation of persistent homology (0-dimensional only)"""
    print("Computing persistent homology using simple method...")
    
    n = distance_matrix.shape[0]
    
    # Create list of edges with distances
    edges = []
    for i in range(n):
        for j in range(i+1, n):
            if distance_matrix[i, j] <= max_edge_length:
                edges.append((distance_matrix[i, j], i, j))
    
    # Sort edges by distance
    edges.sort()
    
    # Union-Find for connected components
    parent = list(range(n))
    
    def find(x):
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]
    
    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py
            return True
        return False
    
    # Track connected components (0-dimensional homology)
    births = [0.0] * n  # Each point is born at distance 0
    deaths = []
    
    # Process edges
    for dist, i, j in edges:
        if union(i, j):
            # A component died (merged with another)
            deaths.append(dist)
    
    # Remaining components live forever
    n_remaining = len(set(find(i) for i in range(n)))
    deaths.extend([np.inf] * n_remaining)
    
    # Create persistence diagram for dimension 0
    persistence_0 = np.column_stack([births, deaths])
    
    return [persistence_0]  # Only 0-dimensional homology

def compute_persistent_homology(distance_matrix, maxdim=1):
    """Compute persistent homology using Ripser or simple method"""
    
    if HAS_RIPSER:
        print("Computing persistent homology using Ripser...")
        try:
            diagrams = ripser(distance_matrix, distance_matrix=True, maxdim=maxdim)
            persistence_diagrams = diagrams['dgms']
            return persistence_diagrams, diagrams
        except Exception as e:
            print(f"Ripser failed: {e}, using simple method")
    
    # Fallback to simple method
    persistence_diagrams = simple_persistent_homology(distance_matrix)
    return persistence_diagrams, {'dgms': persistence_diagrams}

def compute_persistence_features(persistence_diagrams):
    """Extract features from persistence diagrams"""
    features = {}
    
    for dim, diagram in enumerate(persistence_diagrams):
        if len(diagram) > 0:
            # Birth and death times
            births = diagram[:, 0]
            deaths = diagram[:, 1]
            
            # Persistence (lifetime)
            persistence = deaths - births
            
            # Remove infinite persistence points for dimension > 0
            if dim > 0:
                finite_mask = np.isfinite(persistence)
                persistence = persistence[finite_mask]
                births = births[finite_mask]
                deaths = deaths[finite_mask]
            
            if len(persistence) > 0:
                features[f'dim_{dim}_num_features'] = len(persistence)
                features[f'dim_{dim}_max_persistence'] = float(np.max(persistence))
                features[f'dim_{dim}_mean_persistence'] = float(np.mean(persistence))
                features[f'dim_{dim}_std_persistence'] = float(np.std(persistence))
                features[f'dim_{dim}_sum_persistence'] = float(np.sum(persistence))
                
                # Birth and death statistics
                features[f'dim_{dim}_mean_birth'] = float(np.mean(births))
                features[f'dim_{dim}_mean_death'] = float(np.mean(deaths[np.isfinite(deaths)]))
                features[f'dim_{dim}_max_birth'] = float(np.max(births))
                max_finite_death = np.max(deaths[np.isfinite(deaths)]) if np.any(np.isfinite(deaths)) else 0.0
                features[f'dim_{dim}_max_death'] = float(max_finite_death)
            else:
                # No features in this dimension
                for stat in ['num_features', 'max_persistence', 'mean_persistence', 
                           'std_persistence', 'sum_persistence', 'mean_birth', 
                           'mean_death', 'max_birth', 'max_death']:
                    features[f'dim_{dim}_{stat}'] = 0.0
    
    return features

def plot_persistence_diagrams(persistence_diagrams, output_dir):
    """Plot persistence diagrams"""
    print("Plotting persistence diagrams...")
    
    fig, axes = plt.subplots(1, len(persistence_diagrams), 
                            figsize=(5*len(persistence_diagrams), 5))
    
    if len(persistence_diagrams) == 1:
        axes = [axes]
    
    colors = ['red', 'blue', 'green', 'orange']
    
    for dim, (ax, diagram) in enumerate(zip(axes, persistence_diagrams)):
        if len(diagram) > 0:
            births = diagram[:, 0]
            deaths = diagram[:, 1]
            
            # Remove infinite points for plotting
            finite_mask = np.isfinite(deaths)
            births_finite = births[finite_mask]
            deaths_finite = deaths[finite_mask]
            
            # Plot finite points
            if len(births_finite) > 0:
                ax.scatter(births_finite, deaths_finite, 
                          c=colors[dim % len(colors)], alpha=0.7, s=50)
            
            # Plot infinite points
            if np.any(~finite_mask):
                births_inf = births[~finite_mask]
                max_death = np.max(deaths_finite) if len(deaths_finite) > 0 else 1
                ax.scatter(births_inf, [max_death * 1.1] * len(births_inf),
                          c=colors[dim % len(colors)], marker='^', s=100, alpha=0.7)
            
            # Diagonal line
            max_val = max(np.max(births), np.max(deaths_finite) if len(deaths_finite) > 0 else 1)
            ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.5)
            
        ax.set_xlabel('Birth')
        ax.set_ylabel('Death')
        ax.set_title(f'H_{dim} Persistence Diagram')
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'persistence_diagrams.png'), 
                dpi=300, bbox_inches='tight')
    plt.close()

def save_results(persistence_diagrams, features, connectivity_matrix, output_dir):
    """Save results to output directory"""
    print("Saving results...")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Save persistence diagrams as numpy arrays
    for dim, diagram in enumerate(persistence_diagrams):
        np.save(os.path.join(output_dir, f'persistence_diagram_dim{dim}.npy'), diagram)
    
    # Save features as JSON and CSV
    with open(os.path.join(output_dir, 'persistence_features.json'), 'w') as f:
        json.dump(features, f, indent=2)
    
    # Convert features to DataFrame and save as CSV
    features_df = pd.DataFrame([features])
    features_df.to_csv(os.path.join(output_dir, 'persistence_features.csv'), index=False)
    
    # Save connectivity matrix
    np.save(os.path.join(output_dir, 'connectivity_matrix.npy'), connectivity_matrix)
    
    # Plot connectivity matrix
    plt.figure(figsize=(10, 8))
    plt.imshow(connectivity_matrix, cmap='RdBu_r', vmin=-1, vmax=1)
    plt.colorbar(label='Correlation')
    plt.title('Functional Connectivity Matrix')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'connectivity_matrix.png'), 
                dpi=300, bbox_inches='tight')
    plt.close()

def create_product_json(features, output_dir):
    """Create product.json file for Brainlife"""
    product = {
        'brainlife': [
            {
                'type': 'image',
                'name': 'Persistence Diagrams',
                'base64': None,
                'filename': 'persistence_diagrams.png'
            },
            {
                'type': 'image', 
                'name': 'Connectivity Matrix',
                'base64': None,
                'filename': 'connectivity_matrix.png'
            }
        ],
        'persistent_homology_summary': {
            'total_features_dim0': features.get('dim_0_num_features', 0),
            'total_features_dim1': features.get('dim_1_num_features', 0),
            'max_persistence_dim0': features.get('dim_0_max_persistence', 0),
            'max_persistence_dim1': features.get('dim_1_max_persistence', 0)
        }
    }
    
    with open(os.path.join(output_dir, 'product.json'), 'w') as f:
        json.dump(product, f, indent=2)

def main():
    print("Starting fMRI Persistent Homology Analysis...")
    
    # Load configuration
    config = load_config()
    print(f"Configuration: {config}")
    
    # Get input parameters
    fmri_path = config.get('fmri')
    if not fmri_path:
        print("Error: 'fmri' path not specified in config.json")
        sys.exit(1)
    
    atlas_path = config.get('atlas', None)
    n_regions = config.get('n_regions', 100)
    maxdim = config.get('maxdim', 1)
    connectivity_method = config.get('connectivity_method', 'correlation')
    
    # Create output directory
    output_dir = 'output'
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        # Load fMRI data
        fmri_img = load_fmri_data(fmri_path)
        
        # Extract time series
        time_series = extract_time_series_advanced(fmri_img, atlas_path, n_regions)
        
        # Compute functional connectivity
        connectivity_matrix = compute_functional_connectivity(time_series, connectivity_method)
        
        # Convert correlation to distance
        distance_matrix = 1 - np.abs(connectivity_matrix)
        
        # Ensure distance matrix is valid
        np.fill_diagonal(distance_matrix, 0)
        distance_matrix = np.maximum(distance_matrix, 0)
        
        # Compute persistent homology
        persistence_diagrams, diagrams = compute_persistent_homology(distance_matrix, maxdim)
        
        # Extract features
        features = compute_persistence_features(persistence_diagrams)
        
        # Plot results
        plot_persistence_diagrams(persistence_diagrams, output_dir)
        
        # Save results
        save_results(persistence_diagrams, features, connectivity_matrix, output_dir)
        
        # Create product.json for Brainlife
        create_product_json(features, output_dir)
        
        print("Analysis completed successfully!")
        print(f"Results saved to: {output_dir}")
        print(f"Extracted {len(features)} features from persistent homology")
        
        # Print summary
        for key, value in features.items():
            if 'num_features' in key:
                print(f"{key}: {value}")
        
    except Exception as e:
        print(f"Error during analysis: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
