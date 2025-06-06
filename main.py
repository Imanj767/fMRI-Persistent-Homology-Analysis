#!/usr/bin/env python3

import os
import json
import numpy as np
import nibabel as nib
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Import topological data analysis libraries
try:
    import gudhi
    from ripser import ripser
    from persim import plot_diagrams
    TOPOLOGY_AVAILABLE = True
except ImportError:
    TOPOLOGY_AVAILABLE = False
    print("Warning: Topology libraries not available. Installing...")

def load_config():
    """Load configuration from config.json"""
    try:
        with open('config.json', 'r') as f:
            config = json.load(f)
    except FileNotFoundError:
        # Default configuration
        config = {
            'max_dimension': 2,
            'window_size': 50,
            'overlap': 0.5,
            'distance_metric': 'correlation',
            'roi_extraction': 'parcellation'
        }
    return config

def load_fmri_data(func_path):
    """Load fMRI data from NIfTI file"""
    print(f"Loading fMRI data from: {func_path}")
    img = nib.load(func_path)
    data = img.get_fdata()
    
    # Reshape to (timepoints, voxels)
    original_shape = data.shape
    if len(original_shape) == 4:
        data = data.reshape(-1, original_shape[-1]).T
    else:
        raise ValueError("Expected 4D fMRI data (x, y, z, time)")
    
    print(f"fMRI data shape: {data.shape} (timepoints, voxels)")
    return data, img.affine, original_shape

def apply_brain_mask(data, mask_path, original_shape):
    """Apply brain mask to fMRI data"""
    if mask_path and os.path.exists(mask_path):
        print(f"Applying brain mask: {mask_path}")
        mask_img = nib.load(mask_path)
        mask_data = mask_img.get_fdata()
        mask_flat = mask_data.flatten().astype(bool)
        
        # Apply mask to spatial dimension
        data_masked = data[:, mask_flat]
        print(f"After masking: {data_masked.shape}")
        return data_masked
    else:
        print("No mask provided, using all voxels")
        return data

def extract_roi_timeseries(data, method='parcellation', n_components=100):
    """Extract ROI time series using different methods"""
    print(f"Extracting ROI time series using method: {method}")
    
    if method == 'pca':
        # Use PCA to reduce dimensionality
        scaler = StandardScaler()
        data_scaled = scaler.fit_transform(data.T).T
        
        pca = PCA(n_components=min(n_components, data.shape[1]))
        roi_timeseries = pca.fit_transform(data)
        print(f"PCA extracted {roi_timeseries.shape[1]} components")
        
    elif method == 'parcellation':
        # Simple spatial parcellation - divide brain into regions
        n_regions = min(100, data.shape[1] // 100)
        if n_regions == 0:
            n_regions = min(20, data.shape[1])
        
        regions_per_voxel = data.shape[1] // n_regions
        roi_timeseries = []
        
        for i in range(n_regions):
            start_idx = i * regions_per_voxel
            end_idx = min((i + 1) * regions_per_voxel, data.shape[1])
            region_mean = np.mean(data[:, start_idx:end_idx], axis=1)
            roi_timeseries.append(region_mean)
        
        roi_timeseries = np.array(roi_timeseries).T
        print(f"Parcellation extracted {roi_timeseries.shape[1]} regions")
        
    else:  # voxel
        # Use subset of voxels directly
        n_voxels = min(200, data.shape[1])
        step = max(1, data.shape[1] // n_voxels)
        roi_timeseries = data[:, ::step]
        print(f"Voxel method selected {roi_timeseries.shape[1]} voxels")
    
    return roi_timeseries

def sliding_window_analysis(data, window_size, overlap):
    """Create sliding windows for temporal analysis"""
    n_timepoints = data.shape[0]
    step_size = int(window_size * (1 - overlap))
    
    windows = []
    window_starts = []
    
    for start in range(0, n_timepoints - window_size + 1, step_size):
        end = start + window_size
        windows.append(data[start:end, :])
        window_starts.append(start)
    
    print(f"Created {len(windows)} sliding windows")
    return windows, window_starts

def compute_distance_matrix(data, metric='correlation'):
    """Compute distance matrix for persistence computation"""
    if metric == 'correlation':
        # Correlation distance: 1 - |correlation|
        corr_matrix = np.corrcoef(data.T)
        # Handle NaN values
        corr_matrix = np.nan_to_num(corr_matrix, nan=0.0)
        distance_matrix = 1 - np.abs(corr_matrix)
        
    elif metric == 'euclidean':
        # Euclidean distance between time series
        distance_matrix = squareform(pdist(data.T, metric='euclidean'))
        # Normalize to [0, 1]
        distance_matrix = distance_matrix / np.max(distance_matrix)
        
    elif metric == 'cosine':
        # Cosine distance
        distance_matrix = squareform(pdist(data.T, metric='cosine'))
        distance_matrix = np.nan_to_num(distance_matrix, nan=1.0)
    
    return distance_matrix

def compute_persistent_homology(distance_matrix, max_dimension=2):
    """Compute persistent homology using Ripser"""
    print("Computing persistent homology...")
    
    if not TOPOLOGY_AVAILABLE:
        print("Error: Topology libraries not available")
        return None
    
    try:
        # Compute persistence diagrams
        result = ripser(distance_matrix, metric='precomputed', maxdim=max_dimension)
        diagrams = result['dgms']
        
        print(f"Computed persistence diagrams for dimensions 0-{max_dimension}")
        return diagrams, result
        
    except Exception as e:
        print(f"Error computing persistent homology: {e}")
        return None

def extract_topological_features(diagrams):
    """Extract numerical features from persistence diagrams"""
    features = {}
    
    for dim, diagram in enumerate(diagrams):
        if len(diagram) == 0:
            continue
            
        # Birth and death times
        births = diagram[:, 0]
        deaths = diagram[:, 1]
        
        # Persistence (lifetime) of each topological feature
        persistence = deaths - births
        
        # Statistical features
        features[f'dim_{dim}_num_features'] = len(diagram)
        features[f'dim_{dim}_total_persistence'] = np.sum(persistence)
        features[f'dim_{dim}_max_persistence'] = np.max(persistence) if len(persistence) > 0 else 0
        features[f'dim_{dim}_mean_persistence'] = np.mean(persistence) if len(persistence) > 0 else 0
        features[f'dim_{dim}_std_persistence'] = np.std(persistence) if len(persistence) > 0 else 0
        
        # Birth and death statistics
        features[f'dim_{dim}_mean_birth'] = np.mean(births) if len(births) > 0 else 0
        features[f'dim_{dim}_mean_death'] = np.mean(deaths) if len(deaths) > 0 else 0
        
    return features

def save_results(diagrams, features, window_features_list, config):
    """Save results to output files"""
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)
    
    # Save persistence diagrams
    if diagrams is not None:
        np.save(os.path.join(output_dir, "persistence_diagrams.npy"), diagrams, allow_pickle=True)
        
        # Plot persistence diagrams
        try:
            fig, axes = plt.subplots(1, len(diagrams), figsize=(4*len(diagrams), 4))
            if len(diagrams) == 1:
                axes = [axes]
            
            for dim, (ax, diagram) in enumerate(zip(axes, diagrams)):
                if len(diagram) > 0:
                    ax.scatter(diagram[:, 0], diagram[:, 1], alpha=0.6)
                    ax.plot([0, np.max(diagram)], [0, np.max(diagram)], 'k--', alpha=0.5)
                    ax.set_xlabel('Birth')
                    ax.set_ylabel('Death')
                    ax.set_title(f'Dimension {dim}')
                else:
                    ax.text(0.5, 0.5, 'No features', ha='center', va='center', transform=ax.transAxes)
                    ax.set_title(f'Dimension {dim}')
            
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "persistence_diagrams.png"), dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Error plotting diagrams: {e}")
    
    # Save topological features
    if features:
        features_df = pd.DataFrame([features])
        features_df.to_csv(os.path.join(output_dir, "topological_features.csv"), index=False)
    
    # Save window-wise features
    if window_features_list:
        windows_df = pd.DataFrame(window_features_list)
        windows_df.to_csv(os.path.join(output_dir, "window_features.csv"), index=False)
        
        # Plot feature evolution over windows
        try:
            fig, axes = plt.subplots(2, 2, figsize=(12, 8))
            axes = axes.flatten()
            
            feature_names = ['dim_0_num_features', 'dim_0_total_persistence', 
                           'dim_1_num_features', 'dim_1_total_persistence']
            
            for i, feature in enumerate(feature_names):
                if feature in windows_df.columns:
                    axes[i].plot(windows_df[feature])
                    axes[i].set_title(feature.replace('_', ' ').title())
                    axes[i].set_xlabel('Window')
                    axes[i].set_ylabel('Value')
            
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "feature_evolution.png"), dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Error plotting feature evolution: {e}")
    
    # Save configuration and summary
    summary = {
        'config': config,
        'n_windows': len(window_features_list) if window_features_list else 0,
        'global_features': features if features else {},
        'processing_completed': True
    }
    
    with open(os.path.join(output_dir, "analysis_summary.json"), 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"Results saved to {output_dir}/")

def main():
    """Main processing function"""
    print("Starting fMRI Persistent Homology Analysis")
    
    # Load configuration
    config = load_config()
    print(f"Configuration: {config}")
    
    # Check for required libraries
    if not TOPOLOGY_AVAILABLE:
        print("Installing required topology libraries...")
        os.system("pip install gudhi ripser persim")
        try:
            import gudhi
            from ripser import ripser
            from persim import plot_diagrams
            TOPOLOGY_AVAILABLE = True
        except ImportError:
            print("Failed to install topology libraries. Exiting.")
            return
    
    # Load fMRI data
    func_path = "func/func.nii.gz"
    if not os.path.exists(func_path):
        func_path = "func.nii.gz"
    
    if not os.path.exists(func_path):
        print("Error: fMRI data file not found")
        return
    
    data, affine, original_shape = load_fmri_data(func_path)
    
    # Apply brain mask if available
    mask_path = "mask/mask.nii.gz" if os.path.exists("mask/mask.nii.gz") else None
    data = apply_brain_mask(data, mask_path, original_shape)
    
    # Extract ROI time series
    roi_data = extract_roi_timeseries(data, method=config['roi_extraction'])
    
    # Remove NaN and infinite values
    roi_data = np.nan_to_num(roi_data, nan=0.0, posinf=0.0, neginf=0.0)
    
    # Global analysis on entire time series
    print("Computing global persistent homology...")
    global_distance_matrix = compute_distance_matrix(roi_data, config['distance_metric'])
    global_result = compute_persistent_homology(global_distance_matrix, config['max_dimension'])
    
    global_features = {}
    if global_result:
        global_diagrams, _ = global_result
        global_features = extract_topological_features(global_diagrams)
    
    # Sliding window analysis
    print("Starting sliding window analysis...")
    windows, window_starts = sliding_window_analysis(
        roi_data, config['window_size'], config['overlap']
    )
    
    window_features_list = []
    
    for i, (window_data, start_time) in enumerate(zip(windows, window_starts)):
        print(f"Processing window {i+1}/{len(windows)}")
        
        # Compute distance matrix for this window
        window_distance_matrix = compute_distance_matrix(window_data, config['distance_metric'])
        
        # Compute persistent homology
        window_result = compute_persistent_homology(window_distance_matrix, config['max_dimension'])
        
        if window_result:
            window_diagrams, _ = window_result
            window_features = extract_topological_features(window_diagrams)
            window_features['window_id'] = i
            window_features['start_time'] = start_time
            window_features_list.append(window_features)
    
    # Save all results
    diagrams_to_save = global_diagrams if global_result else None
    save_results(diagrams_to_save, global_features, window_features_list, config)
    
    print("Analysis completed successfully!")

if __name__ == "__main__":
    main()
