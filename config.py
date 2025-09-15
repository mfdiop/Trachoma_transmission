#!/usr/bin/env python3
"""
Configuration file for Genomic Transmission Dynamics Model

Adjust these parameters based on your specific pathogen and study requirements.
"""

# Model Parameters
class ModelConfig:
    """Configuration class for model parameters"""
    
    # Evolution parameters
    MUTATION_RATE = 1e-6          # Per-site per-day mutation rate
    RECOMBINATION_RATE = 1e-7     # Per-site per-day recombination rate
    
    # Transmission parameters
    BASE_TRANSMISSION_RATE = 0.1   # Base transmission probability per day
    SPATIAL_DECAY_DISTANCE = 10.0  # Characteristic distance for spatial decay (km)
    TEMPORAL_DECAY_TIME = 30.0     # Characteristic time for temporal decay (days)
    
    # Recombination analysis parameters
    RECOMBINATION_SPATIAL_WINDOW = 50.0   # Max distance for potential donors (km)
    RECOMBINATION_TEMPORAL_WINDOW = 60    # Max time for potential donors (days)
    
    # Network analysis parameters
    TRANSMISSION_PROBABILITY_THRESHOLD = 0.01  # Min probability for network edges
    
    # Simulation parameters
    SIMULATION_STEPS = 100         # Number of simulation time steps
    INITIAL_INFECTED_FRACTION = 0.1  # Fraction of samples to start as infected

# Pathogen-specific configurations
class PathogenConfigs:
    """Pre-configured parameter sets for different pathogens"""
    
    # SARS-CoV-2 (example parameters)
    COVID19 = {
        'mutation_rate': 1e-6,
        'recombination_rate': 1e-8,  # Lower recombination for RNA viruses
        'transmission_rate': 0.15,
        'spatial_decay_distance': 5.0,  # Shorter range transmission
        'temporal_decay_time': 14.0     # Shorter infectious period
    }
    
    # Influenza (example parameters)
    INFLUENZA = {
        'mutation_rate': 2e-6,      # Higher mutation rate
        'recombination_rate': 5e-7,  # Higher recombination (segmented genome)
        'transmission_rate': 0.2,
        'spatial_decay_distance': 8.0,
        'temporal_decay_time': 10.0
    }
    
    # Malaria parasite (example parameters)
    MALARIA = {
        'mutation_rate': 5e-7,
        'recombination_rate': 1e-6,  # High recombination
        'transmission_rate': 0.05,   # Vector-borne, lower direct transmission
        'spatial_decay_distance': 2.0,  # Localized transmission
        'temporal_decay_time': 60.0  # Longer persistence
    }
    
    # Tuberculosis (example parameters)
    TUBERCULOSIS = {
        'mutation_rate': 3e-7,      # Slower evolution
        'recombination_rate': 1e-8,  # Very low recombination
        'transmission_rate': 0.03,   # Slower transmission
        'spatial_decay_distance': 15.0,  # Household/close contact
        'temporal_decay_time': 90.0  # Long infectious period
    }
    
    # Bacterial pathogen (example parameters)
    BACTERIAL = {
        'mutation_rate': 1e-6,
        'recombination_rate': 1e-6,  # High recombination potential
        'transmission_rate': 0.08,
        'spatial_decay_distance': 12.0,
        'temporal_decay_time': 21.0
    }

# Visualization parameters
class VisualizationConfig:
    """Configuration for visualization settings"""
    
    # Figure settings
    FIGURE_DPI = 300
    FIGURE_FORMAT = 'png'
    FIGURE_SIZE_LARGE = (12, 10)
    FIGURE_SIZE_MEDIUM = (12, 8)
    FIGURE_SIZE_SMALL = (8, 6)
    
    # Color schemes
    VILLAGE_COLORS = {
        'Village_A': '#FF6B6B',  # Red
        'Village_B': '#4ECDC4',  # Teal
        'Village_C': '#45B7D1',  # Blue
        'Village_D': '#96CEB4',  # Green
        'Village_E': '#FFEAA7',  # Yellow
        'Village_F': '#DDA0DD',  # Plum
    }
    
    # Network visualization
    NODE_SIZE = 300
    EDGE_WIDTH_MULTIPLIER = 5
    NETWORK_ALPHA = 0.8
    LABEL_FONT_SIZE = 8
    
    # Heatmap settings
    HEATMAP_COLORMAP = 'viridis'
    GENOMIC_DISTANCE_COLORMAP = 'plasma'

# Analysis parameters
class AnalysisConfig:
    """Configuration for analysis settings"""
    
    # Statistical thresholds
    SIGNIFICANCE_LEVEL = 0.05
    CORRELATION_THRESHOLD = 0.3
    
    # Clustering parameters
    GENOMIC_DISTANCE_THRESHOLD = 10  # Max distance for clustering
    SPATIAL_DISTANCE_THRESHOLD = 50  # Max distance for spatial clustering (km)
    TEMPORAL_DISTANCE_THRESHOLD = 30  # Max time for temporal clustering (days)
    
    # Simulation parameters
    N_BOOTSTRAP_SAMPLES = 1000
    N_SIMULATION_REPLICATES = 100
    
    # Output settings
    SAVE_INTERMEDIATE_RESULTS = True
    OUTPUT_PRECISION = 6  # Decimal places for probabilities

# File I/O configuration
class IOConfig:
    """Configuration for file input/output"""
    
    # Input file formats
    SUPPORTED_SEQUENCE_FORMATS = ['.fasta', '.fa', '.fas', '.fna']
    SUPPORTED_METADATA_FORMATS = ['.csv', '.tsv', '.xlsx']
    
    # Output file naming
    OUTPUT_PREFIX = 'transmission_dynamics'
    TIMESTAMP_FORMAT = '%Y%m%d_%H%M%S'
    
    # Default file paths
    DEFAULT_OUTPUT_DIR = './output'
    DEFAULT_PLOTS_DIR = './plots'
    DEFAULT_DATA_DIR = './data'
    
    # File encoding
    DEFAULT_ENCODING = 'utf-8'

def get_pathogen_config(pathogen_name: str) -> dict:
    """
    Get configuration for a specific pathogen
    
    Args:
        pathogen_name: Name of pathogen ('covid19', 'influenza', etc.)
        
    Returns:
        Dictionary of parameters for the pathogen
    """
    pathogen_name = pathogen_name.upper()
    
    if hasattr(PathogenConfigs, pathogen_name):
        return getattr(PathogenConfigs, pathogen_name)
    else:
        print(f"Warning: Unknown pathogen '{pathogen_name}'. Using default parameters.")
        return {
            'mutation_rate': ModelConfig.MUTATION_RATE,
            'recombination_rate': ModelConfig.RECOMBINATION_RATE,
            'transmission_rate': ModelConfig.BASE_TRANSMISSION_RATE,
            'spatial_decay_distance': ModelConfig.SPATIAL_DECAY_DISTANCE,
            'temporal_decay_time': ModelConfig.TEMPORAL_DECAY_TIME
        }

def create_custom_config(**kwargs) -> dict:
    """
    Create a custom configuration by overriding default parameters
    
    Args:
        **kwargs: Parameter overrides
        
    Returns:
        Dictionary of parameters
    """
    config = {
        'mutation_rate': ModelConfig.MUTATION_RATE,
        'recombination_rate': ModelConfig.RECOMBINATION_RATE,
        'transmission_rate': ModelConfig.BASE_TRANSMISSION_RATE,
        'spatial_decay_distance': ModelConfig.SPATIAL_DECAY_DISTANCE,
        'temporal_decay_time': ModelConfig.TEMPORAL_DECAY_TIME
    }
    
    # Override with provided parameters
    config.update(kwargs)
    
    return config

# Example usage:
if __name__ == "__main__":
    print("Available pathogen configurations:")
    
    for attr_name in dir(PathogenConfigs):
        if not attr_name.startswith('_'):
            config = getattr(PathogenConfigs, attr_name)
            if isinstance(config, dict):
                print(f"\n{attr_name}:")
                for key, value in config.items():
                    print(f"  {key}: {value}")
    
    print("\nExample custom configuration:")
    custom_config = create_custom_config(
        mutation_rate=2e-6,
        transmission_rate=0.15
    )
    for key, value in custom_config.items():
        print(f"  {key}: {value}")