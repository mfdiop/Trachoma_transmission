# Genomic Transmission Dynamics Model for Tanzanian Villages

This repository contains a comprehensive Python script for modeling the transmission dynamics of pathogens across villages in Tanzania, incorporating genomic data, spatial-temporal constraints, and recombination events.

## Features

### Core Functionality
- **Person-to-person transmission modeling** with spatial and temporal decay functions
- **Recombination modeling** that considers available donor sequences in the population
- **Evolution probability calculation** for genome A → genome B transitions over time
- **Spatial-temporal analysis** using GPS coordinates and collection dates
- **Network-based transmission modeling** with NetworkX integration

### Key Models Implemented

1. **Transmission Probability Model**
   - Exponential decay with spatial distance (characteristic distance: 10 km)
   - Exponential decay with temporal distance (characteristic time: 30 days)
   - Configurable base transmission rate

2. **Evolution Probability Model**
   - Mutation component using Poisson distribution
   - Recombination component based on donor sequence availability
   - Combined probability calculation with spatial-temporal constraints

3. **Recombination Model**
   - Identifies potential donor sequences within 50 km and 60 days
   - Calculates explanation scores for observed differences
   - Incorporates recombination rate and donor availability

## Installation

```bash
pip install -r requirements.txt
```

## Usage

### Basic Usage

```python
from genomic_transmission_dynamics import TransmissionDynamicsModel, GenomicSample, TransmissionVisualizer
from datetime import datetime

# Create your genomic samples
samples = [
    GenomicSample(
        sample_id="sample_001",
        genome_sequence="ATCG...",  # Your actual genome sequence
        village="Village_A",
        coordinates=(-6.1659, 35.7516),  # (latitude, longitude)
        collection_date=datetime(2024, 1, 15),
        individual_id="IND_001"
    ),
    # ... more samples
]

# Initialize model
model = TransmissionDynamicsModel(
    samples=samples,
    mutation_rate=1e-6,      # per site per day
    recombination_rate=1e-7,  # per site per day
    transmission_rate=0.1     # base transmission probability
)

# Calculate evolution probability
prob = model.calculate_evolution_probability(
    genome_a_idx=0,
    genome_b_idx=1,
    time_period=30  # days
)

print(f"Evolution probability: {prob}")
```

### Running the Demo

```bash
python genomic_transmission_dynamics.py
```

This will:
1. Generate sample data for 3 Tanzanian villages (30 samples total)
2. Calculate transmission and evolution probabilities
3. Build transmission networks
4. Run epidemic simulations
5. Generate visualization plots

## Data Structure

### GenomicSample
Each genomic sample should include:
- `sample_id`: Unique identifier
- `genome_sequence`: DNA/RNA sequence string
- `village`: Village name
- `coordinates`: (latitude, longitude) tuple
- `collection_date`: datetime object
- `individual_id`: Individual identifier

### Model Parameters
- `mutation_rate`: Per-site per-day mutation rate (default: 1e-6)
- `recombination_rate`: Per-site per-day recombination rate (default: 1e-7)
- `transmission_rate`: Base transmission probability per day (default: 0.1)

## Key Methods

### TransmissionDynamicsModel

- `calculate_transmission_probability(source_idx, target_idx)`: Calculate transmission probability between samples
- `calculate_evolution_probability(genome_a_idx, genome_b_idx, time_period)`: Calculate probability of genome evolution
- `build_transmission_network(probability_threshold)`: Build NetworkX transmission network
- `simulate_transmission_dynamics(n_steps, initial_infected)`: Run epidemic simulation

### TransmissionVisualizer

- `plot_spatial_distribution()`: Plot sample locations and spatial distances
- `plot_genomic_relationships()`: Plot genomic distance matrices and correlations
- `plot_transmission_network()`: Visualize transmission network
- `plot_simulation_results()`: Plot epidemic simulation results

## Output Files

The demo generates several visualization files:
- `spatial_distribution.png`: Sample locations and spatial distance matrix
- `genomic_relationships.png`: Genomic distances and spatial correlations
- `transmission_network.png`: Network visualization of transmission probabilities
- `simulation_results.png`: Epidemic simulation results over time

## Model Assumptions

1. **Spatial Decay**: Transmission probability decreases exponentially with distance
2. **Temporal Decay**: Transmission probability decreases exponentially with time
3. **Recombination**: Requires donor sequences within spatial-temporal windows
4. **Mutation**: Follows Poisson distribution based on mutation rate and time
5. **Related Genomes**: Show correlation between genetic and spatial-temporal proximity

## Customization

### Adjusting Model Parameters
```python
model = TransmissionDynamicsModel(
    samples=your_samples,
    mutation_rate=2e-6,      # Higher mutation rate
    recombination_rate=5e-8,  # Lower recombination rate
    transmission_rate=0.05    # Lower transmission rate
)
```

### Custom Distance Functions
The model uses haversine distance for spatial calculations and Hamming distance for genomic differences. These can be modified by overriding the respective methods.

### Network Analysis
```python
network = model.build_transmission_network(probability_threshold=0.01)
print(f"Network density: {nx.density(network)}")
print(f"Connected components: {nx.number_weakly_connected_components(network)}")
```

## Scientific Background

This model implements concepts from:
- Phylogeographic analysis
- Epidemic modeling
- Population genetics
- Spatial epidemiology
- Recombination detection in pathogen genomes

The approach is particularly suitable for:
- Viral transmission studies
- Bacterial outbreak investigations
- Parasite population dynamics
- Any pathogen with measurable genomic variation

## Contributing

Feel free to extend this model by:
- Adding more sophisticated recombination models
- Implementing phylogenetic tree-based approaches
- Adding demographic and environmental factors
- Incorporating selection pressure models

## License

This code is provided for research and educational purposes.