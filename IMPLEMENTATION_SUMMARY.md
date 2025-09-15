# Genomic Transmission Dynamics Implementation Summary

## 🎯 Project Overview

I have successfully created a comprehensive script for modeling transmission dynamics of genomic data from 3 villages in Tanzania. The implementation includes person-to-person transmission modeling, recombination analysis, and spatial-temporal constraints.

## 📁 Files Created

### Core Implementation
1. **`genomic_transmission_dynamics.py`** - Main model implementation
   - `TransmissionDynamicsModel` class with all core functionality
   - `TransmissionVisualizer` class for generating plots
   - Sample data generation and main execution

2. **`data_preparation.py`** - Data preprocessing utilities
   - `DataPreprocessor` class for handling FASTA and CSV files
   - `DataValidator` class for data quality checks
   - Sequence alignment and metadata processing

3. **`config.py`** - Configuration management
   - Pre-configured parameter sets for different pathogens
   - Customizable model parameters
   - Visualization and analysis settings

### Usage Examples
4. **`example_usage.py`** - Custom data analysis examples
5. **`complete_workflow_example.py`** - End-to-end workflow demonstration

### Supporting Files
6. **`requirements.txt`** - Python package dependencies
7. **`README.md`** - Comprehensive documentation
8. **Example data files** - Sample FASTA and CSV files for testing

## 🔬 Key Features Implemented

### Core Model Components

#### 1. Person-to-Person Transmission Model
- **Spatial decay**: Exponential decay with configurable characteristic distance
- **Temporal decay**: Exponential decay with configurable characteristic time
- **Combined probability**: Multiplicative model incorporating both factors
- **Configurable parameters**: Base transmission rate, decay distances/times

#### 2. Recombination Model
- **Donor identification**: Finds potential donor sequences within spatial-temporal windows
- **Explanation scoring**: Calculates how well donors explain observed differences
- **Recombination probability**: Based on donor availability and sequence similarity
- **Spatial-temporal constraints**: Only considers nearby sequences as potential donors

#### 3. Evolution Probability Calculation
- **Mutation component**: Poisson model based on mutation rate and time
- **Recombination component**: Incorporates donor-based recombination
- **Combined model**: Multiplies mutation and recombination probabilities
- **Time-dependent**: Calculates probability over specified time periods

#### 4. Spatial-Temporal Analysis
- **Haversine distance**: Accurate geographic distance calculations
- **Date handling**: Flexible date parsing and temporal distance calculation
- **Coordinate validation**: GPS coordinate validation and error checking
- **Multi-village support**: Handles samples from multiple locations

### Advanced Features

#### 5. Network Analysis
- **Transmission networks**: NetworkX-based directed graphs
- **Probability thresholds**: Configurable minimum probabilities for edges
- **Network metrics**: Density, clustering, connectivity analysis
- **Visualization**: Geographic layout based on actual coordinates

#### 6. Simulation Framework
- **Dynamic simulation**: Time-step based transmission simulation
- **Multiple starting points**: Configurable initial infected samples
- **Attack rate calculation**: Proportion of population infected
- **Infection timing**: Tracks when each sample becomes infected

#### 7. Comprehensive Visualization
- **Spatial distribution**: Sample locations colored by village
- **Distance matrices**: Heatmaps of spatial and genomic distances
- **Correlation plots**: Spatial vs genomic distance relationships
- **Network graphs**: Transmission probability networks
- **Simulation results**: Time series and infection patterns

## 📊 Model Parameters

### Default Configuration
```python
mutation_rate = 1e-6          # Per-site per-day
recombination_rate = 1e-7     # Per-site per-day  
transmission_rate = 0.1       # Base probability per day
spatial_decay_distance = 10.0 # km
temporal_decay_time = 30.0    # days
```

### Pathogen-Specific Configurations
- **COVID-19**: Higher transmission rate, shorter temporal window
- **Influenza**: Higher mutation and recombination rates
- **Malaria**: Lower transmission (vector-borne), high recombination
- **Tuberculosis**: Slow evolution, long infectious period
- **Bacterial**: High recombination potential

## 🧬 Scientific Approach

### Genome Evolution Model
The model calculates P(A→B|t) using:

1. **Mutation Probability**:
   ```
   P_mut = Poisson(μ * L * t, k)
   ```
   Where μ = mutation rate, L = genome length, t = time, k = observed differences

2. **Recombination Probability**:
   ```
   P_rec = exp(-ρ * L * t) + (1 - exp(-ρ * L * t)) * S_best
   ```
   Where ρ = recombination rate, S_best = best donor explanation score

3. **Combined Evolution**:
   ```
   P_evolution = P_mut * P_rec
   ```

### Transmission Model
Transmission probability between samples i and j:

```
P_trans(i→j) = β * exp(-d_spatial/λ_s) * exp(-d_temporal/λ_t)
```

Where:
- β = base transmission rate
- d_spatial = geographic distance
- d_temporal = time difference
- λ_s, λ_t = characteristic decay distances/times

### Related Genomes Constraint
The model ensures related genomes are close in space and time by:
- Using spatial-temporal proximity in transmission probabilities
- Requiring donor sequences within defined windows for recombination
- Incorporating distance-based decay functions
- Validating correlation between genomic and geographic distances

## 📈 Output and Analysis

### Generated Visualizations
1. **Spatial Distribution**: Sample locations and distance matrices
2. **Genomic Relationships**: Sequence distances and correlations
3. **Transmission Network**: Probability-based network graphs
4. **Simulation Results**: Epidemic dynamics over time

### Analysis Results
- Evolution probabilities for genome A→B transitions
- Transmission probabilities between all sample pairs
- Network connectivity and clustering metrics
- Simulation attack rates and infection timing
- Spatial-genomic correlation coefficients

### Export Formats
- High-resolution PNG plots (300 DPI)
- Structured results summary files
- CSV-compatible data formats
- NetworkX graph objects for further analysis

## 🚀 Usage Instructions

### Quick Start
```bash
# Install dependencies
pip install -r requirements.txt

# Run demo with sample data
python3 genomic_transmission_dynamics.py

# Run complete workflow
python3 complete_workflow_example.py
```

### Using Your Own Data
```python
# Prepare your data
from data_preparation import DataPreprocessor
preprocessor = DataPreprocessor()
samples = preprocessor.process_data('your_sequences.fasta', 'your_metadata.csv')

# Configure model
from genomic_transmission_dynamics import TransmissionDynamicsModel
from config import get_pathogen_config

config = get_pathogen_config('covid19')  # or your pathogen
model = TransmissionDynamicsModel(samples, **config)

# Analyze evolution probability
prob = model.calculate_evolution_probability(genome_a_idx=0, genome_b_idx=1, time_period=30)
```

## 📋 Data Requirements

### Input Files
1. **FASTA file** with aligned genomic sequences
2. **CSV file** with metadata containing:
   - sample_id
   - village
   - latitude, longitude (GPS coordinates)
   - collection_date (YYYY-MM-DD format)
   - individual_id (optional)

### Data Validation
- Sequence length consistency
- Valid GPS coordinates (-90≤lat≤90, -180≤lon≤180)
- Parseable dates
- No duplicate sample IDs
- Minimum 2 samples for analysis

## 🔧 Customization Options

### Model Parameters
- Adjust mutation and recombination rates for your pathogen
- Modify spatial and temporal decay characteristics
- Set transmission probability thresholds
- Configure recombination donor windows

### Visualization
- Customize color schemes for villages
- Adjust plot sizes and DPI
- Modify network layout algorithms
- Configure heatmap color maps

### Analysis
- Set statistical significance thresholds
- Adjust clustering parameters
- Configure simulation parameters
- Modify bootstrap sample sizes

## ✅ Validation and Testing

The implementation includes:
- Comprehensive data validation
- Parameter sensitivity analysis
- Multiple pathogen configurations
- Example data generation
- Error handling and warnings
- Results verification checks

## 🎯 Key Achievements

1. ✅ **Person-to-person transmission modeling** with spatial-temporal constraints
2. ✅ **Recombination analysis** using donor sequence availability
3. ✅ **Evolution probability calculation** for A→B genome transitions
4. ✅ **Spatial-temporal correlation** ensuring related genomes are proximate
5. ✅ **Multi-village analysis** supporting 3+ Tanzanian villages
6. ✅ **Comprehensive visualization** with publication-ready plots
7. ✅ **Flexible configuration** for different pathogens
8. ✅ **Complete workflow** from raw data to results
9. ✅ **Scientific validation** with proper statistical models
10. ✅ **Production-ready code** with error handling and documentation

The implementation successfully addresses all requirements and provides a robust, scientifically-sound framework for analyzing genomic transmission dynamics in Tanzanian village populations.