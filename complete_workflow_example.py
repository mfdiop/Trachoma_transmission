#!/usr/bin/env python3
"""
Complete Workflow Example for Genomic Transmission Dynamics

This script demonstrates the complete workflow from raw data to analysis results.
"""

import os
from datetime import datetime

# Import our modules
from genomic_transmission_dynamics import TransmissionDynamicsModel, TransmissionVisualizer
from data_preparation import DataPreprocessor, DataValidator
from config import get_pathogen_config, create_custom_config

def complete_workflow_demo():
    """Demonstrate the complete workflow"""
    
    print("🧬 Genomic Transmission Dynamics - Complete Workflow")
    print("=" * 60)
    
    # Step 1: Data Preparation
    print("\n📁 Step 1: Data Preparation")
    print("-" * 30)
    
    preprocessor = DataPreprocessor()
    
    # Process example data (in real use, replace with your files)
    samples = preprocessor.process_data(
        fasta_file='example_sequences.fasta',
        metadata_file='example_metadata.csv',
        alignment_method='simple'
    )
    
    if not samples:
        print("❌ No samples loaded. Please check your data files.")
        return
    
    # Step 2: Data Validation
    print("\n✅ Step 2: Data Validation")
    print("-" * 30)
    
    validator = DataValidator()
    issues = validator.validate_samples(samples)
    validator.print_validation_report(issues)
    
    if issues['errors']:
        print("❌ Data validation failed. Please fix errors before proceeding.")
        return
    
    # Step 3: Model Configuration
    print("\n⚙️  Step 3: Model Configuration")
    print("-" * 30)
    
    # Example: Configure for COVID-19
    covid_config = get_pathogen_config('covid19')
    print("Using COVID-19 configuration:")
    for key, value in covid_config.items():
        print(f"  {key}: {value}")
    
    # Step 4: Model Initialization
    print("\n🔬 Step 4: Model Initialization")
    print("-" * 30)
    
    model = TransmissionDynamicsModel(
        samples=samples,
        **covid_config  # Unpack configuration
    )
    
    print(f"Model initialized with {len(samples)} samples")
    
    # Step 5: Analysis
    print("\n📊 Step 5: Analysis")
    print("-" * 30)
    
    # Calculate key metrics
    n_samples = len(samples)
    
    # Evolution probabilities
    if n_samples >= 2:
        evolution_prob_7d = model.calculate_evolution_probability(0, 1, 7)
        evolution_prob_30d = model.calculate_evolution_probability(0, 1, 30)
        
        print(f"Evolution probabilities ({samples[0].sample_id} → {samples[1].sample_id}):")
        print(f"  7 days:  {evolution_prob_7d:.8f}")
        print(f"  30 days: {evolution_prob_30d:.8f}")
    
    # Transmission network
    network = model.build_transmission_network(probability_threshold=0.001)
    print(f"Transmission network: {network.number_of_nodes()} nodes, {network.number_of_edges()} edges")
    
    # Step 6: Simulation
    print("\n🎯 Step 6: Simulation")
    print("-" * 30)
    
    simulation_results = model.simulate_transmission_dynamics(
        n_steps=50,
        initial_infected=[0]
    )
    
    final_infected = simulation_results['total_infected'][-1]
    attack_rate = final_infected / len(samples) * 100
    
    print(f"Simulation results:")
    print(f"  Final infected: {final_infected}/{len(samples)}")
    print(f"  Attack rate: {attack_rate:.1f}%")
    
    # Step 7: Visualization
    print("\n📈 Step 7: Visualization")
    print("-" * 30)
    
    visualizer = TransmissionVisualizer(model)
    
    output_files = []
    
    try:
        # Create output directory
        os.makedirs('workflow_output', exist_ok=True)
        
        # Generate visualizations
        fig1 = visualizer.plot_spatial_distribution()
        fig1.savefig('workflow_output/spatial_analysis.png', dpi=300, bbox_inches='tight')
        output_files.append('workflow_output/spatial_analysis.png')
        
        fig2 = visualizer.plot_genomic_relationships()
        fig2.savefig('workflow_output/genomic_analysis.png', dpi=300, bbox_inches='tight')
        output_files.append('workflow_output/genomic_analysis.png')
        
        fig3 = visualizer.plot_transmission_network()
        fig3.savefig('workflow_output/network_analysis.png', dpi=300, bbox_inches='tight')
        output_files.append('workflow_output/network_analysis.png')
        
        fig4 = visualizer.plot_simulation_results(simulation_results)
        fig4.savefig('workflow_output/simulation_analysis.png', dpi=300, bbox_inches='tight')
        output_files.append('workflow_output/simulation_analysis.png')
        
        print(f"Generated {len(output_files)} visualization files:")
        for file in output_files:
            print(f"  ✅ {file}")
            
    except Exception as e:
        print(f"❌ Visualization error: {e}")
    
    # Step 8: Results Summary
    print("\n📋 Step 8: Results Summary")
    print("-" * 30)
    
    # Calculate summary statistics
    spatial_distances = model.spatial_distances
    genomic_distances = model.genomic_distances
    
    # Get upper triangle (avoid diagonal and duplicates)
    import numpy as np
    upper_tri = np.triu_indices_from(spatial_distances, k=1)
    spatial_flat = spatial_distances[upper_tri]
    genomic_flat = genomic_distances[upper_tri]
    
    correlation = np.corrcoef(spatial_flat, genomic_flat)[0, 1] if len(spatial_flat) > 1 else 0
    
    print("Summary Statistics:")
    print(f"  Samples analyzed: {len(samples)}")
    print(f"  Villages: {len(set(s.village for s in samples))}")
    print(f"  Date range: {min(s.collection_date for s in samples).strftime('%Y-%m-%d')} to "
          f"{max(s.collection_date for s in samples).strftime('%Y-%m-%d')}")
    print(f"  Average spatial distance: {np.mean(spatial_flat):.2f} km")
    print(f"  Average genomic distance: {np.mean(genomic_flat):.2f} substitutions")
    print(f"  Spatial-genomic correlation: {correlation:.3f}")
    
    # Step 9: Export Results
    print("\n💾 Step 9: Export Results")
    print("-" * 30)
    
    # Create results summary
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    results_file = f'workflow_output/results_summary_{timestamp}.txt'
    
    with open(results_file, 'w') as f:
        f.write("Genomic Transmission Dynamics Analysis Results\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Configuration: COVID-19 parameters\n\n")
        
        f.write("Data Summary:\n")
        f.write(f"  Samples: {len(samples)}\n")
        f.write(f"  Villages: {len(set(s.village for s in samples))}\n")
        f.write(f"  Average spatial distance: {np.mean(spatial_flat):.2f} km\n")
        f.write(f"  Average genomic distance: {np.mean(genomic_flat):.2f} substitutions\n")
        f.write(f"  Spatial-genomic correlation: {correlation:.3f}\n\n")
        
        f.write("Simulation Results:\n")
        f.write(f"  Final infected: {final_infected}/{len(samples)}\n")
        f.write(f"  Attack rate: {attack_rate:.1f}%\n\n")
        
        f.write("Network Analysis:\n")
        f.write(f"  Nodes: {network.number_of_nodes()}\n")
        f.write(f"  Edges: {network.number_of_edges()}\n")
        
        if network.number_of_edges() > 0:
            import networkx as nx
            f.write(f"  Network density: {nx.density(network):.3f}\n")
    
    print(f"✅ Results exported to: {results_file}")
    
    print("\n🎉 Workflow completed successfully!")
    print("\nNext steps:")
    print("1. Review the generated visualizations")
    print("2. Examine the results summary file")
    print("3. Adjust model parameters if needed")
    print("4. Run with your real genomic data")

def demonstrate_parameter_comparison():
    """Demonstrate analysis with different parameter sets"""
    
    print("\n🔬 Parameter Comparison Analysis")
    print("=" * 40)
    
    # Load data
    preprocessor = DataPreprocessor()
    samples = preprocessor.process_data(
        fasta_file='example_sequences.fasta',
        metadata_file='example_metadata.csv'
    )
    
    if len(samples) < 2:
        print("Need at least 2 samples for comparison")
        return
    
    # Test different pathogen configurations
    pathogens = ['covid19', 'influenza', 'malaria', 'tuberculosis']
    results = {}
    
    for pathogen in pathogens:
        config = get_pathogen_config(pathogen)
        
        model = TransmissionDynamicsModel(samples=samples, **config)
        
        # Run short simulation
        sim_results = model.simulate_transmission_dynamics(n_steps=30, initial_infected=[0])
        final_infected = sim_results['total_infected'][-1]
        
        # Calculate evolution probability
        evolution_prob = model.calculate_evolution_probability(0, 1, 14) if len(samples) >= 2 else 0
        
        results[pathogen] = {
            'final_infected': final_infected,
            'attack_rate': final_infected / len(samples) * 100,
            'evolution_prob': evolution_prob
        }
    
    print("Pathogen comparison results:")
    print("Pathogen\t\tFinal Infected\tAttack Rate\tEvolution Prob (14d)")
    print("-" * 70)
    
    for pathogen, result in results.items():
        print(f"{pathogen.upper():<15}\t{result['final_infected']}\t\t"
              f"{result['attack_rate']:.1f}%\t\t{result['evolution_prob']:.2e}")

if __name__ == "__main__":
    # Run complete workflow
    complete_workflow_demo()
    
    # Run parameter comparison
    demonstrate_parameter_comparison()
    
    print("\n" + "=" * 60)
    print("For production use with your data:")
    print("1. Replace 'example_sequences.fasta' with your FASTA file")
    print("2. Replace 'example_metadata.csv' with your metadata CSV")
    print("3. Choose appropriate pathogen configuration")
    print("4. Adjust parameters based on your specific study")
    print("5. Validate results with known transmission patterns")
    print("=" * 60)