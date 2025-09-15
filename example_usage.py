#!/usr/bin/env python3
"""
Example Usage of Genomic Transmission Dynamics Model

This script demonstrates how to use the model with your own genomic data
from Tanzanian villages.
"""

from genomic_transmission_dynamics import (
    TransmissionDynamicsModel, 
    GenomicSample, 
    TransmissionVisualizer
)
from datetime import datetime
import numpy as np

def load_your_genomic_data():
    """
    Replace this function with your actual data loading logic.
    
    Your data should include:
    - Genomic sequences (DNA/RNA)
    - GPS coordinates (latitude, longitude)
    - Collection dates
    - Village/location information
    """
    
    # Example of how your data structure should look
    samples = [
        GenomicSample(
            sample_id="TZ_001",
            genome_sequence="ATCGATCGATCG...",  # Your actual sequence
            village="Dodoma_Village",
            coordinates=(-6.1659, 35.7516),  # Actual GPS coordinates
            collection_date=datetime(2024, 3, 15),
            individual_id="Patient_001"
        ),
        GenomicSample(
            sample_id="TZ_002", 
            genome_sequence="ATCGATCGATCG...",  # Different sequence
            village="Dodoma_Village",
            coordinates=(-6.1670, 35.7520),  # Nearby location
            collection_date=datetime(2024, 3, 20),
            individual_id="Patient_002"
        ),
        # Add more samples...
    ]
    
    return samples

def analyze_custom_data():
    """Analyze your custom genomic data"""
    
    print("Loading your genomic data...")
    samples = load_your_genomic_data()
    
    if len(samples) < 2:
        print("Error: Need at least 2 samples for analysis")
        return
    
    print(f"Loaded {len(samples)} samples")
    
    # Initialize model with custom parameters
    model = TransmissionDynamicsModel(
        samples=samples,
        mutation_rate=1e-6,       # Adjust based on your pathogen
        recombination_rate=1e-7,   # Adjust based on your pathogen  
        transmission_rate=0.1      # Adjust based on transmission characteristics
    )
    
    # Example analyses
    print("\n--- Transmission Analysis ---")
    
    # Calculate transmission probabilities between all pairs
    n_samples = len(samples)
    high_prob_pairs = []
    
    for i in range(n_samples):
        for j in range(n_samples):
            if i != j:
                prob = model.calculate_transmission_probability(i, j)
                if prob > 0.01:  # High probability threshold
                    high_prob_pairs.append((i, j, prob))
    
    print(f"Found {len(high_prob_pairs)} high-probability transmission pairs")
    
    # Show top 5 transmission pairs
    high_prob_pairs.sort(key=lambda x: x[2], reverse=True)
    print("\nTop transmission probability pairs:")
    for i, (source, target, prob) in enumerate(high_prob_pairs[:5]):
        source_sample = samples[source]
        target_sample = samples[target]
        print(f"{i+1}. {source_sample.sample_id} → {target_sample.sample_id}: {prob:.4f}")
    
    print("\n--- Evolution Analysis ---")
    
    # Calculate evolution probabilities for different time periods
    time_periods = [7, 14, 30, 60]  # days
    
    if len(samples) >= 2:
        sample_a_idx = 0
        sample_b_idx = 1
        
        print(f"Evolution probabilities for {samples[sample_a_idx].sample_id} → {samples[sample_b_idx].sample_id}:")
        
        for days in time_periods:
            prob = model.calculate_evolution_probability(sample_a_idx, sample_b_idx, days)
            print(f"  {days} days: {prob:.8f}")
    
    print("\n--- Network Analysis ---")
    
    # Build transmission network
    network = model.build_transmission_network(probability_threshold=0.001)
    
    print(f"Network statistics:")
    print(f"  Nodes: {network.number_of_nodes()}")
    print(f"  Edges: {network.number_of_edges()}")
    
    if network.number_of_edges() > 0:
        import networkx as nx
        print(f"  Average clustering: {nx.average_clustering(network):.3f}")
        print(f"  Network density: {nx.density(network):.3f}")
    
    print("\n--- Simulation ---")
    
    # Run transmission simulation
    simulation_results = model.simulate_transmission_dynamics(
        n_steps=100,
        initial_infected=[0]  # Start with first sample
    )
    
    final_infected = simulation_results['total_infected'][-1]
    print(f"Simulation results:")
    print(f"  Final infected: {final_infected}/{len(samples)}")
    print(f"  Attack rate: {final_infected/len(samples)*100:.1f}%")
    
    # Generate visualizations
    print("\n--- Generating Visualizations ---")
    
    visualizer = TransmissionVisualizer(model)
    
    try:
        # Spatial distribution
        fig1 = visualizer.plot_spatial_distribution()
        fig1.savefig('custom_spatial_distribution.png', dpi=300, bbox_inches='tight')
        print("  Saved: custom_spatial_distribution.png")
        
        # Genomic relationships  
        fig2 = visualizer.plot_genomic_relationships()
        fig2.savefig('custom_genomic_relationships.png', dpi=300, bbox_inches='tight')
        print("  Saved: custom_genomic_relationships.png")
        
        # Transmission network
        fig3 = visualizer.plot_transmission_network()
        fig3.savefig('custom_transmission_network.png', dpi=300, bbox_inches='tight')
        print("  Saved: custom_transmission_network.png")
        
        # Simulation results
        fig4 = visualizer.plot_simulation_results(simulation_results)
        fig4.savefig('custom_simulation_results.png', dpi=300, bbox_inches='tight')
        print("  Saved: custom_simulation_results.png")
        
    except Exception as e:
        print(f"  Visualization error: {e}")
    
    return model, simulation_results

def parameter_sensitivity_analysis():
    """Demonstrate parameter sensitivity analysis"""
    
    print("\n=== Parameter Sensitivity Analysis ===")
    
    # Create test data
    samples = load_your_genomic_data()
    if len(samples) < 2:
        print("Need at least 2 samples for sensitivity analysis")
        return
    
    # Test different parameter combinations
    mutation_rates = [1e-7, 1e-6, 1e-5]
    transmission_rates = [0.01, 0.1, 0.2]
    
    results = []
    
    for mut_rate in mutation_rates:
        for trans_rate in transmission_rates:
            model = TransmissionDynamicsModel(
                samples=samples,
                mutation_rate=mut_rate,
                recombination_rate=1e-7,
                transmission_rate=trans_rate
            )
            
            # Run short simulation
            sim_results = model.simulate_transmission_dynamics(n_steps=50)
            final_infected = sim_results['total_infected'][-1]
            
            results.append({
                'mutation_rate': mut_rate,
                'transmission_rate': trans_rate,
                'final_infected': final_infected,
                'attack_rate': final_infected / len(samples)
            })
    
    print("Parameter sensitivity results:")
    print("Mut_Rate\tTrans_Rate\tFinal_Infected\tAttack_Rate")
    print("-" * 60)
    
    for result in results:
        print(f"{result['mutation_rate']:.1e}\t{result['transmission_rate']:.2f}\t\t"
              f"{result['final_infected']}\t\t{result['attack_rate']:.3f}")

def village_comparison_analysis():
    """Compare transmission dynamics between villages"""
    
    print("\n=== Village Comparison Analysis ===")
    
    samples = load_your_genomic_data()
    if len(samples) < 3:
        print("Need at least 3 samples for village comparison")
        return
    
    model = TransmissionDynamicsModel(samples)
    
    # Group samples by village
    villages = {}
    for i, sample in enumerate(samples):
        if sample.village not in villages:
            villages[sample.village] = []
        villages[sample.village].append(i)
    
    print(f"Villages found: {list(villages.keys())}")
    
    # Calculate within-village vs between-village transmission probabilities
    within_village_probs = []
    between_village_probs = []
    
    for village1, indices1 in villages.items():
        for village2, indices2 in villages.items():
            for i in indices1:
                for j in indices2:
                    if i != j:
                        prob = model.calculate_transmission_probability(i, j)
                        
                        if village1 == village2:
                            within_village_probs.append(prob)
                        else:
                            between_village_probs.append(prob)
    
    if within_village_probs and between_village_probs:
        print(f"\nTransmission probabilities:")
        print(f"  Within villages - Mean: {np.mean(within_village_probs):.6f}, "
              f"Max: {np.max(within_village_probs):.6f}")
        print(f"  Between villages - Mean: {np.mean(between_village_probs):.6f}, "
              f"Max: {np.max(between_village_probs):.6f}")
        
        ratio = np.mean(within_village_probs) / np.mean(between_village_probs) if np.mean(between_village_probs) > 0 else float('inf')
        print(f"  Within/Between ratio: {ratio:.2f}")

if __name__ == "__main__":
    print("Genomic Transmission Dynamics - Custom Analysis Example")
    print("=" * 60)
    
    # Note: This example uses placeholder data
    # Replace load_your_genomic_data() with your actual data loading logic
    
    try:
        # Main analysis
        model, simulation_results = analyze_custom_data()
        
        # Additional analyses
        parameter_sensitivity_analysis()
        village_comparison_analysis()
        
        print("\nAnalysis complete!")
        print("\nTo use with your real data:")
        print("1. Modify load_your_genomic_data() function")
        print("2. Ensure your sequences are aligned")
        print("3. Provide accurate GPS coordinates")
        print("4. Include precise collection dates")
        print("5. Adjust model parameters for your pathogen")
        
    except Exception as e:
        print(f"Analysis failed: {e}")
        print("This is expected with placeholder data. Please provide real genomic data.")