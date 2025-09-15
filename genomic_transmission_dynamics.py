#!/usr/bin/env python3
"""
Genomic Transmission Dynamics Model for Tanzanian Villages

This script models the transmission dynamics of a pathogen across 3 villages in Tanzania,
incorporating person-to-person transmission, recombination events, and spatial-temporal
constraints. The model calculates the likelihood of genome A evolving to genome B over time.

Author: AI Assistant
Date: September 15, 2025
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from scipy.optimize import minimize
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional
import networkx as nx
from datetime import datetime, timedelta
import warnings
warnings.filterwarnings('ignore')

@dataclass
class GenomicSample:
    """Represents a genomic sample with metadata"""
    sample_id: str
    genome_sequence: str
    village: str
    coordinates: Tuple[float, float]  # (latitude, longitude)
    collection_date: datetime
    individual_id: str

class TransmissionDynamicsModel:
    """
    Models genomic transmission dynamics with recombination and spatial-temporal constraints
    """
    
    def __init__(self, samples: List[GenomicSample], mutation_rate: float = 1e-6,
                 recombination_rate: float = 1e-7, transmission_rate: float = 0.1,
                 spatial_decay_distance: float = 10.0, temporal_decay_time: float = 30.0,
                 **kwargs):
        """
        Initialize the transmission dynamics model
        
        Args:
            samples: List of genomic samples
            mutation_rate: Per-site per-day mutation rate
            recombination_rate: Per-site per-day recombination rate
            transmission_rate: Base transmission probability per day
            spatial_decay_distance: Characteristic distance for spatial decay (km)
            temporal_decay_time: Characteristic time for temporal decay (days)
            **kwargs: Additional parameters (ignored)
        """
        self.samples = samples
        self.mutation_rate = mutation_rate
        self.recombination_rate = recombination_rate
        self.transmission_rate = transmission_rate
        self.spatial_decay_distance = spatial_decay_distance
        self.temporal_decay_time = temporal_decay_time
        
        # Create distance matrices
        self._create_distance_matrices()
        
        # Calculate genomic distances
        self._calculate_genomic_distances()
    
    def _create_distance_matrices(self):
        """Create spatial and temporal distance matrices"""
        n_samples = len(self.samples)
        
        # Spatial distance matrix (in km)
        coords = np.array([s.coordinates for s in self.samples])
        self.spatial_distances = self._haversine_distance_matrix(coords)
        
        # Temporal distance matrix (in days)
        dates = np.array([s.collection_date for s in self.samples])
        self.temporal_distances = np.abs(dates[:, None] - dates[None, :])
        self.temporal_distances = np.array([[td.days for td in row] for row in self.temporal_distances])
    
    def _haversine_distance_matrix(self, coords: np.ndarray) -> np.ndarray:
        """Calculate haversine distance matrix between coordinates"""
        n = len(coords)
        distances = np.zeros((n, n))
        
        for i in range(n):
            for j in range(n):
                if i != j:
                    distances[i, j] = self._haversine_distance(
                        coords[i][0], coords[i][1], coords[j][0], coords[j][1]
                    )
        return distances
    
    def _haversine_distance(self, lat1: float, lon1: float, lat2: float, lon2: float) -> float:
        """Calculate haversine distance between two points in km"""
        R = 6371  # Earth's radius in km
        
        lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
        dlat = lat2 - lat1
        dlon = lon2 - lon1
        
        a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
        c = 2 * np.arcsin(np.sqrt(a))
        
        return R * c
    
    def _calculate_genomic_distances(self):
        """Calculate pairwise genomic distances (Hamming distance)"""
        n_samples = len(self.samples)
        self.genomic_distances = np.zeros((n_samples, n_samples))
        
        for i in range(n_samples):
            for j in range(n_samples):
                if i != j:
                    seq1 = self.samples[i].genome_sequence
                    seq2 = self.samples[j].genome_sequence
                    self.genomic_distances[i, j] = self._hamming_distance(seq1, seq2)
    
    def _hamming_distance(self, seq1: str, seq2: str) -> int:
        """Calculate Hamming distance between two sequences"""
        if len(seq1) != len(seq2):
            raise ValueError("Sequences must have equal length")
        return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))
    
    def calculate_transmission_probability(self, source_idx: int, target_idx: int) -> float:
        """
        Calculate transmission probability between two samples based on spatial-temporal proximity
        
        Args:
            source_idx: Index of source sample
            target_idx: Index of target sample
            
        Returns:
            Transmission probability
        """
        spatial_dist = self.spatial_distances[source_idx, target_idx]
        temporal_dist = self.temporal_distances[source_idx, target_idx]
        
        # Spatial decay function (exponential decay with distance)
        spatial_factor = np.exp(-spatial_dist / self.spatial_decay_distance)
        
        # Temporal decay function (exponential decay with time)
        temporal_factor = np.exp(-temporal_dist / self.temporal_decay_time)
        
        # Combined transmission probability
        transmission_prob = self.transmission_rate * spatial_factor * temporal_factor
        
        return min(transmission_prob, 1.0)  # Cap at 1.0
    
    def calculate_evolution_probability(self, genome_a_idx: int, genome_b_idx: int, 
                                     time_period: int) -> float:
        """
        Calculate probability of genome A evolving to genome B over a given time period
        
        Args:
            genome_a_idx: Index of source genome
            genome_b_idx: Index of target genome
            time_period: Time period in days
            
        Returns:
            Evolution probability
        """
        seq_a = self.samples[genome_a_idx].genome_sequence
        seq_b = self.samples[genome_b_idx].genome_sequence
        
        if len(seq_a) != len(seq_b):
            raise ValueError("Sequences must have equal length")
        
        genome_length = len(seq_a)
        differences = self._hamming_distance(seq_a, seq_b)
        
        # Mutation probability
        mutation_prob = self._calculate_mutation_probability(differences, genome_length, time_period)
        
        # Recombination probability
        recombination_prob = self._calculate_recombination_probability(
            genome_a_idx, genome_b_idx, time_period
        )
        
        # Combined evolution probability
        evolution_prob = mutation_prob * recombination_prob
        
        return min(evolution_prob, 1.0)
    
    def _calculate_mutation_probability(self, differences: int, genome_length: int, 
                                      time_period: int) -> float:
        """Calculate mutation probability using Poisson model"""
        expected_mutations = self.mutation_rate * genome_length * time_period
        
        # Poisson probability of observing exactly 'differences' mutations
        mutation_prob = stats.poisson.pmf(differences, expected_mutations)
        
        return mutation_prob
    
    def _calculate_recombination_probability(self, genome_a_idx: int, genome_b_idx: int, 
                                           time_period: int) -> float:
        """Calculate recombination probability based on available donor sequences"""
        # Find potential donor sequences (other samples collected around the same time/location)
        potential_donors = []
        
        for i, sample in enumerate(self.samples):
            if i != genome_a_idx and i != genome_b_idx:
                spatial_dist = self.spatial_distances[genome_a_idx, i]
                temporal_dist = abs(
                    (self.samples[genome_a_idx].collection_date - sample.collection_date).days
                )
                
                # Consider samples within 50 km and 60 days as potential donors
                if spatial_dist <= 50 and temporal_dist <= 60:
                    potential_donors.append(i)
        
        if not potential_donors:
            return 1.0  # No recombination, pure mutation
        
        # Calculate recombination probability based on sequence similarity to donors
        seq_a = self.samples[genome_a_idx].genome_sequence
        seq_b = self.samples[genome_b_idx].genome_sequence
        genome_length = len(seq_a)
        
        # Expected recombination events
        expected_recombinations = self.recombination_rate * genome_length * time_period
        
        # Calculate how well donors can explain the differences
        best_explanation_score = 0
        
        for donor_idx in potential_donors:
            donor_seq = self.samples[donor_idx].genome_sequence
            explanation_score = self._calculate_recombination_explanation(
                seq_a, seq_b, donor_seq
            )
            best_explanation_score = max(best_explanation_score, explanation_score)
        
        # Recombination probability increases with better donor explanations
        recombination_prob = np.exp(-expected_recombinations) + \
                           (1 - np.exp(-expected_recombinations)) * best_explanation_score
        
        return recombination_prob
    
    def _calculate_recombination_explanation(self, seq_a: str, seq_b: str, 
                                           donor_seq: str) -> float:
        """Calculate how well a donor sequence explains differences between A and B"""
        if len(seq_a) != len(seq_b) or len(seq_a) != len(donor_seq):
            return 0.0
        
        explained_positions = 0
        total_differences = 0
        
        for i, (a, b, d) in enumerate(zip(seq_a, seq_b, donor_seq)):
            if a != b:
                total_differences += 1
                if b == d:  # Donor explains this position
                    explained_positions += 1
        
        if total_differences == 0:
            return 1.0
        
        return explained_positions / total_differences
    
    def build_transmission_network(self, probability_threshold: float = 0.01) -> nx.DiGraph:
        """
        Build a transmission network based on calculated probabilities
        
        Args:
            probability_threshold: Minimum probability to include edge
            
        Returns:
            NetworkX directed graph representing transmission network
        """
        G = nx.DiGraph()
        
        # Add nodes
        for i, sample in enumerate(self.samples):
            G.add_node(i, 
                      sample_id=sample.sample_id,
                      village=sample.village,
                      coordinates=sample.coordinates,
                      collection_date=sample.collection_date)
        
        # Add edges based on transmission probabilities
        for i in range(len(self.samples)):
            for j in range(len(self.samples)):
                if i != j:
                    trans_prob = self.calculate_transmission_probability(i, j)
                    if trans_prob >= probability_threshold:
                        G.add_edge(i, j, weight=trans_prob)
        
        return G
    
    def simulate_transmission_dynamics(self, n_steps: int = 100, 
                                     initial_infected: Optional[List[int]] = None) -> Dict:
        """
        Simulate transmission dynamics over time
        
        Args:
            n_steps: Number of simulation steps
            initial_infected: Initial infected sample indices
            
        Returns:
            Dictionary containing simulation results
        """
        if initial_infected is None:
            initial_infected = [0]  # Start with first sample
        
        n_samples = len(self.samples)
        infected_status = np.zeros((n_steps, n_samples), dtype=bool)
        
        # Initialize
        for idx in initial_infected:
            infected_status[0, idx] = True
        
        # Simulate transmission
        for step in range(1, n_steps):
            infected_status[step] = infected_status[step-1].copy()
            
            for i in range(n_samples):
                if infected_status[step-1, i]:  # If infected
                    for j in range(n_samples):
                        if not infected_status[step, j]:  # If susceptible
                            trans_prob = self.calculate_transmission_probability(i, j)
                            if np.random.random() < trans_prob:
                                infected_status[step, j] = True
        
        return {
            'infected_status': infected_status,
            'total_infected': np.sum(infected_status, axis=1),
            'infection_times': self._calculate_infection_times(infected_status)
        }
    
    def _calculate_infection_times(self, infected_status: np.ndarray) -> np.ndarray:
        """Calculate first infection time for each sample"""
        infection_times = np.full(infected_status.shape[1], -1)
        
        for i in range(infected_status.shape[1]):
            infected_steps = np.where(infected_status[:, i])[0]
            if len(infected_steps) > 0:
                infection_times[i] = infected_steps[0]
        
        return infection_times

class TransmissionVisualizer:
    """Visualization tools for transmission dynamics"""
    
    def __init__(self, model: TransmissionDynamicsModel):
        self.model = model
    
    def plot_spatial_distribution(self, figsize: Tuple[int, int] = (12, 8)):
        """Plot spatial distribution of samples"""
        fig, axes = plt.subplots(1, 2, figsize=figsize)
        
        # Extract coordinates and villages
        coords = np.array([s.coordinates for s in self.model.samples])
        villages = [s.village for s in self.model.samples]
        
        # Plot 1: Sample locations by village
        village_colors = {'Village_A': 'red', 'Village_B': 'blue', 'Village_C': 'green'}
        
        for village in set(villages):
            mask = np.array(villages) == village
            axes[0].scatter(coords[mask, 1], coords[mask, 0], 
                          c=village_colors.get(village, 'gray'), 
                          label=village, s=100, alpha=0.7)
        
        axes[0].set_xlabel('Longitude')
        axes[0].set_ylabel('Latitude')
        axes[0].set_title('Sample Locations by Village')
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)
        
        # Plot 2: Spatial distance heatmap
        im = axes[1].imshow(self.model.spatial_distances, cmap='viridis')
        axes[1].set_title('Spatial Distance Matrix (km)')
        axes[1].set_xlabel('Sample Index')
        axes[1].set_ylabel('Sample Index')
        plt.colorbar(im, ax=axes[1])
        
        plt.tight_layout()
        return fig
    
    def plot_genomic_relationships(self, figsize: Tuple[int, int] = (12, 8)):
        """Plot genomic relationships"""
        fig, axes = plt.subplots(1, 2, figsize=figsize)
        
        # Plot 1: Genomic distance heatmap
        im1 = axes[0].imshow(self.model.genomic_distances, cmap='plasma')
        axes[0].set_title('Genomic Distance Matrix')
        axes[0].set_xlabel('Sample Index')
        axes[0].set_ylabel('Sample Index')
        plt.colorbar(im1, ax=axes[0])
        
        # Plot 2: Genomic vs spatial distance correlation
        spatial_flat = self.model.spatial_distances[np.triu_indices_from(
            self.model.spatial_distances, k=1)]
        genomic_flat = self.model.genomic_distances[np.triu_indices_from(
            self.model.genomic_distances, k=1)]
        
        axes[1].scatter(spatial_flat, genomic_flat, alpha=0.6)
        axes[1].set_xlabel('Spatial Distance (km)')
        axes[1].set_ylabel('Genomic Distance')
        axes[1].set_title('Spatial vs Genomic Distance')
        
        # Add correlation coefficient
        corr_coef = np.corrcoef(spatial_flat, genomic_flat)[0, 1]
        axes[1].text(0.05, 0.95, f'r = {corr_coef:.3f}', 
                    transform=axes[1].transAxes, fontsize=12,
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
        
        plt.tight_layout()
        return fig
    
    def plot_transmission_network(self, probability_threshold: float = 0.01, 
                                figsize: Tuple[int, int] = (12, 10)):
        """Plot transmission network"""
        G = self.model.build_transmission_network(probability_threshold)
        
        fig, ax = plt.subplots(figsize=figsize)
        
        # Create layout based on coordinates
        pos = {}
        for node in G.nodes():
            coords = self.model.samples[node].coordinates
            pos[node] = (coords[1], coords[0])  # (lon, lat)
        
        # Color nodes by village
        villages = [self.model.samples[i].village for i in G.nodes()]
        village_colors = {'Village_A': 'red', 'Village_B': 'blue', 'Village_C': 'green'}
        node_colors = [village_colors.get(v, 'gray') for v in villages]
        
        # Draw network
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, 
                              node_size=300, alpha=0.8, ax=ax)
        
        # Draw edges with thickness proportional to probability
        edges = G.edges(data=True)
        edge_weights = [d['weight'] for _, _, d in edges]
        nx.draw_networkx_edges(G, pos, width=[w*5 for w in edge_weights], 
                              alpha=0.6, edge_color='gray', ax=ax)
        
        # Add labels
        labels = {i: self.model.samples[i].sample_id for i in G.nodes()}
        nx.draw_networkx_labels(G, pos, labels, font_size=8, ax=ax)
        
        ax.set_title('Transmission Network\n(Edge thickness ∝ transmission probability)')
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        
        # Add legend
        legend_elements = [plt.Line2D([0], [0], marker='o', color='w', 
                                    markerfacecolor=color, markersize=10, label=village)
                          for village, color in village_colors.items()]
        ax.legend(handles=legend_elements, loc='upper right')
        
        return fig
    
    def plot_simulation_results(self, simulation_results: Dict, 
                              figsize: Tuple[int, int] = (12, 8)):
        """Plot simulation results"""
        fig, axes = plt.subplots(2, 2, figsize=figsize)
        
        # Plot 1: Infection spread over time
        axes[0, 0].plot(simulation_results['total_infected'])
        axes[0, 0].set_xlabel('Time Step')
        axes[0, 0].set_ylabel('Total Infected')
        axes[0, 0].set_title('Infection Spread Over Time')
        axes[0, 0].grid(True, alpha=0.3)
        
        # Plot 2: Infection status heatmap
        im = axes[0, 1].imshow(simulation_results['infected_status'].T, 
                              aspect='auto', cmap='Reds')
        axes[0, 1].set_xlabel('Time Step')
        axes[0, 1].set_ylabel('Sample Index')
        axes[0, 1].set_title('Infection Status Over Time')
        
        # Plot 3: Infection times histogram
        infection_times = simulation_results['infection_times']
        valid_times = infection_times[infection_times >= 0]
        axes[1, 0].hist(valid_times, bins=20, alpha=0.7)
        axes[1, 0].set_xlabel('Infection Time')
        axes[1, 0].set_ylabel('Frequency')
        axes[1, 0].set_title('Distribution of Infection Times')
        
        # Plot 4: Village-wise infection progression
        villages = [s.village for s in self.model.samples]
        village_names = list(set(villages))
        
        for village in village_names:
            village_indices = [i for i, v in enumerate(villages) if v == village]
            village_infections = np.sum(simulation_results['infected_status'][:, village_indices], axis=1)
            axes[1, 1].plot(village_infections, label=village)
        
        axes[1, 1].set_xlabel('Time Step')
        axes[1, 1].set_ylabel('Infected in Village')
        axes[1, 1].set_title('Infection by Village')
        axes[1, 1].legend()
        axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        return fig

def create_sample_data() -> List[GenomicSample]:
    """Create sample genomic data for 3 Tanzanian villages"""
    np.random.seed(42)
    
    # Village coordinates (approximate locations in Tanzania)
    villages = {
        'Village_A': (-6.1659, 35.7516),  # Near Dodoma
        'Village_B': (-6.7924, 39.2083),  # Near Dar es Salaam  
        'Village_C': (-3.3869, 36.6830)   # Near Arusha
    }
    
    samples = []
    base_date = datetime(2024, 1, 1)
    
    # Generate reference genome (1000 bp for demonstration)
    reference_genome = ''.join(np.random.choice(['A', 'T', 'G', 'C'], 1000))
    
    sample_id = 1
    for village, (lat, lon) in villages.items():
        # Generate 10 samples per village
        for i in range(10):
            # Add some spatial noise around village center
            sample_lat = lat + np.random.normal(0, 0.01)
            sample_lon = lon + np.random.normal(0, 0.01)
            
            # Add temporal variation (samples collected over 6 months)
            collection_date = base_date + timedelta(days=np.random.randint(0, 180))
            
            # Generate genome with mutations from reference
            genome = list(reference_genome)
            n_mutations = np.random.poisson(5)  # Average 5 mutations
            
            for _ in range(n_mutations):
                pos = np.random.randint(0, len(genome))
                genome[pos] = np.random.choice(['A', 'T', 'G', 'C'])
            
            # Add some village-specific mutations for clustering
            if village == 'Village_B':
                # Village B has some specific mutations
                for pos in [100, 200, 300]:
                    if pos < len(genome):
                        genome[pos] = 'T'
            elif village == 'Village_C':
                # Village C has different specific mutations
                for pos in [150, 250, 350]:
                    if pos < len(genome):
                        genome[pos] = 'G'
            
            sample = GenomicSample(
                sample_id=f"{village}_{i+1:02d}",
                genome_sequence=''.join(genome),
                village=village,
                coordinates=(sample_lat, sample_lon),
                collection_date=collection_date,
                individual_id=f"IND_{sample_id:03d}"
            )
            
            samples.append(sample)
            sample_id += 1
    
    return samples

def main():
    """Main function demonstrating the transmission dynamics model"""
    print("Genomic Transmission Dynamics Model for Tanzanian Villages")
    print("=" * 60)
    
    # Create sample data
    print("1. Creating sample genomic data...")
    samples = create_sample_data()
    print(f"   Created {len(samples)} samples from 3 villages")
    
    # Initialize model
    print("2. Initializing transmission dynamics model...")
    model = TransmissionDynamicsModel(samples)
    
    # Example: Calculate evolution probability
    print("3. Calculating evolution probabilities...")
    genome_a_idx = 0
    genome_b_idx = 5
    time_period = 30  # 30 days
    
    evolution_prob = model.calculate_evolution_probability(
        genome_a_idx, genome_b_idx, time_period
    )
    
    print(f"   Probability of genome {samples[genome_a_idx].sample_id} "
          f"evolving to {samples[genome_b_idx].sample_id} "
          f"over {time_period} days: {evolution_prob:.6f}")
    
    # Calculate transmission probability
    trans_prob = model.calculate_transmission_probability(genome_a_idx, genome_b_idx)
    print(f"   Transmission probability: {trans_prob:.6f}")
    
    # Build transmission network
    print("4. Building transmission network...")
    network = model.build_transmission_network(probability_threshold=0.001)
    print(f"   Network has {network.number_of_nodes()} nodes and {network.number_of_edges()} edges")
    
    # Run simulation
    print("5. Running transmission simulation...")
    simulation_results = model.simulate_transmission_dynamics(n_steps=50, initial_infected=[0])
    final_infected = simulation_results['total_infected'][-1]
    print(f"   Final number of infected: {final_infected}/{len(samples)}")
    
    # Create visualizations
    print("6. Creating visualizations...")
    visualizer = TransmissionVisualizer(model)
    
    # Spatial distribution
    fig1 = visualizer.plot_spatial_distribution()
    fig1.savefig('/workspace/spatial_distribution.png', dpi=300, bbox_inches='tight')
    print("   Saved: spatial_distribution.png")
    
    # Genomic relationships
    fig2 = visualizer.plot_genomic_relationships()
    fig2.savefig('/workspace/genomic_relationships.png', dpi=300, bbox_inches='tight')
    print("   Saved: genomic_relationships.png")
    
    # Transmission network
    fig3 = visualizer.plot_transmission_network()
    fig3.savefig('/workspace/transmission_network.png', dpi=300, bbox_inches='tight')
    print("   Saved: transmission_network.png")
    
    # Simulation results
    fig4 = visualizer.plot_simulation_results(simulation_results)
    fig4.savefig('/workspace/simulation_results.png', dpi=300, bbox_inches='tight')
    print("   Saved: simulation_results.png")
    
    print("\nAnalysis complete! Check the generated plots for insights into transmission dynamics.")
    
    # Summary statistics
    print("\nSummary Statistics:")
    print("-" * 30)
    
    villages = list(set([s.village for s in samples]))
    for village in villages:
        village_samples = [s for s in samples if s.village == village]
        print(f"{village}: {len(village_samples)} samples")
    
    # Calculate average distances
    spatial_distances = model.spatial_distances[np.triu_indices_from(model.spatial_distances, k=1)]
    genomic_distances = model.genomic_distances[np.triu_indices_from(model.genomic_distances, k=1)]
    
    print(f"\nAverage spatial distance: {np.mean(spatial_distances):.2f} km")
    print(f"Average genomic distance: {np.mean(genomic_distances):.2f} substitutions")
    print(f"Spatial-genomic correlation: {np.corrcoef(spatial_distances, genomic_distances)[0,1]:.3f}")

if __name__ == "__main__":
    main()