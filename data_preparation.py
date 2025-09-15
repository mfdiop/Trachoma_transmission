#!/usr/bin/env python3
"""
Data Preparation Utilities for Genomic Transmission Dynamics

This module provides utilities to prepare your genomic data for analysis,
including sequence alignment, metadata processing, and data validation.
"""

import pandas as pd
import numpy as np
from datetime import datetime
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass
import os
import warnings

from genomic_transmission_dynamics import GenomicSample

@dataclass
class RawGenomicData:
    """Raw genomic data before processing"""
    sample_id: str
    sequence: str
    village: str
    latitude: Optional[float] = None
    longitude: Optional[float] = None
    collection_date: Optional[str] = None
    individual_id: Optional[str] = None
    additional_metadata: Optional[Dict] = None

class DataPreprocessor:
    """Preprocessing utilities for genomic data"""
    
    def __init__(self, reference_sequence: Optional[str] = None):
        """
        Initialize preprocessor
        
        Args:
            reference_sequence: Reference sequence for alignment (optional)
        """
        self.reference_sequence = reference_sequence
        self.processed_samples = []
    
    def load_sequences_from_fasta(self, fasta_file: str) -> Dict[str, str]:
        """
        Load sequences from FASTA file
        
        Args:
            fasta_file: Path to FASTA file
            
        Returns:
            Dictionary mapping sample IDs to sequences
        """
        sequences = {}
        current_id = None
        current_seq = []
        
        try:
            with open(fasta_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_id is not None:
                            sequences[current_id] = ''.join(current_seq)
                        current_id = line[1:]  # Remove '>'
                        current_seq = []
                    else:
                        current_seq.append(line.upper())
                
                # Add last sequence
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq)
                    
        except FileNotFoundError:
            print(f"Error: FASTA file '{fasta_file}' not found")
            return {}
        
        print(f"Loaded {len(sequences)} sequences from {fasta_file}")
        return sequences
    
    def load_metadata_from_csv(self, csv_file: str) -> pd.DataFrame:
        """
        Load metadata from CSV file
        
        Expected columns:
        - sample_id: Unique sample identifier
        - village: Village name
        - latitude: GPS latitude
        - longitude: GPS longitude  
        - collection_date: Date in YYYY-MM-DD format
        - individual_id: Individual identifier (optional)
        
        Args:
            csv_file: Path to CSV file
            
        Returns:
            DataFrame with metadata
        """
        try:
            df = pd.read_csv(csv_file)
            
            # Validate required columns
            required_cols = ['sample_id', 'village', 'latitude', 'longitude', 'collection_date']
            missing_cols = [col for col in required_cols if col not in df.columns]
            
            if missing_cols:
                print(f"Warning: Missing required columns: {missing_cols}")
            
            print(f"Loaded metadata for {len(df)} samples from {csv_file}")
            return df
            
        except FileNotFoundError:
            print(f"Error: Metadata file '{csv_file}' not found")
            return pd.DataFrame()
    
    def validate_sequences(self, sequences: Dict[str, str]) -> Dict[str, str]:
        """
        Validate and clean genomic sequences
        
        Args:
            sequences: Dictionary of sequences
            
        Returns:
            Dictionary of validated sequences
        """
        validated = {}
        valid_bases = set('ATCGN-')
        
        for sample_id, sequence in sequences.items():
            # Remove whitespace
            clean_seq = sequence.replace(' ', '').replace('\n', '').replace('\t', '')
            
            # Check for valid DNA bases
            invalid_bases = set(clean_seq) - valid_bases
            if invalid_bases:
                print(f"Warning: Sample {sample_id} contains invalid bases: {invalid_bases}")
                # Replace invalid bases with N
                for base in invalid_bases:
                    clean_seq = clean_seq.replace(base, 'N')
            
            # Check sequence length
            if len(clean_seq) == 0:
                print(f"Warning: Sample {sample_id} has empty sequence")
                continue
            
            validated[sample_id] = clean_seq
        
        return validated
    
    def align_sequences(self, sequences: Dict[str, str], 
                       method: str = 'simple') -> Dict[str, str]:
        """
        Align sequences (simple implementation)
        
        For production use, consider using proper alignment tools like MAFFT, MUSCLE, etc.
        
        Args:
            sequences: Dictionary of sequences
            method: Alignment method ('simple' or 'reference')
            
        Returns:
            Dictionary of aligned sequences
        """
        if method == 'simple':
            # Simple alignment: pad to same length
            max_length = max(len(seq) for seq in sequences.values())
            
            aligned = {}
            for sample_id, sequence in sequences.items():
                # Pad with N's to reach max length
                padded_seq = sequence + 'N' * (max_length - len(sequence))
                aligned[sample_id] = padded_seq
            
            print(f"Simple alignment completed. Aligned length: {max_length}")
            return aligned
            
        elif method == 'reference' and self.reference_sequence:
            # Reference-based alignment (simplified)
            ref_length = len(self.reference_sequence)
            
            aligned = {}
            for sample_id, sequence in sequences.items():
                if len(sequence) == ref_length:
                    aligned[sample_id] = sequence
                elif len(sequence) < ref_length:
                    # Pad to reference length
                    aligned[sample_id] = sequence + 'N' * (ref_length - len(sequence))
                else:
                    # Truncate to reference length
                    aligned[sample_id] = sequence[:ref_length]
                    print(f"Warning: Truncated sequence {sample_id} to reference length")
            
            print(f"Reference alignment completed. Aligned length: {ref_length}")
            return aligned
        
        else:
            print("Warning: Unknown alignment method or missing reference. Using simple alignment.")
            return self.align_sequences(sequences, method='simple')
    
    def parse_date(self, date_str: str) -> Optional[datetime]:
        """
        Parse date string to datetime object
        
        Args:
            date_str: Date string in various formats
            
        Returns:
            datetime object or None if parsing fails
        """
        if pd.isna(date_str) or date_str == '':
            return None
        
        # Try different date formats
        date_formats = [
            '%Y-%m-%d',
            '%Y/%m/%d', 
            '%d-%m-%Y',
            '%d/%m/%Y',
            '%m-%d-%Y',
            '%m/%d/%Y',
            '%Y-%m-%d %H:%M:%S',
            '%Y/%m/%d %H:%M:%S'
        ]
        
        for fmt in date_formats:
            try:
                return datetime.strptime(str(date_str), fmt)
            except ValueError:
                continue
        
        print(f"Warning: Could not parse date '{date_str}'")
        return None
    
    def create_genomic_samples(self, sequences: Dict[str, str], 
                              metadata: pd.DataFrame) -> List[GenomicSample]:
        """
        Create GenomicSample objects from sequences and metadata
        
        Args:
            sequences: Dictionary of aligned sequences
            metadata: DataFrame with sample metadata
            
        Returns:
            List of GenomicSample objects
        """
        samples = []
        
        for sample_id, sequence in sequences.items():
            # Find metadata for this sample
            sample_metadata = metadata[metadata['sample_id'] == sample_id]
            
            if sample_metadata.empty:
                print(f"Warning: No metadata found for sample {sample_id}")
                continue
            
            row = sample_metadata.iloc[0]
            
            # Parse coordinates
            try:
                latitude = float(row['latitude'])
                longitude = float(row['longitude'])
            except (ValueError, KeyError):
                print(f"Warning: Invalid coordinates for sample {sample_id}")
                continue
            
            # Parse date
            collection_date = self.parse_date(row['collection_date'])
            if collection_date is None:
                print(f"Warning: Invalid date for sample {sample_id}")
                continue
            
            # Create sample
            sample = GenomicSample(
                sample_id=sample_id,
                genome_sequence=sequence,
                village=str(row['village']),
                coordinates=(latitude, longitude),
                collection_date=collection_date,
                individual_id=str(row.get('individual_id', sample_id))
            )
            
            samples.append(sample)
        
        print(f"Created {len(samples)} GenomicSample objects")
        return samples
    
    def process_data(self, fasta_file: str, metadata_file: str, 
                    alignment_method: str = 'simple') -> List[GenomicSample]:
        """
        Complete data processing pipeline
        
        Args:
            fasta_file: Path to FASTA file with sequences
            metadata_file: Path to CSV file with metadata
            alignment_method: Method for sequence alignment
            
        Returns:
            List of processed GenomicSample objects
        """
        print("Starting data processing pipeline...")
        
        # Load data
        sequences = self.load_sequences_from_fasta(fasta_file)
        metadata = self.load_metadata_from_csv(metadata_file)
        
        if not sequences or metadata.empty:
            print("Error: Could not load sequences or metadata")
            return []
        
        # Validate sequences
        sequences = self.validate_sequences(sequences)
        
        # Align sequences
        aligned_sequences = self.align_sequences(sequences, method=alignment_method)
        
        # Create samples
        samples = self.create_genomic_samples(aligned_sequences, metadata)
        
        print("Data processing completed!")
        return samples

class DataValidator:
    """Validate processed genomic data"""
    
    @staticmethod
    def validate_samples(samples: List[GenomicSample]) -> Dict[str, List[str]]:
        """
        Validate list of genomic samples
        
        Args:
            samples: List of GenomicSample objects
            
        Returns:
            Dictionary of validation issues
        """
        issues = {
            'errors': [],
            'warnings': []
        }
        
        if not samples:
            issues['errors'].append("No samples provided")
            return issues
        
        # Check sequence lengths
        seq_lengths = [len(s.genome_sequence) for s in samples]
        if len(set(seq_lengths)) > 1:
            issues['warnings'].append(f"Sequences have different lengths: {set(seq_lengths)}")
        
        # Check for duplicate sample IDs
        sample_ids = [s.sample_id for s in samples]
        duplicates = set([x for x in sample_ids if sample_ids.count(x) > 1])
        if duplicates:
            issues['errors'].append(f"Duplicate sample IDs: {duplicates}")
        
        # Check coordinate validity
        for sample in samples:
            lat, lon = sample.coordinates
            if not (-90 <= lat <= 90):
                issues['errors'].append(f"Invalid latitude for {sample.sample_id}: {lat}")
            if not (-180 <= lon <= 180):
                issues['errors'].append(f"Invalid longitude for {sample.sample_id}: {lon}")
        
        # Check date validity
        dates = [s.collection_date for s in samples]
        min_date = min(dates)
        max_date = max(dates)
        date_range = (max_date - min_date).days
        
        if date_range > 365 * 5:  # More than 5 years
            issues['warnings'].append(f"Large date range: {date_range} days")
        
        # Check village distribution
        villages = [s.village for s in samples]
        village_counts = pd.Series(villages).value_counts()
        
        if len(village_counts) < 2:
            issues['warnings'].append("Only one village represented")
        
        single_sample_villages = village_counts[village_counts == 1].index.tolist()
        if single_sample_villages:
            issues['warnings'].append(f"Villages with single samples: {single_sample_villages}")
        
        return issues
    
    @staticmethod
    def print_validation_report(issues: Dict[str, List[str]]):
        """Print validation report"""
        print("\n=== Data Validation Report ===")
        
        if issues['errors']:
            print("ERRORS:")
            for error in issues['errors']:
                print(f"  ❌ {error}")
        
        if issues['warnings']:
            print("WARNINGS:")
            for warning in issues['warnings']:
                print(f"  ⚠️  {warning}")
        
        if not issues['errors'] and not issues['warnings']:
            print("✅ All validation checks passed!")

def create_example_data_files():
    """Create example data files for testing"""
    
    print("Creating example data files...")
    
    # Create example FASTA file
    fasta_content = """>Sample_001
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>Sample_002  
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>Sample_003
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
"""
    
    with open('/workspace/example_sequences.fasta', 'w') as f:
        f.write(fasta_content)
    
    # Create example metadata CSV
    metadata_df = pd.DataFrame({
        'sample_id': ['Sample_001', 'Sample_002', 'Sample_003'],
        'village': ['Village_A', 'Village_A', 'Village_B'],
        'latitude': [-6.1659, -6.1660, -6.7924],
        'longitude': [35.7516, 35.7517, 39.2083],
        'collection_date': ['2024-01-15', '2024-01-20', '2024-02-01'],
        'individual_id': ['IND_001', 'IND_002', 'IND_003']
    })
    
    metadata_df.to_csv('/workspace/example_metadata.csv', index=False)
    
    print("Created example_sequences.fasta and example_metadata.csv")

if __name__ == "__main__":
    print("Genomic Data Preparation Utilities")
    print("=" * 40)
    
    # Create example files
    create_example_data_files()
    
    # Demonstrate data processing
    preprocessor = DataPreprocessor()
    
    try:
        samples = preprocessor.process_data(
            fasta_file='/workspace/example_sequences.fasta',
            metadata_file='/workspace/example_metadata.csv'
        )
        
        # Validate samples
        validator = DataValidator()
        issues = validator.validate_samples(samples)
        validator.print_validation_report(issues)
        
        print(f"\nProcessed {len(samples)} samples successfully!")
        
        # Show sample details
        if samples:
            print("\nSample details:")
            for sample in samples[:3]:  # Show first 3
                print(f"  {sample.sample_id}: {sample.village}, "
                      f"{len(sample.genome_sequence)} bp, "
                      f"{sample.collection_date.strftime('%Y-%m-%d')}")
        
    except Exception as e:
        print(f"Error during processing: {e}")
    
    print("\nTo use with your data:")
    print("1. Prepare FASTA file with your sequences")
    print("2. Prepare CSV file with metadata (sample_id, village, latitude, longitude, collection_date)")
    print("3. Run: samples = preprocessor.process_data('your_sequences.fasta', 'your_metadata.csv')")
    print("4. Use samples with TransmissionDynamicsModel")