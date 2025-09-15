## Transmission likelihood with recombination and spatiotemporal kernel

This repository provides a simple Python CLI to estimate pairwise transmission log-likelihoods between genomes assuming person-to-person spread, with a mutation-plus-recombination model and a spatiotemporal kernel favoring proximity in space and time. It is designed for genomic data from multiple villages (e.g., 3 in Tanzania), but will work for any dataset with coordinates and collection dates.

### Inputs
- FASTA: an aligned nucleotide FASTA. All sequences must have the same length. Headers' first token must match `sample_id` in metadata.
- Metadata CSV with columns:
  - `sample_id` (required)
  - `collection_date` in YYYY-MM-DD (required)
  - `latitude`, `longitude` in decimal degrees (optional but recommended)
  - `village` (optional; if present, pairs from the same village get a configurable bonus)

### Install
No external dependencies beyond Python 3.9+. Clone this folder and run the script directly.

### Usage
```bash
python /workspace/transmit_likelihood.py \
  --fasta /absolute/path/genomes.fasta \
  --meta /absolute/path/metadata.csv \
  --mu 1e-6 \
  --r 1e-3 \
  --length-scale-km 10 \
  --time-scale-days 30 \
  --cluster-bonus 2.0 \
  --out /absolute/path/pairwise_likelihoods.csv
```

Optional village filter:
```bash
python /workspace/transmit_likelihood.py --fasta ... --meta ... --filter-villages VillageA VillageB --out ...
```

### Model
- Mutation: Poisson on number of differences `d` with rate `λ_mut = μ * L * Δt`.
- Recombination: expected events `λ_rec = r * Δt`. We approximate effects via a three-component mixture over recombination counts `k ∈ {0, 1, ≥2}` that reduce the effective distance by fixed fractions (configurable in code as `alpha` and `beta`). This captures that recombination can import locally common alleles, reducing the apparent divergence due to pure mutation.
- Space-time kernel: multiplicative weight `exp(-distance_km/ℓ_space) * exp(-Δt/ℓ_time) * cluster_bonus`, where `cluster_bonus` > 1 for same-village pairs.

The output is a CSV with pairwise rows and columns: `source,target,days,space_km,hamming,log_like,log_like_seq,log_like_mut_only,lam_mut,lam_rec`.

### Notes
- Ensure sequences are aligned; ambiguous bases `N` are ignored in Hamming distance.
- If coordinates are missing, the spatial term defaults to neutral (no penalty or bonus).
- The parameters `alpha` and `beta` in `transmit_model/model.py` can be tuned or replaced by a more detailed recombination tract model if desired.

### Limitations
- This is a pragmatic approximation; it does not reconstruct transmission trees, only pairwise likelihoods.
- Recombination reduction factors are coarse; consider calibrating to your organism.
- Does not model within-host diversity explicitly.
