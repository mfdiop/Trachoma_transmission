import argparse
import csv
from typing import Dict, Tuple, List
from datetime import datetime

from transmit_model.io import read_fasta, read_metadata_csv, get_sequence_length
from transmit_model.model import (
	hamming_distance,
	RecombinationMutationParams,
	SpaceTimeKernelParams,
	transmission_log_likelihood,
)


def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(description="Compute pairwise transmission log-likelihoods with recombination and spatiotemporal kernel")
	parser.add_argument("--fasta", required=True, help="Path to aligned FASTA of genomes")
	parser.add_argument("--meta", required=True, help="Path to metadata CSV with sample_id,collection_date,[latitude,longitude,village]")
	parser.add_argument("--mu", type=float, default=1e-6, help="Mutation rate per site per day (default 1e-6)")
	parser.add_argument("--r", type=float, default=1e-3, help="Recombination rate per genome per day (default 1e-3)")
	parser.add_argument("--length-scale-km", type=float, default=10.0, help="Spatial kernel length scale in km")
	parser.add_argument("--time-scale-days", type=float, default=30.0, help="Temporal kernel time scale in days")
	parser.add_argument("--cluster-bonus", type=float, default=1.0, help="Multiplicative bonus for same-village pairs")
	parser.add_argument("--out", required=True, help="Output CSV path for pairwise log-likelihoods")
	parser.add_argument("--filter-villages", nargs="*", default=None, help="Optional list of village names to include")
	return parser.parse_args()


def main():
	args = parse_args()

	seqs: Dict[str, str] = read_fasta(args.fasta)
	meta = read_metadata_csv(args.meta)
	L = get_sequence_length(seqs)

	params = RecombinationMutationParams(
		mutation_rate_per_site_per_day=args.mu,
		recombination_rate_per_genome_per_day=args.r,
		genome_length=L,
	)
	kernel = SpaceTimeKernelParams(
		length_scale_km=args.length_scale_km,
		time_scale_days=args.time_scale_days,
		cluster_bonus=args.cluster_bonus,
	)

	# Build ordered list of samples present in both FASTA and metadata, optionally filtered
	sample_ids: List[str] = []
	for sid in seqs.keys():
		if sid not in meta:
			continue
		if args.filter_villages is not None:
			v = (meta[sid].get("village") or "").strip()
			if v not in args.filter_villages:
				continue
		sample_ids.append(sid)

	# Pre-extract fields
	def get_latlon(sid: str) -> Tuple[float | None, float | None]:
		lat_raw = meta[sid].get("latitude")
		lon_raw = meta[sid].get("longitude")
		try:
			lat_val = float(lat_raw) if lat_raw not in (None, "") else None
		except ValueError:
			lat_val = None
		try:
			lon_val = float(lon_raw) if lon_raw not in (None, "") else None
		except ValueError:
			lon_val = None
		return lat_val, lon_val

	def get_date(sid: str) -> datetime:
		return meta[sid]["parsed_date"]

	def same_village_bonus(a: str, b: str) -> float:
		va = (meta[a].get("village") or "").strip()
		vb = (meta[b].get("village") or "").strip()
		return kernel.cluster_bonus if va and vb and va == vb else 1.0

	rows: List[Dict[str, object]] = []
	for i, a in enumerate(sample_ids):
		for j, b in enumerate(sample_ids):
			if i == j:
				continue
			seq_a = seqs[a]
			seq_b = seqs[b]
			date_a = get_date(a)
			date_b = get_date(b)
			lat_a, lon_a = get_latlon(a)
			lat_b, lon_b = get_latlon(b)

			# Adjust kernel with per-pair cluster bonus
			pair_kernel = SpaceTimeKernelParams(
				length_scale_km=kernel.length_scale_km,
				time_scale_days=kernel.time_scale_days,
				cluster_bonus=same_village_bonus(a, b),
			)

			log_like, details = transmission_log_likelihood(
				seq_a, seq_b, date_a, date_b, lat_a, lon_a, lat_b, lon_b, params, pair_kernel
			)

			rows.append({
				"source": a,
				"target": b,
				"days": details["days"],
				"space_km": details["space_km"],
				"hamming": details["hamming"],
				"log_like": log_like,
				"log_like_seq": details["log_like_seq"],
				"log_like_mut_only": details["log_like_mut_only"],
				"lam_mut": details["lam_mut"],
				"lam_rec": details["lam_rec"],
			})

	with open(args.out, "w", newline="", encoding="utf-8") as handle:
		fieldnames = [
			"source", "target", "days", "space_km", "hamming", "log_like", "log_like_seq", "log_like_mut_only", "lam_mut", "lam_rec"
		]
		writer = csv.DictWriter(handle, fieldnames=fieldnames)
		writer.writeheader()
		for row in rows:
			writer.writerow(row)


if __name__ == "__main__":
	main()