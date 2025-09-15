from __future__ import annotations

import math
from dataclasses import dataclass
from datetime import datetime
from typing import Dict, Tuple


def hamming_distance(seq_a: str, seq_b: str) -> int:
	if len(seq_a) != len(seq_b):
		raise ValueError("Sequences must be same length for Hamming distance")
	return sum(1 for a, b in zip(seq_a, seq_b) if a != b and a != 'N' and b != 'N')


@dataclass
class RecombinationMutationParams:
	mutation_rate_per_site_per_day: float
	recombination_rate_per_genome_per_day: float
	genome_length: int


@dataclass
class SpaceTimeKernelParams:
	length_scale_km: float
	time_scale_days: float
	cluster_bonus: float = 1.0


EARTH_RADIUS_KM = 6371.0


def haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
	"""
	Great-circle distance between two points on Earth in kilometers.
	"""
	phi1 = math.radians(lat1)
	phi2 = math.radians(lat2)
	dphi = phi2 - phi1
	dlambda = math.radians(lon2 - lon1)
	a = math.sin(dphi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda / 2) ** 2
	c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
	return EARTH_RADIUS_KM * c


def poisson_log_pmf(k: int, lam: float) -> float:
	if lam <= 0:
		return float('-inf') if k > 0 else 0.0
	return k * math.log(lam) - lam - math.lgamma(k + 1)


def logsumexp(values):
	m = max(values)
	return m + math.log(sum(math.exp(v - m) for v in values))


def transmission_log_likelihood(
	seq_a: str,
	seq_b: str,
	date_a: datetime,
	date_b: datetime,
	lat_a: float | None,
	lon_a: float | None,
	lat_b: float | None,
	lon_b: float | None,
	rec_mut_params: RecombinationMutationParams,
	kernel_params: SpaceTimeKernelParams,
	prior_same_cluster: float | None = None,
) -> Tuple[float, Dict[str, float]]:
	"""
	Returns log-likelihood that genome A evolves to genome B over the elapsed time,
	combining mutation and recombination, and weighted by a spatiotemporal kernel.

	We model mutations as Poisson with rate mu * L * dt on the number of nucleotide changes
	observed (Hamming distance ignoring Ns). Recombination events per genome are Poisson with
	mean r * dt. Each recombination draws a random donor from the local diversity, which we
	approximate by allowing a mixture: with probability p_rec_effect the observed distance is
	reduced compared to pure mutation. For computational simplicity, we integrate over a small
	number of recombination counts (0,1,2+) using a collapsed approximation.
	"""
	if len(seq_a) != len(seq_b):
		raise ValueError("Sequences must be same length")
	if rec_mut_params.genome_length != len(seq_a):
		raise ValueError("genome_length in params must match sequence length")

	# Time component
	days = (date_b - date_a).days
	if days < 0:
		return float('-inf'), {"days": days, "space_km": float('nan'), "mut_log_like": float('-inf'), "rec_effect": 0.0, "kernel_log_weight": float('-inf')}

	# Spatial component
	space_km: float
	if None in (lat_a, lon_a, lat_b, lon_b):
		space_km = float('nan')
		space_log_weight = 0.0
	else:
		space_km = haversine_km(float(lat_a), float(lon_a), float(lat_b), float(lon_b))
		space_log_weight = - (space_km / kernel_params.length_scale_km)

	# Temporal kernel
	time_log_weight = - (days / kernel_params.time_scale_days) if kernel_params.time_scale_days > 0 else 0.0

	kernel_log_weight = space_log_weight + time_log_weight + math.log(max(kernel_params.cluster_bonus, 1e-12))

	# Sequence distance
	d = hamming_distance(seq_a, seq_b)
	mu = rec_mut_params.mutation_rate_per_site_per_day
	r = rec_mut_params.recombination_rate_per_genome_per_day
	L = rec_mut_params.genome_length

	# Pure mutation likelihood: Poisson on observed differences
	lam_mut = max(mu * L * days, 1e-12)
	log_like_mut_only = poisson_log_pmf(d, lam_mut)

	# Recombination-adjusted mixture approximation
	# Assume expected number of recombination events: lam_rec = r * days
	lam_rec = max(r * days, 0.0)
	# We approximate effects:
	# - 0 recomb: same as pure mutation
	# - 1 recomb: reduces effective differences by fraction alpha (borrowed from donor proximity)
	# - 2+ recomb: capped reduction beta
	alpha = 0.30  # single recombination reduces effective distance by 30%
	beta = 0.50   # two or more recombinations reduce by 50%
	eps = 1e-12
	p0 = math.exp(-lam_rec)
	p1 = lam_rec * math.exp(-lam_rec)
	p2plus = max(1.0 - p0 - p1, 0.0)

	eff_d_0 = d
	eff_d_1 = max(int(round((1.0 - alpha) * d)), 0)
	eff_d_2 = max(int(round((1.0 - beta) * d)), 0)

	log_mix_components = [
		math.log(max(p0, eps)) + poisson_log_pmf(eff_d_0, lam_mut),
		math.log(max(p1, eps)) + poisson_log_pmf(eff_d_1, lam_mut),
		math.log(max(p2plus, eps)) + poisson_log_pmf(eff_d_2, lam_mut),
	]
	log_like_seq = logsumexp(log_mix_components)

	# Combine with spatiotemporal kernel as weight (acts like a prior)
	log_like_total = log_like_seq + kernel_log_weight

	return log_like_total, {
		"days": float(days),
		"space_km": space_km,
		"hamming": float(d),
		"lam_mut": lam_mut,
		"lam_rec": lam_rec,
		"kernel_log_weight": kernel_log_weight,
		"log_like_mut_only": log_like_mut_only,
		"log_like_seq": log_like_seq,
	}