import csv
from datetime import datetime
from typing import Dict, List, Tuple


def read_fasta(filepath: str) -> Dict[str, str]:
	"""
	Reads an aligned FASTA file and returns a mapping of sample_id to sequence (uppercased).
	Assumes all sequences are the same length. Does not require external libraries.
	"""
	sequences: Dict[str, List[str]] = {}
	current_id: str = ""
	try:
		with open(filepath, "r", encoding="utf-8") as handle:
			for raw_line in handle:
				line = raw_line.strip()
				if not line:
					continue
				if line.startswith(">"):
					current_id = line[1:].strip().split()[0]
					if current_id in sequences:
						raise ValueError(f"Duplicate FASTA record id: {current_id}")
					sequences[current_id] = []
				else:
					if not current_id:
						raise ValueError("FASTA sequence data encountered before any header line starting with '>'")
					sequences[current_id].append(line.upper())
	except FileNotFoundError as exc:
		raise FileNotFoundError(f"FASTA file not found: {filepath}") from exc

	joined: Dict[str, str] = {seq_id: "".join(parts) for seq_id, parts in sequences.items()}
	if not joined:
		raise ValueError(f"No sequences found in FASTA: {filepath}")

	# Validate lengths
	lengths = {len(seq) for seq in joined.values()}
	if len(lengths) != 1:
		raise ValueError(f"FASTA sequences must be aligned with equal lengths, found lengths: {sorted(lengths)}")

	return joined


def read_metadata_csv(filepath: str) -> Dict[str, Dict[str, str]]:
	"""
	Reads metadata CSV with required columns: sample_id, collection_date (YYYY-MM-DD).
	Optional columns: latitude, longitude, village.
	Returns mapping from sample_id to its metadata dict including parsed_date (datetime).
	"""
	records: Dict[str, Dict[str, str]] = {}
	try:
		with open(filepath, "r", encoding="utf-8") as handle:
			reader = csv.DictReader(handle)
			required = {"sample_id", "collection_date"}
			missing = [c for c in required if c not in reader.fieldnames]
			if missing:
				raise ValueError(f"Metadata CSV missing required columns: {missing}")
			for row in reader:
				sample_id = (row.get("sample_id") or "").strip()
				if not sample_id:
					continue
				date_str = (row.get("collection_date") or "").strip()
				if not date_str:
					raise ValueError(f"Missing collection_date for sample_id={sample_id}")
				try:
					parsed_date = datetime.strptime(date_str, "%Y-%m-%d")
				except ValueError as exc:
					raise ValueError(
						f"Invalid date format for sample_id={sample_id}, collection_date='{date_str}', expected YYYY-MM-DD"
					) from exc
				row["parsed_date"] = parsed_date
				# Normalize optional fields
				for key in ["latitude", "longitude", "village"]:
					if key in row and isinstance(row[key], str):
						row[key] = row[key].strip()
				records[sample_id] = row
	except FileNotFoundError as exc:
		raise FileNotFoundError(f"Metadata CSV not found: {filepath}") from exc

	return records


def get_sequence_length(sequences: Dict[str, str]) -> int:
	lengths = {len(seq) for seq in sequences.values()}
	if not lengths:
		return 0
	if len(lengths) != 1:
		raise ValueError(f"Sequences are not all the same length: {sorted(lengths)}")
	return next(iter(lengths))