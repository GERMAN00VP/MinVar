import pysam
import pandas as pd
import glob
import os


def find_file(pattern):
    files = glob.glob(pattern)
    if not files:
        raise FileNotFoundError(f"No file found for pattern: {pattern}")
    return files[0]


def load_gff_regions(gff_path):
    regions = set()

    with open(gff_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                # ignorar líneas malformadas
                continue
            try:
                start = int(parts[3])
                end = int(parts[4])
            except ValueError:
                continue  # saltar si start/end no son números

            regions.update(range(start, end + 1))

    return regions


def parse_gff(gff_path):
    """
    Parse GFF and return list of features
    """
    features = []

    with open(gff_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue  # saltar líneas que no tengan 9 columnas

            chrom = parts[0]
            feature_type = parts[2]
            try:
                start = int(parts[3])
                end = int(parts[4])
            except ValueError:
                continue

            attributes = parts[8]
            attr_dict = {}
            for attr in attributes.split(";"):
                if "=" in attr:
                    k, v = attr.split("=", 1)
                    attr_dict[k] = v

            features.append({
                "start": start,
                "end": end,
                "type": feature_type,
                "id": attr_dict.get("ID", ""),
                "gene": attr_dict.get("gene", attr_dict.get("Name", "")),
                "product": attr_dict.get("product", "")
            })

    return features


def annotate_position(pos, features):
    """
    Return annotation for a genomic position
    """

    hits = []

    for f in features:
        if f["start"] <= pos <= f["end"]:
            hits.append(f)

    if not hits:
        return {
            "Feature_type": "intergenic",
            "Feature_ID": "",
            "Gene_name": "",
            "Product": "",
            "Region_start": "",
            "Region_end": ""
        }

    # si hay múltiples features (ej: CDS + gene), los concatenamos
    return {
        "Feature_type": ",".join([h["type"] for h in hits]),
        "Feature_ID": ",".join([h["id"] for h in hits]),
        "Gene_name": ",".join([h["gene"] for h in hits]),
        "Product": ",".join([h["product"] for h in hits]),
        "Region_start": ",".join([str(h["start"]) for h in hits]),
        "Region_end": ",".join([str(h["end"]) for h in hits])
    }



def analyze_position_all_alleles(bam, positions):
    """
    Optimized single-pass pileup for all alleles in target positions.
    Returns dict of dicts: {pos: {allele: metrics}}
    """
    min_pos = min(positions)
    max_pos = max(positions)
    target_positions = set(positions)
    pileup_cache = {}

    for pileupcolumn in bam.pileup(start=min_pos-1, stop=max_pos, truncate=True):
        pos = pileupcolumn.reference_pos + 1  # 1-based
        if pos not in target_positions:
            continue

        metrics = {}
        for pileupread in pileupcolumn.pileups:
            if pileupread.is_del or pileupread.is_refskip:
                continue
            read = pileupread.alignment
            read_pos = pileupread.query_position
            base = read.query_sequence[read_pos]

            if base not in metrics:
                metrics[base] = {
                    "total": 0,
                    "extreme": 0,
                    "forward": 0,
                    "reverse": 0,
                    "qualities": []
                }

            m = metrics[base]
            m["total"] += 1
            rel_pos = read_pos / read.query_length
            if rel_pos < 0.1 or rel_pos > 0.9:
                m["extreme"] += 1
            if read.is_reverse:
                m["reverse"] += 1
            else:
                m["forward"] += 1
            m["qualities"].append(read.query_qualities[read_pos])

        # store metrics for this position
        pileup_cache[pos] = {
            base: {
                "depth": m["total"],
                "extreme_fraction": m["extreme"] / m["total"],
                "strand_bias": abs(m["forward"] - m["reverse"]) / m["total"],
                "mean_base_quality": sum(m["qualities"]) / len(m["qualities"])
            }
            for base, m in metrics.items() if m["total"] > 0
        }

    return pileup_cache

