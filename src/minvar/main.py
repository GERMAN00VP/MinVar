import os
import pandas as pd
import pysam
import logging
import time
from .utils import find_file, analyze_position_all_alleles,load_gff_regions, parse_gff, annotate_position

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] [%(processName)s] %(message)s",
)



def run_minvar(folder, min_freq, gff_path=None, min_depth=20, min_base_quality=25):

    start_time = time.time()
    warnings = []

    logging.info(f"START processing folder: {folder}")

    try:
        bam_path = find_file(os.path.join(folder, "*.bam"))
        alleles_path = find_file(os.path.join(folder, "tables", "*-allAlleles.txt"))
    except FileNotFoundError as e:
        logging.warning(f"Skipping folder {folder}: {e}")
        return {"folder": folder, "status": "FAILED", "warnings": [str(e)]}

    output_dir = os.path.join(folder, "MinVar_results")
    os.makedirs(output_dir, exist_ok=True)

    logging.info(f"Reading allele table: {alleles_path}")
    df = pd.read_csv(alleles_path, sep="\t")
    df = df[df.Frequency != 1]

    features = None
    if gff_path:
        logging.info("Loading GFF annotations")
        regions = load_gff_regions(gff_path)
        features = parse_gff(gff_path)
        df = df[df.Position.isin(regions)]

    logging.info("Filtering candidate minority variants")
    df_min = df[
        (df.Frequency > min_freq) &
        (df.Frequency < 0.5) &
        (df.Count > min_depth)
    ]

    df_min = df_min[
        (df_min.ConfidenceNotMacErr > 0.95) &
        (df_min.PairedUB < 0.05) &
        (df_min.QualityUB < 0.05)
    ]

    df_cons = df[
        (df.Position.isin(df_min.Position)) &
        (df.Allele_Type == "Consensus")
    ]

    bam = pysam.AlignmentFile(bam_path, "rb")

    logging.info(f"Processing pileups for {len(df_min.Position.unique())} positions")
    pileup_cache = analyze_position_all_alleles(bam, df_min.Position.unique())

    results = []
    total_variants = len(df_min)
    logging.info(f"Processing {total_variants} candidate variants")

    for i, (_, row) in enumerate(df_min.iterrows(), start=1):
        if i % 50 == 0 or i == total_variants:
            logging.info(f"Progress: {i}/{total_variants} variants")

        pos = row.Position
        alt = row.Allele

        cons_row = df_cons[df_cons.Position == pos]
        if cons_row.empty:
            warnings.append(f"No consensus for position {pos}")
            continue
        cons = cons_row.iloc[0].Allele

        all_metrics = pileup_cache.get(pos, {})
        alt_metrics = all_metrics.get(alt)
        cons_metrics = all_metrics.get(cons)

        if not alt_metrics or not cons_metrics:
            warnings.append(f"No metrics for position {pos}")
            continue

        status = "PASS"
        if alt_metrics["extreme_fraction"] > 0.5:
            if abs(alt_metrics["extreme_fraction"] - cons_metrics["extreme_fraction"]) > 0.2:
                status = "FAIL_position_bias"
        if alt_metrics["strand_bias"] > 0.8:
            if abs(alt_metrics["strand_bias"] - cons_metrics["strand_bias"]) > 0.2:
                status = "FAIL_strand_bias"
        if alt_metrics["mean_base_quality"] < min_base_quality:
            status = "FAIL_low_quality"

        annotation = annotate_position(pos, features) if features else {
            "Feature_type": "",
            "Feature_ID": "",
            "Gene_name": "",
            "Product": "",
            "Region_start": "",
            "Region_end": ""
        }

        results.append({
            "Position": pos,
            "Alt": alt,
            "Alt_freq": row.Frequency,
            "Consensus": cons,
            "Alt_extreme_fraction": alt_metrics["extreme_fraction"],
            "Consensus_extreme_fraction": cons_metrics["extreme_fraction"],
            "Consensus_strand_bias": cons_metrics["strand_bias"],
            "Alt_strand_bias": alt_metrics["strand_bias"],
            "Mean_base_quality": alt_metrics["mean_base_quality"],
            "Status": status,
            **annotation
        })

    bam.close()

    df_out = pd.DataFrame(results)
    out_file = os.path.join(output_dir, "minority_variants.tsv")
    df_out.to_csv(out_file, sep="\t", index=False)

    write_metrics_description(output_dir)

    elapsed = time.time() - start_time

    if warnings:
        logging.warning(f"Completed {folder} with {len(warnings)} warnings in {elapsed:.2f}s")
        logging.warning("Warnings: " + "; ".join(warnings))
    else:
        logging.info(f"Completed {folder} successfully in {elapsed:.2f}s")

    return {
        "folder": folder,
        "status": "OK" if not warnings else "WARNING",
        "warnings": warnings,
        "time": elapsed
    }


def write_metrics_description(outdir):

    text = """MINVAR OUTPUT METRICS DESCRIPTION

Position:
Genomic coordinate of the variant.

Alt:
Minority allele detected.

Consensus:
Consensus allele at that position.

Alt_freq:
Frequency of the minority allele.

Alt_extreme_fraction:
Fraction of reads supporting the variant located within the first or last 10% of the read length.
High values suggest positional bias (sequencing artifacts at read ends).

Consensus_extreme_fraction:
Same metric for the consensus allele at the same position.
Used as internal control to distinguish technical bias from true variant-specific bias.

Strand_bias:
Absolute imbalance between forward and reverse reads supporting the variant:
|forward - reverse| / total.
Values close to 1 indicate strong strand bias.

Consensus_strand_bias:
Same strand bias metric calculated for the consensus allele.
Used to detect whether strand bias is due to library preparation or specific to the variant.

Mean_base_quality:
Average Phred quality score of bases supporting the variant.
Low values indicate unreliable base calls.

Status:
PASS
→ Variant passes all filters and is considered reliable.

FAIL_position_bias
→ Alt_extreme_fraction > threshold (default 0.5)
→ AND differs significantly from Consensus_extreme_fraction (default difference > 0.2)
Interpretation: variant is enriched at read ends beyond what is observed in the consensus → likely positional artifact.

FAIL_strand_bias
→ Strand_bias > threshold (default 0.8)
→ AND differs significantly from Consensus_strand_bias (default difference > 0.2)
Interpretation: variant shows strand imbalance not explained by overall library bias → likely strand-specific artifact.

FAIL_low_quality
→ Mean_base_quality < threshold (default 25)
Interpretation: supporting bases have low sequencing quality → unreliable variant call.

NOTE:
Filters are applied sequentially. A variant is marked as FAIL if any condition is met.
Consensus-normalized filtering allows distinguishing true biological variation from systematic sequencing biases.
"""

    with open(os.path.join(outdir, "README_metrics.txt"), "w") as f:
        f.write(text)