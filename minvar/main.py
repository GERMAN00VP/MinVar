import os
import pandas as pd
import pysam
from .utils import find_file, analyze_variant_position,load_gff_regions, parse_gff, annotate_position


def run_minvar(folder, min_freq, gff_path=None):

    bam_path = find_file(os.path.join(folder, "*.bam"))
    alleles_path = find_file(os.path.join(folder, "tables", "*-allAlleles.txt"))

    output_dir = os.path.join(folder, "MinVar_results")
    os.makedirs(output_dir, exist_ok=True)

    df = pd.read_csv(alleles_path, sep="\t")

    df = df[df.Frequency != 1]

    features = None
    if gff_path:
        regions = load_gff_regions(gff_path)
        features = parse_gff(gff_path)
        df = df[df.Position.isin(regions)]

    df_min = df[
        (df.Frequency > min_freq) &
        (df.Frequency < 0.5) &
        (df.Count > 20)
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

    results = []

    for _, row in df_min.iterrows():

        pos = row.Position
        alt = row.Allele

        cons_row = df_cons[df_cons.Position == pos]
        if cons_row.empty:
            continue

        cons = cons_row.iloc[0].Allele

        alt_metrics = analyze_variant_position(bam, pos, alt)
        cons_metrics = analyze_variant_position(bam, pos, cons)

        if not alt_metrics or not cons_metrics:
            continue

        status = "PASS"

        if alt_metrics["extreme_fraction"] > 0.5:
            if abs(alt_metrics["extreme_fraction"] - cons_metrics["extreme_fraction"]) > 0.2:
                status = "FAIL_position_bias"

        if alt_metrics["strand_bias"] > 0.8:
            status = "FAIL_strand_bias"

        if alt_metrics["mean_base_quality"] < 25:
            status = "FAIL_low_quality"
        
        # -------- annotation --------
        if features:
            annotation = annotate_position(pos, features)
        else:
            annotation = {
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
            "Strand_bias": alt_metrics["strand_bias"],
            "Mean_base_quality": alt_metrics["mean_base_quality"],
            "Status": status,

            # GFF FEATURES
            "Feature_type": annotation["Feature_type"],
            "Feature_ID": annotation["Feature_ID"],
            "Gene_name": annotation["Gene_name"],
            "Product": annotation["Product"],
            "Region_start": annotation["Region_start"],
            "Region_end": annotation["Region_end"]
        })


    bam.close()

    df_out = pd.DataFrame(results)

    df_out.to_csv(os.path.join(output_dir, "minority_variants.tsv"), sep="\t", index=False)

    write_metrics_description(output_dir)

    return df_out


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
Fraction of reads supporting the variant located at the first or last 10% of reads.
High values (>0.5) suggest sequencing artifacts.

Consensus_extreme_fraction:
Same metric for the consensus allele. Used as internal control.

Strand_bias:
Difference between forward and reverse reads supporting the variant.
Values close to 1 indicate strong bias (artifact).

Mean_base_quality:
Average Phred quality of bases supporting the variant.
Low values (<25) indicate unreliable calls.

Status:
PASS → likely real variant
FAIL_position_bias → positional artifact
FAIL_strand_bias → strand artifact
FAIL_low_quality → poor quality support
"""

    with open(os.path.join(outdir, "README_metrics.txt"), "w") as f:
        f.write(text)