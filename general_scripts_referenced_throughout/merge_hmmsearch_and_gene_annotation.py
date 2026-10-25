#!/usr/bin/env python3
"""
merge_gff3_hmmsearch.py
-----------------------
Merges a Bakta GFF3 file with an hmmsearch tblout file on locus_tag.

Usage:
    python merge_gff3_hmmsearch.py \
        --gff your_file.gff3 \
        --tblout your_file.tblout \
        --out merged_output.tsv \
        [--join inner|left|right|outer]

Join types:
    inner  – only locus_tags present in BOTH files (default)
    left   – all GFF3 rows, hmmsearch columns filled with NA where no match
    right  – all hmmsearch rows, GFF3 columns filled with NA where no match
    outer  – all rows from both files
"""

import argparse
import csv
import sys
import pandas as pd


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def parse_gff3(path: str) -> pd.DataFrame:
    """
    Parse a GFF3 file, skipping lines that start with '#'.
    Extracts the 9 standard GFF3 columns plus selected attributes
    (ID, Name, locus_tag, product) from column 9.
    Returns a DataFrame with one row per feature.
    """
    GFF3_COLS = [
        "seqid", "source", "type", "start", "end",
        "score", "strand", "phase", "attributes_raw"
    ]

    rows = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 9:
                continue          # skip malformed lines
            rows.append(parts)

    df = pd.DataFrame(rows, columns=GFF3_COLS)

    # Parse key=value pairs from the attributes column
    def parse_attrs(attr_str):
        attrs = {}
        for kv in attr_str.split(";"):
            kv = kv.strip()
            if "=" in kv:
                k, v = kv.split("=", 1)
                attrs[k] = v
        return attrs

    attr_dicts = df["attributes_raw"].apply(parse_attrs)

    for key in ("ID", "Name", "locus_tag", "product"):
        df[f"gff3_{key}"] = attr_dicts.apply(lambda d: d.get(key, pd.NA))

    # Drop the raw attributes column (already unpacked above)
    df = df.drop(columns=["attributes_raw"])

    return df


def parse_tblout(path: str) -> pd.DataFrame:
    """
    Parse an hmmsearch --tblout file.
    Lines starting with '#' are comments/headers and are skipped.
    Columns are whitespace-delimited; the description field (col 19+)
    may contain spaces so we split into at most 19 parts.
    """
    TBLOUT_COLS = [
        "target_name",       # col  1  – locus tag of the hit
        "target_accession",  # col  2
        "query_name",        # col  3  – HMM profile name
        "query_accession",   # col  4
        "full_evalue",       # col  5
        "full_score",        # col  6
        "full_bias",         # col  7
        "dom_evalue",        # col  8
        "dom_score",         # col  9
        "dom_bias",          # col 10
        "exp",               # col 11
        "reg",               # col 12
        "clu",               # col 13
        "ov",                # col 14
        "env",               # col 15
        "dom",               # col 16
        "rep",               # col 17
        "inc",               # col 18
        "description",       # col 19+ (may contain spaces)
    ]

    rows = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            # Split into max 19 parts so description stays together
            parts = line.split(None, 18)
            if len(parts) < 18:
                continue      # malformed row
            if len(parts) == 18:
                parts.append("")   # no description
            rows.append(parts)

    df = pd.DataFrame(rows, columns=TBLOUT_COLS)

    # Convert numeric columns
    for col in ("full_evalue", "full_score", "full_bias",
                "dom_evalue", "dom_score", "dom_bias", "exp"):
        df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


# ---------------------------------------------------------------------------
# Merge
# ---------------------------------------------------------------------------

def merge(gff3_df: pd.DataFrame, tblout_df: pd.DataFrame,
          join: str) -> pd.DataFrame:
    """
    Merge on locus_tag (GFF3) == target_name (tblout).
    """
    merged = gff3_df.merge(
        tblout_df,
        left_on="gff3_locus_tag",
        right_on="target_name",
        how=join,
    )
    # Put the locus_tag/target_name first for readability
    cols = merged.columns.tolist()
    priority = ["gff3_locus_tag", "target_name", "query_name",
                "query_accession", "full_evalue", "full_score",
                "gff3_Name", "gff3_product",
                "seqid", "type", "start", "end", "strand"]
    front = [c for c in priority if c in cols]
    rest  = [c for c in cols if c not in front]
    return merged[front + rest]


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Merge Bakta GFF3 with hmmsearch tblout on locus_tag."
    )
    parser.add_argument("--gff",   required=True, help="Path to .gff3 file")
    parser.add_argument("--tblout", required=True, help="Path to .tblout file")
    parser.add_argument("--out",    required=True, help="Output TSV path")
    parser.add_argument(
        "--join",
        default="inner",
        choices=["inner", "left", "right", "outer"],
        help="Join type (default: inner)"
    )
    args = parser.parse_args()

    print(f"[1/4] Parsing GFF3:    {args.gff}")
    gff3_df = parse_gff3(args.gff)
    print(f"      → {len(gff3_df)} features loaded")

    print(f"[2/4] Parsing tblout:  {args.tblout}")
    tblout_df = parse_tblout(args.tblout)
    print(f"      → {len(tblout_df)} hits loaded")

    print(f"[3/4] Merging ({args.join} join) on locus_tag / target_name …")
    merged = merge(gff3_df, tblout_df, args.join)
    print(f"      → {len(merged)} rows in merged table")

    print(f"[4/4] Writing output:  {args.out}")
    merged.to_csv(args.out, sep="\t", index=False)
    print("Done.")

    # Quick summary
    if "query_name" in merged.columns:
        n_matched = merged["target_name"].notna().sum()
        print(f"\nSummary: {n_matched} GFF3 features matched an hmmsearch hit.")
        print(f"Unique HMM profiles matched: "
              f"{merged['query_name'].nunique()}")
    
if __name__ == "__main__":
    main()
