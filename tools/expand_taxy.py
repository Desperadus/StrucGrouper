#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
import pandas as pd


def sanitize_type(t: str) -> str:
    t = t.strip().lower()
    t = re.sub(r"[^\w]+", "_", t)   # non-alphanum -> underscore
    t = re.sub(r"_+", "_", t).strip("_")
    return t or "unknown"


def parse_lineage(line: str):
    if pd.isna(line) or not str(line).strip():
        return []
    parts = [p.strip() for p in str(line).split(",")]
    parsed = []
    for p in parts:
        m = re.match(r"^(.*)\s+\(([^)]+)\)$", p)
        if m:
            taxon_name = m.group(1).strip()
            taxon_type = sanitize_type(m.group(2))
            parsed.append((taxon_name, taxon_type))
        else:
            parsed.append((p, "unknown"))
    return parsed


def expand_lineage(row_lineage: str):
    parsed = parse_lineage(row_lineage)
    type_counts = {}
    for _, t in parsed:
        type_counts[t] = type_counts.get(t, 0) + 1

    type_seen = {t: 0 for t in type_counts}
    out = {}
    for name, t in parsed:
        type_seen[t] += 1
        col = f"{t}{type_seen[t]}" if type_counts[t] > 1 else t
        out[col] = name
    return out


def build_output_path(input_path: Path) -> Path:
    stem = input_path.stem  # drops .tsv (or whatever extension)
    return input_path.with_name(f"{stem}_extended.tsv")


def main():
    ap = argparse.ArgumentParser(
        description="Expand 'Taxonomic lineage' into typed columns (e.g., clade1, clade2).")
    ap.add_argument("tsv", type=Path, help="Input TSV file")
    ap.add_argument("--lineage-col", default="Taxonomic lineage",
                    help="Column name containing the lineage string (default: 'Taxonomic lineage')")
    ap.add_argument("--encoding", default="utf-8",
                    help="File encoding (default: utf-8)")
    args = ap.parse_args()

    if not args.tsv.exists():
        raise SystemExit(f"Input file not found: {args.tsv}")

    df = pd.read_csv(args.tsv, sep="\t", encoding=args.encoding)

    if args.lineage_col not in df.columns:
        raise SystemExit(
            f"Column '{args.lineage_col}' not found. Available columns: {', '.join(df.columns)}"
        )

    expanded = df[args.lineage_col].apply(expand_lineage)
    expanded_df = pd.json_normalize(expanded)

    df_out = pd.concat([df, expanded_df], axis=1)

    out_path = build_output_path(args.tsv)
    df_out.to_csv(out_path, sep="\t", index=False, encoding=args.encoding)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
