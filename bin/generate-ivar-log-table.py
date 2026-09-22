#!/usr/bin/env python3

import argparse
import csv
import re
from pathlib import Path


def parse_ivar_log(log_file):
    """Parse primer counts and summary statistics from an iVar log."""

    primer_counts = {}
    summary = {}

    # Primer lines, e.g.
    # Uni12-PB1       0
    # PB1-683         5850
    primer_pattern = re.compile(r"^(\S+)\s+(\d+)\s*$")

    # Summary patterns
    trimmed_pattern = re.compile(
        r"Trimmed primers from ([\d.]+)% \((\d+)\) of reads"
    )

    below_min_length_pattern = re.compile(
        r"([\d.]+)% \((\d+)\) of reads were quality trimmed "
        r"below the minimum length"
    )

    outside_primer_pattern = re.compile(
        r"([\d.]+)% \((\d+)\) of reads started outside of primer regions"
    )

    unmapped_pattern = re.compile(
        r"(\d+) unmapped reads"
    )

    insert_smaller_pattern = re.compile(
        r"([\d.]+)% \((\d+)\) of reads had their insert size "
        r"smaller than their read length"
    )

    with log_file.open() as f:
        for line in f:
            line = line.strip()

            # Primer counts
            match = primer_pattern.match(line)
            if match:
                primer, count = match.groups()
                primer_counts[primer] = int(count)
                continue

            # Percentage/count of reads with primers trimmed
            match = trimmed_pattern.search(line)
            if match:
                summary["percentage_trimmed"] = float(match.group(1))
                summary["count_trimmed"] = int(match.group(2))
                continue

            # Percentage/count quality trimmed below minimum length
            match = below_min_length_pattern.search(line)
            if match:
                summary["percentage_below_min_length"] = float(match.group(1))
                summary["count_below_min_length"] = int(match.group(2))
                continue

            # Percentage/count of reads starting outside primer regions
            match = outside_primer_pattern.search(line)
            if match:
                summary["percentage_outside_primer"] = float(match.group(1))
                summary["count_outside_primer"] = int(match.group(2))
                continue

            # Count of unmapped reads
            match = unmapped_pattern.search(line)
            if match:
                summary["count_unmapped"] = int(match.group(1))
                continue

            # Percentage/count of reads with insert smaller than read length
            match = insert_smaller_pattern.search(line)
            if match:
                summary["percentage_insert_smaller_than_read"] = float(
                    match.group(1)
                )
                summary["count_insert_smaller_than_read"] = int(
                    match.group(2)
                )
                continue

    return primer_counts, summary


def main():
    parser = argparse.ArgumentParser(
        description="Consolidate iVar primer trimming logs into a TSV."
    )

    parser.add_argument(
        "--directory",
        required=True,
        type=Path,
        help="Directory containing iVar log files.",
    )

    parser.add_argument(
        "--extension",
        default="ivar.log",
        help="Extension/suffix of iVar log files (default: ivar.log).",
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output TSV file.",
    )

    args = parser.parse_args()

    # Find log files
    log_files = sorted(
        f
        for f in args.directory.iterdir()
        if f.is_file() and f.name.endswith(args.extension)
    )

    if not log_files:
        parser.error(
            f"No files ending with '{args.extension}' found in "
            f"{args.directory}"
        )

    sample_counts = {}
    sample_summary = {}
    primers = []

    for log_file in log_files:
        # Example:
        # sample001.ivar.log -> sample001
        sample_name = log_file.name[: -len(args.extension)].rstrip(".")

        counts, summary = parse_ivar_log(log_file)

        sample_counts[sample_name] = counts
        sample_summary[sample_name] = summary

        # Keep primers in the order they are first encountered
        for primer in counts:
            if primer not in primers:
                primers.append(primer)

    # Summary metrics, in desired output order
    summary_metrics = [
        "percentage_trimmed",
        "count_trimmed",
        "percentage_below_min_length",
        "count_below_min_length",
        "percentage_outside_primer",
        "count_outside_primer",
        "count_unmapped",
        "percentage_insert_smaller_than_read",
        "count_insert_smaller_than_read",
    ]

    samples = list(sample_counts)

    # Write output
    with args.output.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")

        # Header
        writer.writerow(["Metric"] + samples)

        # Primer rows
        for primer in primers:
            writer.writerow(
                [primer]
                + [
                    sample_counts[sample].get(primer, 0)
                    for sample in samples
                ]
            )

        # Summary rows
        for metric in summary_metrics:
            writer.writerow(
                [metric]
                + [
                    sample_summary[sample].get(metric, "")
                    for sample in samples
                ]
            )

    print(
        f"Wrote {len(primers)} primers + "
        f"{len(summary_metrics)} summary metrics × "
        f"{len(samples)} samples"
    )
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
