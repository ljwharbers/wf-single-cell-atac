#!/usr/bin/env python
"""Cluster UMIs."""

from pathlib import Path

import pandas as pd
import pysam
from .util import wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("tag_bams")

    parser.add_argument(
        "--in_bam",
        help="BAM file for tagging",
        type=Path
    )

    parser.add_argument(
        "--out_bam",
        help="Path for tagged output BAM",
        type=Path
    )

    parser.add_argument(
        "--tags",
        help="Read tags TSV",
        type=Path
    )

    parser.add_argument(
        "--chrom",
        help="Chromosome name"
    )
    return parser


def add_tags(tags_file, in_bam, out_bam, chrom):
    """Add all the required tags to the BAM file."""
    df = pd.read_csv(tags_file, sep='\t', index_col=0)
    with pysam.AlignmentFile(in_bam, "rb") as bam_in:
        with pysam.AlignmentFile(out_bam, "wb", template=bam_in) as bam_out:
            for align in bam_in.fetch(contig=chrom):
                read_id = align.query_name
                try:
                    row = df.loc[read_id]
                except KeyError:
                    continue  # No barcode/umi for this read
                # uncorrectred cell barcode
                align.set_tag('CR', row['CR'], value_type="Z")
                # Correctred cell barcode
                align.set_tag('CB', row['CB'], value_type="Z")
                # barcode qscores
                align.set_tag('CY', row['CY'], value_type="Z")

                # Up to this point in the workflow the reads are in reverse
                # orientation in relation to the mRNA.
                # Flip this for the output BAMs.
                align.flag ^= 16  # reverse read alignment flag

                bam_out.write(align)


def main(args):
    """Entry point."""
    add_tags(args.tags, args.in_bam, args.out_bam, args.chrom)
