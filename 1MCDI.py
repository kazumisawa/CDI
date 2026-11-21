#!/usr/bin/env python3

from cyvcf2 import VCF
import sys

# --------------------------------------------------
# Usage
# --------------------------------------------------
# python CDI_window.py input.vcf.gz 1000000 > output.tsv
#
# input.vcf.gz : bgzip + tabix-ed VCF
# 1000000      : window size (bp)
# --------------------------------------------------

if len(sys.argv) != 3:
    sys.stderr.write(
        f"Usage: {sys.argv[0]} input.vcf.gz window_size\n"
    )
    sys.exit(1)

vcf_file = sys.argv[1]
window_size = int(sys.argv[2])

# --------------------------------------------------
# Open VCF
# --------------------------------------------------

vcf = VCF(vcf_file)

samples = vcf.samples
n_samples = len(samples)

# --------------------------------------------------
# Output header
# --------------------------------------------------

print(
    "chrom\tstart\tend\t" +
    "\t".join(samples)
)

# --------------------------------------------------
# Process each chromosome
# --------------------------------------------------

for chrom, chrom_length in zip(vcf.seqnames, vcf.seqlens):

    if chrom_length is None:
        continue

    # ----------------------------------------------
    # Generate windows
    # ----------------------------------------------

    start = 1

    while start <= chrom_length:

        end = min(
            start + window_size - 1,
            chrom_length
        )

        region = f"{chrom}:{start}-{end}"

        # ------------------------------------------
        # initialize counters
        # ------------------------------------------

        hom_alt = [0] * n_samples
        het = [0] * n_samples

        # ------------------------------------------
        # scan variants in region
        # ------------------------------------------

        for rec in vcf(region):

            # CDI assumes biallelic sites
            if len(rec.ALT) != 1:
                continue

            for i, (a1, a2, phased) in enumerate(rec.genotypes):

                # missing genotype
                if a1 < 0 or a2 < 0:
                    continue

                # ALT/ALT  -> +1
                if a1 == 1 and a2 == 1:
                    hom_alt[i] += 1

                # HET -> -0.5
                elif a1 != a2:
                    het[i] += 1

                # REF/REF contributes 0
                else:
                    pass

        # ------------------------------------------
        # CDI calculation
        # CDI = #(1/1) - 0.5 * #(0/1)
        # ------------------------------------------

        cdi_values = []

        for i in range(n_samples):
            cdi = - hom_alt[i] + 0.5 * het[i]
            cdi_values.append(str(cdi))

        # ------------------------------------------
        # output one line per window
        # ------------------------------------------

        print(
            chrom,
            start,
            end,
            *cdi_values,
            sep="\t",
            flush=True
        )
        print(
            f"{chrom}:{start}-{end}",
            file=sys.stderr,
            flush=True
        )

        start += window_size

