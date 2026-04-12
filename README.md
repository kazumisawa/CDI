# CDI
Coalescent Distortion Index (CDI) computed from VCF genotype data
# CDI: Coalescent Distortion Index

## Definition

For a diploid individual with two homologous chromosomes \(A\) and \(B\), and a reference sequence \(R\), we define the **Coalescent Distortion Index (CDI)** as

\[
\mathrm{CDI}
=
\frac{d(R,A) + d(R,B)}{2}
-
d(A,B),
\]

where \(d(X,Y)\) denotes the genetic distance between sequences \(X\) and \(Y\).

In words:

> **CDI is the mean distance between the reference sequence and the two homologous chromosomes of an individual, minus the distance between the two homologous chromosomes themselves.**

---

## Gene genealogy interpretation

We assume a genomic region with **no recombination**, so that the three sequences \(R\), \(A\), and \(B\) share a single local gene genealogy:
    R
    |
    |
    O
   / \
  A   B
  Here:

- \(O\) is the most recent common ancestor of \(R\), \(A\), and \(B\),
- branch lengths represent accumulated mutations,
- genetic distances are additive along the tree.

Under this additive tree metric, we have

\[
\begin{aligned}
d(R,A) &= d(O,R) + d(O,A), \\
d(R,B) &= d(O,R) + d(O,B), \\
d(A,B) &= d(O,A) + d(O,B).
\end{aligned}
\]

Substituting these relations into the definition of CDI yields

\[
\mathrm{CDI}
=
d(O,R)
-
\frac{d(O,A) + d(O,B)}{2}.
\]

Thus, CDI directly measures how much the reference branch is distorted relative to the *coalescent structure* defined by the two homologous chromosomes of an individual.

---

## Interpretation of CDI = 0

From the genealogical expression above,

\[
\mathrm{CDI} = 0
\quad \Longleftrightarrow \quad
d(O,R) = \frac{d(O,A) + d(O,B)}{2}.
\]

This condition has a clear interpretation:

> **CDI = 0 means that the reference sequence lies at the average coalescent (genealogical) position of the two homologous chromosomes.**

In other words, the reference is neither systematically closer to the ancestral node \(O\) nor systematically more diverged than the individual’s two alleles.

Under a neutral Wright–Fisher model with a representative reference sequence, this situation is expected on average.

---

## Site-wise formulation (biallelic, unphased)

In practice, genetic distances are computed site by site from VCF genotype data.  
We assume:

- biallelic sites only,
- diploid individuals,
- **unphased** genotypes,
- the VCF `REF` allele represents the reference sequence \(R\).

For a single site, the contribution to CDI is:

| Genotype | Contribution to CDI |
|----------|---------------------|
| 0/0 (REF/REF) | 0 |
| 0/1 (REF/ALT) | \(-\dfrac{1}{2}\) |
| 1/1 (ALT/ALT) | \(+1\) |

Summing over all sites in a region, CDI can be written as

\[
\boxed{
\mathrm{CDI}
=
\#(1/1)
-
\frac{1}{2}\,\#(0/1)
}
\]

where:

- \(\#(1/1)\) is the **number of sites** at which the individual is homozygous for the alternative allele,
- \(\#(0/1)\) is the **number of heterozygous sites**,
- \(\#(\cdot)\) denotes a **count of sites**, not a frequency.

Sites with genotype 0/0 do not contribute and are therefore omitted from the expression.

---

## Conceptual summary

- \(d(A,B)\) reflects **within-individual genetic diversity** (heterozygosity).
- \(\frac{d(R,A) + d(R,B)}{2}\) reflects **divergence between an individual and the reference**.
- CDI contrasts these two quantities in terms of coalescent branch lengths.

Therefore:

> **The Coalescent Distortion Index quantifies how much the reference sequence deviates from the average coalescent position defined by an individual’s two homologous chromosomes, after accounting for within-individual heterozygosity.**

This interpretation is consistent across
gene genealogies,
additive distance models,
and site-wise genotype counts.
