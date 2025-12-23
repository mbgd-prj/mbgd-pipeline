# DomRefine

Identification of ortholog groups is a crucial step in comparative analysis of multiple genomes. Although several computational methods have been developed to create ortholog groups, most of those methods do not evaluate orthology at the sub-gene level. In our method for domain-based ortholog clustering, DomClust , proteins are split into domains on the basis of alignment boundaries identified by all-against-all pairwise comparison, but it often fails to determine appropriate boundaries.

We developed a method to improve domain-based ortholog classification using multiple alignment information. This method is based on a scoring scheme, the domain-specific sum-of-pairs (DSP) score , which evaluates ortholog clustering results at the domain level as the sum total of domain-based alignment scores. We developed a refinement pipeline to improve domain-based clustering, DomRefine , by optimizing the DSP score.

## Original version
https://mbgd.nibb.ac.jp/domrefine/

### Reference
Chiba, H. and Uchiyama, I.
Improvement of domain-level ortholog clustering by optimizing domain-specific sum-of-pairs score.
[*BMC Bioinformatics*, 15:148 (2014)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-148).
