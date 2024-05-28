# The evolutionary fate of Neanderthal DNA in 30,780 admixed genomes with recent African-like ancestry

This repository contains code used in XX. For a detailed description of methods, see XX. All analyses were implemented in a Snakemake workflows<sup>1</sup> (note that the workflows contained in `snakefiles` do not work independently but only in conjunction). 

#### Workflow overview
We first build reference panels of African, European, and Native American/East Asian genomes. Here, we use 504 African genomes from the 1000 genomes project (1KGP-AFR)<sup>2</sup> + 563 African genomes from All of Us database v7.1 (AOU-AFR)<sup>3,4</sup>, 503 1KGP-EUR genomes + 10,000 AOU-EUR genomes, and 504 1KGP-EAS genomes + 71 AOU-EAS/Native American genomes. We then identify 30,780 recently admixed individuals in the target VCF file with predominantly African-like and European-like ancestry (i.e., "AFR-like" >= 50% & "EUR-like" >= 10% & "AFR- + EUR-like" >= 95%) in All of Us, using previously estimated ancestry proportions<sup>3,4,6</sup>. We then call Neanderthal introgressed segments using [IBDmix v1.0.1](https://github.com/PrincetonUniversity/IBDmix)<sup>6</sup>, and compare observed amounts of Neanderthal ancestry to expected amounts in the admixed individuals. The expected amount of Neanderthal ancestry per admixed genome is calculated based on the mean introgression amounts per individual and the recent ancestry proportions. After phasing genotypes calls with Beagle v5.4<sup>7</sup> and painting local ancestry using [FLARE v0.5.1](https://github.com/browning-lab/flare)<sup>8</sup>, we use similar approach to identiy putatively selected Neanderthal introgressed segments in the admixed individuals, i.e., are regions that are depleted/enriched of Neanderthal ancestry, as well as novel introgression desert-like regions. We control for for false-positive finding using coalescence simulations under a plausible demographic model in msprime <sup>9</sup>. We finally annotate regions of interest with various population genetic statistics (e.g., GERP scores)<sup>10</sup> and functional annotations (e.g., ClinVar)<sup>11</sup> using vcfanno v0.3.3<sup>12</sup>, among others

#### Software dependencies
The Snakemake environment can be set up using the following command as it also contains some additional software requirements:

`conda env create -f snakemake.yaml`

Most other required software will be installed automatically by Snakemake via the provided enviroment files. However, there are a few tools that need to be installed/download manually:
- [IBDmix v1.0.1](https://github.com/PrincetonUniversity/IBDmix)<sup>5</sup>
- [BCFtools v1.14-36-g95760eb](http://www.htslib.org/download/)<sup>14</sup>
- [BEDTools v2.30.0](https://bedtools.readthedocs.io/en/latest/index.html)<sup>15</sup>
- [bwa v0.7.17-r1188](https://github.com/lh3/bwa)<sup>16</sup>
- [FLARE v0.5.1](https://github.com/browning-lab/flare)<sup>5</sup>
- [ENSEMBL API](https://useast.ensembl.org/info/docs/api/api_git.html)<sup>17</sup>
- [Rye v0.1](https://github.com/healthdisparities/rye)<sup>2</sup>

The paths to these tools can be updated in the `config/config.yaml` file.

#### Running the workflow

**Disclaimer: The workflow has not been tested on other datasets. If you have trouble executing it, feel free to reach out and we're happy to assist.**

If analyzing your own data, set **dataset** field in `config/config.yaml` to an arbitrary name (not *test* or *AOU*) and VCFs with phased genotype calls for the target individuals must be provided. The VCFs should be named following this convention: `data_path + target_individuals_chr{chr}.vcf.gz`, where the `data_path` can be updated in the `config/config.yaml` file. A subset of these individuals will be selected based on specified global ancestry proportions (i.e., "AFR-like" >= 50% & "EUR-like" >= 10% & "AFR- + EUR-like" >= 95%). To this end, we build continental reference panels from 1KGP<sup>2</sup> and HGDP data<sup>13</sup> and estimate ancestry proportions in supervised manner using [Rye v0.1](https://github.com/healthdisparities/rye)<sup>5</sup>. Additionally, it is assumed that a file with related individuals that are to be excluded from the analyses is provided as well as a personal LDlinktoken. All other data files will be automatically retrieved or are contained in this repository. The urls to such files can be updated in the `config/config.yaml`. The acronyms of included reference population as well as their mapping to continental super populations should also also be updated in`config/config.yaml`.

The workflow can then be locally executed, using:<br>
```snakemake --use-conda --cores all```

**Also note that the current workflow is tailored to three-way admixed individuals with predominantly recent African- and European-like ancestry and some filter criteria and modeling steps may need to be updated to work for other admixture cases.**

#### Precomputed result files
Due to restriction from All of Us, we cannot make available the callset of introgressed segments on an individual level. However, we have uploaded files with calculated introgression and local ancestry frequencies in 50 kb windows. Note that all coordinates are relative to hg38:

- `results/ibdmix_Vindija33.19/ibdmix_results_masked_denisovan_combined_50kb_4.0LOD_coverage_per_window50000_s_10000.bed`<br> The contains the number of introgressed base pairs per 50 kb across all individuals for all superpopulations (AFR, EUR, EAS(/NA), and AMR). The 4th column indicates the superpopulation and the 5th the number of introgressed base pairs across all individuals. To obtain frequencies divide the number of introgressed base pairs by the window size and the number of individuals in the respective superpopulation (i.e., AFR=1067, EUR=10503, EAS=567, and AMR=30780).

- `results/ibdmix_Vindija33.19/ibdmix_results_masked_denisovan_combined_50kb_4.0LOD_afr_masked_coverage_per_window50000_s_10000.bed`<br> Same as above but this time the number of introgressed base pairs is calculated after applying the African mask, i.e., removing any introgressed segment that overlaps with a segment in an African reference genome.

- `results/ibdmix_Vindija33.19/ibdmix_results_masked_denisovan_combined_50kb_4.0LOD_afr_masked_coverage_per_window50000_s_10000_expectations.bed`<br> Contains the number of introgressed base pairs per 50 kb window across all individuals for all superpopulation normalized by the window size (to obtain frequencies you still need to divide by the number of individuals) and the number of local ancestry haplotypes overlapping a given window in the 30,780 admixed individuals. It also contains the expected number introgressed haplotypes in the admixed individuals as well as the Agresti-Coull correct observed number of archaic haplotypes in the admixed individuals (`freq_updated`).

- `results/ibdmix_Vindija33.19/ibdmix_results_masked_denisovan_combined_50kb_4.0LOD_afr_masked_coverage_per_window50000_s_10000_pvalues.bed`<br> For a subset of the windows in the above file (i.e., is windows with expectations >= 1, have less than 50% masked sites and that pass local ancestry and recombination rate filters), we calculated probabilities for observing the number of observed archaic haplotypes given the expectations. The field `corrected` gives the Bonferroni corrected p-values. The fields `lower_drift` and `upper_drift` provide the confidence intervals for what sort of introgression frequencies (i.e., already divided by the sample size of 30780) you could expect to see after 15 generations of drift. 

- `results/ibdmix_Vindija33.19/<superpop>_introgression_frequencies_and_rank_callable_windows.bed`<br> Contains information about the number of masked sites and introgressed base pairs for large windows of various sizes (8-15 Mb) that had less than 50% masked sites as well as their percentile rank in terms of number of introgressed base pairs. 

- `results/ibdmix_Vindija33.19/<superpop>_introgression_frequencies_and_rank_callable_windows_afr_masked.bed`<br> Same as above but number of introgressed base pairs are calculated based on the African masked call set of Neanderthal introgressed segments. 

- `results/ibdmix_Vindija33.19/AMR_novel_introgression_deserts.bed`<br> Genomic coordinates of novel introgression desert-like regions in the admixed genomes.

- `results/ibdmix_Vindija33.19/AMR_novel_introgression_deserts_pvalues.bed`<br> Genomic coordinates of novel introgression desert-like regions in the admixed genomes. Contains the number of introgressed base pairs and the effective number of local ancestry haplotypes for each large window that is in the bottom 5% in terms of introgressed based pairs. It also contains the expected number introgressed haplotypes in the admixed individuals as well as the Agresti-Coull correct observed number of archaic haplotypes in the admixed individuals (`freq_updated`). We calculated probabilities for observing the number of observed archaic haplotypes given the expectations. The field `corrected` gives the Bonferroni corrected p-values. We also provide the percentile of the window in EUR and EAS in terms of number of introgressed base pairs.

If you would like to access individual level calls of introgressed segments or local ancestry, feel free to reach out and we are happy to share the files within the All of Us research bench.  

#### References
1. Mölder, F., Jablonski, K.P., Letcher, B., et al. (2021) Sustainable data analysis with Snakemake. F1000Res 10, 33. [https://doi.org/10.12688/f1000research.29032.1](https://doi.org/10.12688/f1000research.29032.1)
2. 1000 Genomes Project Consortium, Auton, A., Brooks, L. D., et al. (2015). A global reference for human genetic variation. Nature, 526(7571), 68–74. [https://doi.org/10.1038/nature15393](https://doi.org/10.1038/nature15393)
3.  All of Us Research Program Investigators , Denny J, Rutter J, et al (2019) The
“all of us” research program. New England Journal of Medicine 381(7):668–676. [https://doi.org/10.1056/NEJMsr1809937](https://doi.org/10.1056/NEJMsr1809937)
4. Bick AG, Metcalf GA, Mayo KR, et al (2024) Genomic data in the All of Us Research Program.704
Nature pp 1–7. [https://doi.org/10.1038/s41586-023-06957-x](https://doi.org/10.1038/s41586-023-06957-x)
5. Andrew B Conley, Lavanya Rishishwar, Maria Ahmad, et al. (2023) Rye: genetic ancestry inference at biobank scale, Nucleic Acids Research, Volume 51, Issue 8, 8 May 2023, Page e44, [https://doi.org/10.1093/nar/gkad149](https://doi.org/10.1093/nar/gkad149)
6. Chen, L., Wolf, A. B., Fu, W., et al. (2020) Identifying and Interpreting Apparent Neanderthal Ancestry in African Individuals. Cell, 180(4), 677–687.e16. [https://doi.org/10.1016/j.cell.2020.01.012](https://doi.org/10.1016/j.cell.2020.01.012)
7. B L Browning, X Tian, Y Zhou, and S R Browning (2021) Fast two-stage phasing of large-scale sequence data. Am J Hum Genet 108(10):1880-1890. [https://doi:10.1016/j.ajhg.2021.08.005](https://doi:10.1016/j.ajhg.2021.08.005)
8. Browning, S. R., Waples, R. K., & Browning, B. L. (2023) Fast, accurate local ancestry inference with FLARE. American journal of human genetics, 110(2), 326–335. [https://doi.org/10.1016/j.ajhg.2022.12.010](https://doi.org/10.1016/j.ajhg.2022.12.010)
9. Franz Baumdicker, Gertjan Bisschop, Daniel Goldstein, et al (2022) Efficient ancestry and mutation simulation with msprime 1.0, Genetics, Volume 220, Issue 3, iyab229, [https://doi.org/10.1093/genetics/iyab229](https://doi.org/10.1093/genetics/iyab229)
10. Huber CD, Kim BY, Lohmueller KE (2020) Population genetic models of GERP scores suggest pervasive turnover of constrained sites across mammalian evolution. PLOS Genetics 16(5): e1008827. [https://doi.org/10.1371/journal.pgen.1008827](https://doi.org/10.1371/journal.pgen.1008827)
11. Melissa J Landrum, Jennifer M Lee, Mark Benson, et al. (2018) ClinVar: improving access to variant interpretations and supporting evidence, Nucleic Acids Research, Volume 46, Issue D1, 4 January 2018, Pages D1062–D1067, [https://doi.org/10.1093/nar/gkx1153](https://doi.org/10.1093/nar/gkx1153)
12. Pedersen, B.S., Layer, R.M. & Quinlan, A.R. (2016) Vcfanno: fast, flexible annotation of genetic variants. Genome Biol 17, 118. [https://doi.org/10.1186/s13059-016-0973-5](https://doi.org/10.1186/s13059-016-0973-5)
13. Anders Bergström et al. (2020) Insights into human genetic variation and population history from 929 diverse genomes.Science367,eaay5012. [10.1126/science.aay5012](10.1126/science.aay5012)
14. Petr Danecek, James K Bonfield, Jennifer Liddle, et al. (2021), Twelve years of SAMtools and BCFtools, GigaScience, Volume 10, Issue 2, February 2021, giab008, [https://doi.org/10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008)
15. Quinlan, A. R., & Hall, I. M. (2010). BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics (Oxford, England), 26(6), 841–842. [https://doi.org/10.1093/bioinformatics/btq033](https://doi.org/10.1093/bioinformatics/btq033)
16. Heng Li, Richard Durbin, (2009) Fast and accurate short read alignment with Burrows–Wheeler transform, Bioinformatics, Volume 25, Issue 14, July 2009, Pages 1754–1760, [https://doi.org/10.1093/bioinformatics/btp324](https://doi.org/10.1093/bioinformatics/btp324)
17. Andrew Yates, Kathryn Beal, Stephen Keenan, et al. (2015) The Ensembl REST API: Ensembl Data for Any Language, Bioinformatics, Volume 31, Issue 1, January 2015, Pages 143–145, [https://doi.org/10.1093/bioinformatics/btu613](https://doi.org/10.1093/bioinformatics/btu613)

