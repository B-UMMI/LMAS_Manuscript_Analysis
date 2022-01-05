# LMAS Manuscript Analysis

                      _    __  __   _   ___
       /\︵︵/\      | |  |  \/  | /_\ / __|
      (◕('人')◕)     | |__| |\/| |/ _ \\__ \
         |︶|        |____|_|  |_/_/ \_\___/

         Last Metagenomic Assembler Standing

## Overview

To evaluate the performance of the assemblers in [LMAS](https://github.com/B-UMMI/LMAS), the eight bacterial genomes and four plasmids of the ZymoBIOMICS Microbial Community Standards](https://zenodo.org/record/4588970#.YEeA83X7RhE) as reference. 
It contains complete sequences for the following species:

- *Bacillus subtilis* 
- *Enterococcus faecalis*
- *Escherichia coli*
   - *Escherichia coli* plasmid
- *Lactobacillus fermentum*
- *Listeria monocytogenes*
- *Pseudomonas aeruginosa*
- *Salmonella enterica*
- *Staphylococcus aureus*
   - *Staphylococcus aureus* plasmid 1
   - *Staphylococcus aureus* plasmid 2
   - *Staphylococcus aureus* plasmid 3


As input we used the raw sequence reads of mock communities with an even and logarithmic distribution of species, from real sequencing runs and simulated read datasets, with and without error, matching the distribution of species in each sample. 
Our dataset is composed of samples:

- ENN: in silico generated evenly distributed without error
- EMS: (in silico generated evenly distributed with Illumina MiSeq error model
- ERR2984773: evenly distributed real Illumina MiSeq sample
- LNN: in silico generated logarithmically distributed without error
- LHS: in silico generated logarithmically distributed with Illumina HiSeq error model
- ERR2935805: logarithmically distributed real Illumina HiSeq sample

## Data availability

The datasets analysed during the current study are available in the Zenodo repository, under https://doi.org/10.5281/zenodo.4588969. Real sequencing data of the ZymoBIOMICS Microbial Community Standards is available under accessions [ERR2984773](https://www.ebi.ac.uk/ena/browser/view/ERR2984773?show=reads) and [ERR2935805](https://www.ebi.ac.uk/ena/browser/view/ERR2935805?show=reads).

## Assessment of assembly success

The complete set of results for 3 LMAS runs for the raw sequence reads of mock communities with an even and logarithmic distribution of species, from real sequencing runs and simulated read datasets, with and without error, matching the intended distribution of species in each sample for the eight bacterial genomes and four plasmids of the ZymoBIOMICS Microbial Community Standards as reference is available in the [Results folder](https://github.com/B-UMMI/LMAS_Manuscript_Analysis/tree/main/Results).

For the assessment of the assembly success for each sample, availabe in the [Analysis folder](https://github.com/B-UMMI/LMAS_Manuscript_Analysis/tree/main/Analysis), the different metrics for all LMAS runs were combined and descriptive statistics, such as the average value, standard deviation, minimum and maximum, were obtained through Python’s Pandas `describe` function (https://pandas.pydata.org/, https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.describe.html). 

For descriptive statistics on several assembler by each assembly type (genomic or metagenomic) and each assembler algorithm (single or multiple k-mer), the use of median was preferred due to its higher robustness against outliers and the high range of the distribution of the results. Plotly (https://plotly.com/python/) was used to compute the graphs aggregating the results obtained.

The top result of each assembler for each sample was selected, based on the following criteria:

- For the number of uncalled bases, number of misassembled contigs and number of misassembly events, the lower the value, the better, with the exception of 0 for the number of contigs;
- For the percentage of mapped reads and N50, the higher the value, the better; 
- The number of basepairs, the best results was the one closest to the target value of the number of basepairs in the reference replicons.
- For reference-specific metrics, in addition to the ones stated above when applicable (Number of contigs produced, number of uncalled bases, number of misassembled contigs and number of misassembly events), the following criteria were used:
- For the L90 metric, the lower value was better, with the exception of 0;
- For LSA, NA, NG, breadth of coverage, identity and lowest identity, the higher the value the better;
- For multiplicity, parsimony and validity, the closer to 1, the better. 

To obtain the worst value in each metric, the opposite criteria were used. The normalised score for each metric was obtained from the best result for each assembler in each sample through the equation below. For the assessment of assembler consistency, each contig for each assembler was considered the same as its size was exactly the same in each LMAS run. 

![equation](https://imgur.com/a/N4AvpJx)

## Taxonomic composition

The taxonomic composition of the ZymoBIOMICS standard samples, both real and mocks was determined through [Kraken2]() using the Standard Database (https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20210517.tar.gz). The following command was used:  

    kraken2 --output $sample.kraken --report $sample.kraken_report --memory-mapping --paired --gzip-compressed  $fastq_pair[0] $fastq_pair[1]

where `$sample` is the sample name, `$fastq_pair[0]` contains the forward-facing reads, and `$fastq_pair[1]` the reverse-facing reads. The results are availabe in the [Kraken folder](https://github.com/B-UMMI/LMAS_Manuscript_Analysis/tree/main/Kraken).

The processing of the kraken reports was performed through custom python code (https://github.com/B-UMMI/LMAS_Manuscript_Analysis/blob/main/Analysis/Sample%20taxonomy.ipynb) where all the percentage of reads that matched for the species in the dataset were saved, as well as the percentage of unclassified reads. For the *Lactobacillus fermentum*, as in the Standard Kraken database no general Species level classification is available, the percentage of reads was calculated as the sum of all reads aligning to one of the *L. fermentum* subspecies. The rest of the reads that were classified as any other species were saved conjunctively as “Other”. Supplemental Table S20 contains the percentage of classified reads for each of the species in the community, as well as “other” and unclassified reads. 

## Citation and Contacts

LMAS, and the subsequent data analysis, is developed at the Molecular [Microbiology and Infection Unit (UMMI)](http://darwin.phyloviz.net/wiki/doku.php) at the [Instituto de Medicina Molecular Joao Antunes](https://imm.medicina.ulisboa.pt/en/), in collaboration with [Microbiology, Advanced Genomics and Infection Control Applications Laboratory (MAGICAL)](https://morangiladlab.com) at the [Faculty of Health Sciences, Ben-Gurion University of the Negev](https://in.bgu.ac.il/en/fohs/Pages/default.aspx). 

If you use LMAS please [cite LMAS repository](https://github.com/cimendes/LMAS/blob/main/CITATION.cff).