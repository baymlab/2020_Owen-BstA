# 2020_Owen-BstA

## About
Code corresponing to the paper [Prophage-encoded abortive infection proteins that encode cognate self-immunity]()

## Data
- BstA homologs
  - [Sequences](https://github.com/baymlab/2020_Owen-BstA/blob/master/data/BstAhomologs/bsta_homologs.faa)
  - [Homologs metadata](https://github.com/baymlab/2020_Owen-BstA/blob/master/data/BstAhomologs/bsta_homologs.tsv) <!--TO_DO-->
  - [BstA profile-hmm](https://github.com/baymlab/2020_Owen-BstA/blob/master/data/BstAhomologs/bsta.hmm)
- [Raw annotation results]() <!--TO_DO: Zenodo?-->
- Processed results
  - [Figure 2.A: Annotation counts](https://github.com/baymlab/2020_Owen-BstA/blob/master/data/Figures/2A/annotation_counts.tsv) <!--TO_DO: annot tables and counts-->
  - [Figure 2.B: Annotation of genomic neighborhoods](https://github.com/baymlab/2020_Owen-BstA/tree/master/data/Figures/2B)
  - [Figure 2.C: BstA homologs alignment](https://github.com/baymlab/2020_Owen-BstA/blob/master/data/Figures/2C/bsta_alignment.afa)

## Code
- [`Snakefile`](https://github.com/baymlab/2020_Owen-BstA/blob/master/Snakefile): A Snakemake pipeline which does the core genomic neighborhood analysis. It consists of downloading a specified sequence region from NCBI, annotating the region with `Prokka`, running `hmmscan`, and updating the annotations with BstA and `Pfam` hits. This produces the [raw annotation results](TO_DO). A [selection](https://github.com/baymlab/2020_Owen-BstA/blob/master/data/Figures/2B/selected_genomic_neighborhoods.tsv) of these results was processed to create *Figure 2.B*.
2. [`ipynb_notebook/`](https://github.com/baymlab/2020_Owen-BstA/tree/601ce14f9d81d701d49e474615e261b4d5f28230/ipynb_notebook): Contains how the pipeline was ran as well as additional processing of the raw data. The annotations analysis for *Figure 2.A* can be found in [`ipynb_notebook/load_all_tables.ipynb`](https://github.com/baymlab/2020_Owen-BstA/blob/601ce14f9d81d701d49e474615e261b4d5f28230/ipynb_notebook/load_all_tables.ipynb).
