# 2020_Owen-BstA

## About
Code corresponing to the paper:

[* Prophage-encoded phage defence proteins with cognate self-immunity (Owen, *et. al.*, 2020**)](https://www.biorxiv.org/content/10.1101/2020.07.13.199331v1)

## Data
- BstA homologs
  - [Sequences](https://github.com/baymlab/2020_Owen-BstA/blob/master/data/BstAhomologs/bsta_homologs.faa)
  - [Homologs metadata](https://github.com/baymlab/2020_Owen-BstA/blob/master/data/BstAhomologs/bsta_homologs.tsv)
  - [BstA profile-hmm](https://github.com/baymlab/2020_Owen-BstA/blob/master/data/BstAhomologs/bsta.hmm)
- Figures
  - [Figure 2.A: All annotations](https://github.com/baymlab/2020_Owen-BstA/blob/master/data/Figures/2A/annotation_counts.tsv)
  - [Figure 2.B: Annotation of genomic neighborhoods](https://github.com/baymlab/2020_Owen-BstA/tree/master/data/Figures/2B)
  - [Figure 2.C: BstA homologs alignment](https://github.com/baymlab/2020_Owen-BstA/blob/master/data/Figures/2C/bsta_alignment.afa)
- Raw data
  - The full raw data of the pipeline is generated as specified in [ipynb_notebook/run_snakemake.ipynb](ipynb_notebook/run_snakemake.ipynb).

## Code
- [`Snakefile`](https://github.com/baymlab/2020_Owen-BstA/blob/master/Snakefile): A Snakemake pipeline which does the core genomic neighborhood analysis. It consists of downloading a specified sequence region from NCBI, annotating the region with `Prokka`, running `hmmscan`, and updating the annotations with BstA and `Pfam` hits. This produces the [raw annotation data](ipynb_notebook/run_snakemake.ipynb). A [selection](https://github.com/baymlab/2020_Owen-BstA/blob/master/data/Figures/2B/selected_genomic_neighborhoods.tsv) of these results was processed to create *Figure 2.B*.
2. [`ipynb_notebook/`](https://github.com/baymlab/2020_Owen-BstA/tree/601ce14f9d81d701d49e474615e261b4d5f28230/ipynb_notebook): Contains further details about how the pipeline is run, as well as additional processing of the raw data. The annotations analysis for *Figure 2.A* can be found in [`ipynb_notebook/load_all_tables.ipynb`](https://github.com/baymlab/2020_Owen-BstA/blob/601ce14f9d81d701d49e474615e261b4d5f28230/ipynb_notebook/load_all_tables.ipynb).
