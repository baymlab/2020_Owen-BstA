# 2020_Owen-BstA

## About
Code corresponing to the paper [Prophage-encoded abortive infection proteins that encode cognate self-immunity]()

## Data
- [BstA homologs]()
  - [Sequences]()
  - [Homologs metadata]()
- [BstA profile-hmm]()
- [Raw annotation results]()
- [Processed results]()
  - [Figure 2.A: Annotation counts]() <!--TO_DO: annot tables and counts-->
  - [Figure 2.B: Annotation of genomic neighborhoods]() <!--point at prokka+viking gb-->
  - [Figure 2.C: BstA homologs alignment]() <!--point at alignment file -->

## Code
- [`Snakefile`](https://github.com/baymlab/2020_Owen-BstA/blob/master/Snakefile): A Snakemake pipeline which does the core genomic neighborhood analysis. It consists of downloading a specified sequence region from NCBI, annotating the region with `Prokka`, running `hmmscan`, and updating the annotations with BstA and `Pfam` hits. This produces the [raw annotation results](). A [selection]() of these results was processed to create *Figure 2.B*.
2. [`ipynb_notebook`](https://github.com/baymlab/2020_Owen-BstA/tree/601ce14f9d81d701d49e474615e261b4d5f28230/ipynb_notebook): Contains additional processing of the raw data. The annotations counts for *Figure 2.A* can be found in [`ipynb_notebook/load_all_tables.ipynb`](https://github.com/baymlab/2020_Owen-BstA/blob/601ce14f9d81d701d49e474615e261b4d5f28230/ipynb_notebook/load_all_tables.ipynb).
