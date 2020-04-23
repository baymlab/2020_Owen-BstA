import os
import sys
import pandas as pd
from scripts import geneviking as gv

# configfile: 'config.yaml'

# read dataframe
df = pd.read_table(config['input_file']).set_index("acc",
                                                   drop=False)

df['files'] = df['acc'] + '_' + df['start'].astype(str) + '-' + df['end'].astype(str)

# pipeline
rule all:
    input:
        expand('{direc}/viking/gb/{sample}.gbk', sample=df['files'], direc=config['output_dir'])

rule download:
    # define wildcards
    output: '{direc}/fasta/{acc}_{start}-{end}.fasta'
    
    run:
        # make query object
        query = gv.NCBIQuery(wildcards.acc,
                             int(wildcards.start),
                             int(wildcards.end),
                             config['params']['neighborhood'])
        
        if config['params']['download_fmt'] == 'fasta':
            query.download_fasta(output=f'{wildcards.direc}/fasta/{wildcards.acc}_{wildcards.start}-{wildcards.end}.fasta')
        
        else:
            print('Unrecognized format. Please indicate <fasta>.')

rule prokka:
    output:
        '{direc}/prokka/{acc}_{start}-{end}/{acc}_{start}-{end}.gbk',
        '{direc}/prokka/{acc}_{start}-{end}/{acc}_{start}-{end}.faa'
    
    input: '{direc}/fasta/{acc}_{start}-{end}.fasta'
    
    shell:
        "prokka --outdir {wildcards.direc}/prokka/{wildcards.acc}_{wildcards.start}-{wildcards.end} --force --quiet --prefix {wildcards.acc}_{wildcards.start}-{wildcards.end} --compliant --kingdom Bacteria {wildcards.direc}/fasta/{wildcards.acc}_{wildcards.start}-{wildcards.end}.fasta"

rule hmmer:
    output:
        '{direc}/hmmscan/{acc}_{start}-{end}/{acc}_{start}-{end}.out',
        '{direc}/hmmscan/{acc}_{start}-{end}/{acc}_{start}-{end}.tbl'
    
    input:
        '{direc}/prokka/{acc}_{start}-{end}/{acc}_{start}-{end}.faa'
    
    shell:
        "hmmscan -o {wildcards.direc}/hmmscan/{wildcards.acc}_{wildcards.start}-{wildcards.end}/{wildcards.acc}_{wildcards.start}-{wildcards.end}.out --domtblout {wildcards.direc}/hmmscan/{wildcards.acc}_{wildcards.start}-{wildcards.end}/{wildcards.acc}_{wildcards.start}-{wildcards.end}.tbl data/pfam/Pfam-BstA.hmm {wildcards.direc}/prokka/{wildcards.acc}_{wildcards.start}-{wildcards.end}/{wildcards.acc}_{wildcards.start}-{wildcards.end}.faa"
        
rule replace_annotations:
    output:
        '{direc}/viking/gb/{acc}_{start}-{end}.gbk',
        '{direc}/viking/tables/{acc}_{start}-{end}.tsv'
    input:
        '{direc}/prokka/{acc}_{start}-{end}/{acc}_{start}-{end}.gbk',
        '{direc}/hmmscan/{acc}_{start}-{end}/{acc}_{start}-{end}.tbl'
    run:
        gv.update_gb(f'{wildcards.direc}/hmmscan/{wildcards.acc}_{wildcards.start}-{wildcards.end}/{wildcards.acc}_{wildcards.start}-{wildcards.end}.tbl',
                     f'{wildcards.direc}/prokka/{wildcards.acc}_{wildcards.start}-{wildcards.end}/{wildcards.acc}_{wildcards.start}-{wildcards.end}.gbk',
                     save_gb= f'{wildcards.direc}/viking/gb/{wildcards.acc}_{wildcards.start}-{wildcards.end}.gbk',
                     save_table= f'{wildcards.direc}/viking/tables/{wildcards.acc}_{wildcards.start}-{wildcards.end}.tsv')
        