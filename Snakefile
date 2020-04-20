import os
import sys
import pandas as pd
from scripts import geneviking as gv

configfile: "config.yaml"

# read dataframe
df = pd.read_table(config['input_file']).set_index("acc",
                                                   drop=False)


SAMPLES = df['acc'] + '_' + df['start'].astype('str') + '-' + df['end'].astype('str')

# pipeline
rule all:
    input: expand('{sample}.fasta', sample=SAMPLES)

rule download:
    output: "{ncbiacc}.fasta"
    run: 
        print(config['output_dir'])
        acc = wildcards.ncbiacc.split('_')[0]
        coord = wildcards.ncbiacc.split('_')[1]
        start = int(coord.split('-')[0])
        end = int(coord.split('-')[1])
        
        query = gv.NCBIQuery(acc, start, end, 10)
        print(query.acc, query.rangestart, query.rangestart)
        query.download_fasta(output=f'{wildcards.ncbiacc}.fasta')