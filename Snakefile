import os
import sys
import pandas as pd
from scripts import geneviking as gv

configfile: 'config.yaml'

# read dataframe
df = pd.read_table(config['input_file']).set_index("acc",
                                                   drop=False)

df['files'] = df['acc'] + '_' + df['start'].astype(str) + '-' + df['end'].astype(str)

# pipeline
rule all:
    input: expand('{direc}/{sample}.fasta', sample=df['files'], direc=config['output_dir'])

rule download:
    # define wildcards
    output: '{direc}/{acc}_{start}-{end}.fasta'
    
    run:
        # make query object
        query = gv.NCBIQuery(wildcards.acc,
                             int(wildcards.start),
                             int(wildcards.end),
                             config['params']['neighborhood'])
        
        if config['params']['download_fmt'] == 'fasta':
            query.download_fasta(output=f'{wildcards.direc}/{wildcards.acc}_{wildcards.start}-{wildcards.end}.fasta')
        
        elif config['params']['download_fmt'] == 'genbank':
            query.download_gb(output=f'{wildcards.direc}/{wildcards.acc}_{wildcards.start}-{wildcards.end}.fasta')
        
        else:
            print('Unrecognized format. Please indicate <fasta> or <genbank>.')