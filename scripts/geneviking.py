from Bio import Entrez
import io
from Bio import SeqIO
import pandas as pd

__author__ = 'Natalia Quinones-Olvera'
__email__ = "nquinones@g.harvard.edu"
Entrez.email = __email__


class NCBIQuery:
    
    def __init__(self, acc, start, end, thresh):
        
        # ncbi accession
        self.acc = acc
        
        # range
        self.query_range = set(range(start, end))

        self.start = start
        self.end = end
        
        # neighborhood
        self.rangeend = end + thresh
        self.rangestart = start - thresh
        
        if self.rangestart < 0:
            self.rangestart = 0
            
    
    def __parse_loc(self, loc, ref_start):
        '''
        Function to parse the location object
        of Biopython's .gb parser.
        '''
        
        # start
        if str(type(loc.start)) == "<class 'Bio.SeqFeature.ExactPosition'>":
            start = loc.start + ref_start
            start_str = start
        elif str(type(loc.start)) == "<class 'Bio.SeqFeature.BeforePosition'>":
            start = loc.start + ref_start
            start_str = '<' + str(start)
        elif str(type(loc.start)) == "<class 'Bio.SeqFeature.AfterPosition'>":
            start = loc.start + ref_start
            start_str = '>' + str(start)
        else:
            start_str = 'unknown position type'

        # end
        if str(type(loc.end)) == "<class 'Bio.SeqFeature.ExactPosition'>":
            end = loc.end + ref_start - 1
            end_str = end
        elif str(type(loc.end)) == "<class 'Bio.SeqFeature.BeforePosition'>":
            end = loc.end + ref_start - 1 
            end_str = '<' + str(end)
        elif str(type(loc.end)) == "<class 'Bio.SeqFeature.AfterPosition'>":
            end = loc.end + ref_start -1
            end_str = '>' + str(end)
        else:
            end_str = 'unknown position type'

        if loc.strand == 1:
            strand = '+'
        else:
            strand = '-'

        return [start, end, start_str, end_str, strand]
    
    
    def download_gb(self, output=None):
        '''
        Requests genbank file with Entrez
        '''
        
        epost = Entrez.epost('nuccore', id=self.acc)
        request = Entrez.read(epost)

        response = Entrez.efetch(db='nuccore',
                                 webenv=request['WebEnv'],
                                 query_key=request['QueryKey'],
                                 rettype='gb',
                                 retmode="text",
                                 seq_start=self.rangestart,
                                 seq_stop=self.rangeend)

        response_io = io.StringIO(response.read())
        
        if output:
            for record in SeqIO.parse(response_io, 'genbank'):
                genbank_records = record
            
            SeqIO.write(genbank_records, output, 'genbank')
        
        return response_io
        
    
    def download_fasta(self, output=None):
        '''
        Requests fasta file with Entrez
        '''
        epost = Entrez.epost('nuccore', id=self.acc)
        request = Entrez.read(epost)

        response = Entrez.efetch(db='nuccore',
                                 webenv=request['WebEnv'],
                                 query_key=request['QueryKey'],
                                 rettype='fasta',
                                 retmode="text",
                                 seq_start=self.rangestart,
                                 seq_stop=self.rangeend)

        response_io = io.StringIO(response.read())
        
        if output:
            for record in SeqIO.parse(response_io, 'fasta'):
                fasta_records = record
            
            SeqIO.write(fasta_records, output, 'fasta')
        
        return response_io


def hmmscan_table(file, save=False):
    '''
    '''
    # read dataframe
    df = pd.read_csv(file,
                     comment='#',
                     header=None,
                     usecols=range(22),
                     delim_whitespace=True)
    
    # label columns
    df.columns = ['target_name', 'acc', 'target_len', 'query_name', 'acc_2', 'query_len', 'full_evalue', 'full_score', 'full_bias',
                  'dom_num', 'of', 'dom_c-Evalue', 'dom_i-Evalue', 'dom_score', 'dom_bias',
                  'hmm_from', 'hmm_to',
                  'ali_from', 'ali_to',
                  'env_from', 'env_to', 'acc3']
    
    # 
    all_orfs = []
    
    for query in df['query_name'].value_counts().index:
        df_query = df[df['query_name'] == query]
        df_querybest = df_query[df_query['full_score'] == df_query['full_score'].max()]
        all_orfs.append(df_querybest)
    
    df_calls = pd.concat(all_orfs)
    df_calls = df_calls.sort_values(by=['query_name', 'dom_num'])
    df_calls['coverage'] = ((df_calls['ali_to'] - df_calls['ali_from']) / df_calls['query_len'] * 100).round(0)
    
    # pick important columns
    df_clean = df_calls[['query_name', 'query_len', 'target_name', 'acc', 'coverage', 'full_score', 'ali_from', 'ali_to', 'hmm_from', 'hmm_to']]
    
    if save:
        df.to_csv(save,
                  sep='\t',
                  index=None)
    
    return df_clean


def make_dict(df):
    '''
    '''
    annotations = {}

    for query in sorted(df['query_name'].value_counts().index):

        # subset the dataframe to specific query_name
        df_query = df[df['query_name'] == query]

        query_names = df_query['target_name'].value_counts().index
        accs = df_query['acc'].value_counts().index
        
        ali_coords = []
        scores = []
        
        for column, row in df_query.iterrows():
            loc = '{}-{}'.format(row['ali_from'], row['ali_to'])
            ali_coords.append(loc)

        if len(df_query) is not 1:
            pfam_name = ','.join(query_names)
            pfam_acc = ','.join(accs)
            model_ali = ','.join(ali_coords)
        else:
            pfam_name = query_names[0]
            pfam_acc = accs[0]
            model_ali = ali_coords[0]
            
            

        annotations[query] = {'pfam_name' : pfam_name,
                              'pfam_acc': pfam_acc,
                              'model_alignment': model_ali}
        
    return annotations

def update_gb(hmmer_table, gb_file, save_table=False, save_gb=False):
    '''
    '''
    df = hmmscan_table(hmmer_table, save_table)
    annotation_dict = make_dict(df)
    
    for gb_record in SeqIO.parse(gb_file, "genbank"):
        for feat in gb_record.features:
            if feat.type == 'CDS':
                locus = feat.qualifiers['locus_tag'][0]
                if locus in annotation_dict.keys():
                    feat.qualifiers['product'] = annotation_dict[locus]['pfam_name']
                    feat.qualifiers['protein_id'] = None
                    coords = annotation_dict[locus]['model_alignment']
                    feat.qualifiers['inference'] = f'Annotated with hmmscan. Coords {coords}'
                    feat.qualifiers['db_xref'] = 'Pfam:{}'.format(annotation_dict[locus]['pfam_acc'])
                    feat.qualifiers['gene'] = annotation_dict[locus]['pfam_name']
                else:
                    feat.qualifiers['product'] = 'hypothetical protein'
                    feat.qualifiers['protein_id'] = None
                    feat.qualifiers['inference'] = None
                    feat.qualifiers['gene'] = 'hypothetical protein'
    
    if save_gb:
        SeqIO.write(gb_record, save_gb, "genbank")
    
    return df