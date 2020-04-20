from Bio import Entrez
import io
from Bio import SeqIO

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
