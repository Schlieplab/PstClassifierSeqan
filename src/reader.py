#!/home/ghal/anaconda3/envs/intelpython3/bin/python3.6
from Bio import SeqIO
import gzip
import sys


def main(path):
    Sequence = ''
    fastaRecords = SeqIO.parse(file_path,'fasta')
    for record in fastaRecords:
        if 'chromosome' in record.description:
            print(record.description)
            print(record.id)
            Sequence = Sequence + record.seq
    return Sequence
if __name__ == '__main__':
    try:
        file_path = sys.argv[1]
    except:
        file_path = "/home/ghal/git/DATX05/PyPSTclassifierSeqan/Data/Fasta/Animals/GCA_002180035.3_HG00514_prelim_3.0_genomic.fna";
    return main(file_path)
