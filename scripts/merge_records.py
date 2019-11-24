#!//home/ghal/anaconda3/envs/intelpython3/bin/python3.6
from Bio import SeqIO
import sys

def chunkstring(string, length):
    return (string[0+i:length+i] for i in range(0, len(string), length))

def readFNA(file):
    path = f"{'/'.join(file.split('/')[:-1])}/{file.split('/')[-1].replace('fna', 'fa')}"
    sequence        = ''
    chromosomeCount = 0
    identifier      = file.split('/')[-1].replace('.fna', '')

    fastaRecords = SeqIO.parse(open(file, 'r'),'fasta')
    record       = next(fastaRecords)
    description = record.description
    if 'chromosome' in description:
        chromosomeCount += 1
        sequence = sequence + record.seq
        descList    = description.split(' ')
        descList[0] = f">{identifier}"
        fileHead = ' '.join(descList)
        fileHead = fileHead.replace('chromosome 1', 'all chromosomes')
        with open(path, 'w') as f:
            f.write(fileHead + '\n')
    for record in fastaRecords:
        description = record.description
        if 'chromosome' in description:
            sequence = sequence + record.seq
            chromosomeCount += 1

    with open("/home/ghal/git-data/Fasta/Animals/H2.fa", 'a') as f:
        for line in chunkstring(sequence, 60):
            f.write(str(line.upper()) + '\n')
if __name__ == '__main__':
    readFNA(sys.argv[1])
    #fastaRecords = SeqIO.parse(open(sys.argv[1], 'r'),'fasta')
    #for record in fastaRecords:
#        print(record.id)
