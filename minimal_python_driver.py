import ctypes
import argparse
import gzip
from mimetypes import guess_type
from functools import partial
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
import multiprocessing as mlp
import time
import os
from datetime import datetime

global path

lib = ctypes.cdll.LoadLibrary("./build/libvlmc.so")


## Importing VLMC library training fucntion calls.
train = lib.train_
train.argtypes = [
    ctypes.c_char_p, # ID
    ctypes.c_char_p, # Sequence
    ctypes.c_size_t, # Max Depth
    ctypes.c_size_t, # Min Count
    ctypes.c_float,  # Threshold
    ctypes.c_size_t, # Number of Parameters
    ctypes.c_char_p, # Pruning Method
    ctypes.c_char_p, # Estimator
    ctypes.c_bool,   # Multi-core Execution
    ctypes.c_size_t  # Split Depth

]
train.restype = ctypes.c_char_p



def setPath(identifier):
    global path
    dateTimeObj = datetime.now()
    file        = f'{identifier.decode()}_{dateTimeObj.hour}:{dateTimeObj.minute}:{dateTimeObj.second}'
    resultPath  = f"./Results/{identifier.decode()}"
    if not os.path.isdir(resultPath):
        os.mkdir(resultPath)
    path = f"{resultPath}/{str(file)}.tree"
    return path


def train_call(arg, identifiers, sequence):
    train(identifiers,
            sequence,
            arg.max_depth,
            arg.min_count,
            arg.threshold,
            arg.number_of_parameters,
            arg.pruning_method.encode(),
            arg.estimator.encode(),
            arg.multi_core,
            arg.split_depth)

def PSTConstruction(arg):
    sequence         = ''
    identifiers      = []
    descriptions     = []
    chromosome_count = 0
    scaffold_count   = 0
    other_rec_count  = 0
    encoding = guess_type(arg.file)[1]  # uses file extension
    _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

    start_r = time.time()
    with _open(arg.file) as f:
        #for record in SeqIO.parse(f, 'fasta'):
        for value in SimpleFastaParser(f):
            id, seq = value
            identifiers.append(id.split(' ')[0])
            descriptions.append(id)
            sequence = sequence + seq
    stop_r = time.time()

    for desc in descriptions:
        if 'chromosome' in desc:
            chromosome_count += 1
        elif 'scaffold' in desc:
            scaffold_count   += 1
        else:
            other_rec_count  +=1

    # Example: adapt to extract features you are interested in
    if not arg.quiet:
        print('----------------------------------------------------------')
        print('Processing the file:\n{}'.format(arg.file))
        #print(f'Its description is: \n{descriptions[0]}')
        print(f'The genome assembly contains {chromosome_count} chromosomes, {scaffold_count} scaffolds, and {other_rec_count} other records.')
        print(f'Length of combined sequences is {len(sequence)} nucleotides.')
        print(f'Reading the file took {stop_r - start_r} seconds.')
        print('----------------------------------------------------------')
        print('Generating PST...')
    else:
        print('Tree output quieted.')

    args = (arg, identifiers[0].encode("utf-8"), sequence.encode("utf-8"))
    process = mlp.Process(target=train_call, args=args)              # Start a worker processes.
    process.start()

    sequence = ''
    args     = ''
    process.join()
    return
    #if arg.outputfile:
    #    setPath(identifier)
    #    with open(path, 'w') as f:
    #        f.write("%s\n" % arg)
    #        for line in PST.decode().split('\n'):
    #            f.write("%s\n" % line)

    #if not arg.quiet:
    #    print(PST.decode())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Wrapper for paralell PST construction',
        epilog = 'Example: python vlmc.py Data/Fasta/Virus/NC_001348.1.fa\
                  Command: Generates PST for NC_001348.1, with flags -o for file output,\
                  min-count 10 and number-of-parameters 25.'
    )

    parser.add_argument('file',
                        help = 'FASTA file input for PST construction.' )
    parser.add_argument('--max-depth', '-d',
                        default = '15',
                        dest    = 'max_depth',
                        metavar = '(unsigned 64 bit integer)',
                        type    =  int,
                        help    = 'Max depth of the built probabilistic suffix \
                                   tree.')

    parser.add_argument('--min-count', '-c',
                        default = '100',
                        dest    = 'min_count',
                        type    =  int,
                        metavar = '(unsigned 64 bit integer)',
                        help    = 'Minimum number of time each node/context \
                                   has to appear in the string to be included \
                                   in the probabilistic suffix tree.')

    parser.add_argument('--threshold', '-k',
                        default = '1.2',
                        dest    = 'threshold',
                        type    =  float,
                        metavar = '(float)',
                        help    = 'Threshold for the pruning stage of the \
                                   algorithm. Smaller value gives larger tree.')

    parser.add_argument('--number-of-parameters', '-n',
                        default = '192',
                        dest    = 'number_of_parameters',
                        type    =  int,
                        metavar = '192',
                        help    = "For the 'parameters' puning-method, the \
                                   number of parameters to prune the tree \
                                   until.")

    parser.add_argument('--estimator', '-e',
                        default = 'KL',
                        dest    = 'estimator',
                        type    = str,
                        metavar = '(string)',
                        help    = "Estimator used to determine which states \
                                   should be  pruned. Either 'KL' or 'PS'.")

    parser.add_argument('--pruning-method', '-p',
                        default = 'cutoff',
                        dest    = 'pruning_method',
                        type    = str,
                        metavar = '(string)',
                        help    = "pruning method to use. Either 'cutoff' for \
                                   pruning until the threshold is reached, or \
                                   'parameters' to prune until a certain number\
                                   of parameters have been reached.")

    parser.add_argument('--quiet', '-q',
                        default = False,
                        dest    = 'quiet',
                        action  = 'store_true',
                        help    = 'Use this flag to quiet terminal output.')

    parser.add_argument('--output-file', '-of',
                        default = False,
                        dest    = 'outputfile',
                        action  = 'store_true',
                        help    = 'Use this flag to generate output file into \
                                   ./Data/Result/identifier(HH:mm:ss).tree ')

    parser.add_argument('--multi-core', '-m',
                        default = False,
                        dest    = 'multi_core',
                        action  = 'store_true',
                        help    = 'Use this flag to generate output utilising \
                                   multiple cores. ')

    parser.add_argument('--split-depth', '-s',
                        default = 1,
                        dest    = 'split_depth',
                        metavar = f'{mlp.cpu_count()}',
                        type    = int,
                        help    = 'If using paralelll option, how many processors \
                                   to utilise. Higher core count might lead to lower \
                                   accuracy.')

    args = parser.parse_args()
    start = time.time()
    PSTConstruction(args)
    end = time.time()
    print('----------------------------------')
    print(f'Duration: {round(end-start,3)}s')
    print('----------------------------------')
