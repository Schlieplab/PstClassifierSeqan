from time import sleep
import random
import resource
import multiprocessing
import sys

global k

def chunkIt(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out

def genNumbers(w):
    global k
    sequence = random.choices(['A', 'C', 'G', 'T'], weights=w, k=k)
    return sequence

def merge_subsequences(result):
    sequence = "".join(result)
    return sequence

def pool_handler(weights):

    with multiprocessing.Pool(multiprocessing.cpu_count()) as p:
        result =  p.map(genNumbers, ([w]*multiprocessing.cpu_count()))

    with multiprocessing.Pool(multiprocessing.cpu_count()) as p:
        res =  p.map(merge_subsequences, (result))
    return res

def main(weights):
    sequence = pool_handler(weights)
    tmp = []
    for s in sequence:
        for i in s:
            tmp.append(i)
    sequence = "".join(tmp)
    sequence = [sequence[x:x+60] for x in range(0, len(sequence), 60)]
    DNAnoise = f">XD_424242.42 NOISE chromosome X, Weights: {weights}\n"
    with open(f'noise{weights}.fa', 'w') as f:
            f.write(DNAnoise)
            for s in sequence:
                f.write(f'{s}\n')

if __name__ == '__main__':
    weights = [[0.1, 0.4, 0.3, 0.2], [0.7, 0.1, 0.1, 0.1], [0.3, 0.3, 0.3, 0.1]]#, [0.15, 0.15, 0.4, 0.3]]                                                               ]
    global k
    k = 1666666
    for w in weights:
        main(w)
