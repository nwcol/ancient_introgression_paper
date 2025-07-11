
import pandas 
import numpy as np
import sys


def main():

    if len(sys.argv) == 2:
        chrom = sys.argv[1]
        scale = 5000000
    else:
        chrom = sys.argv[1]
        scale = int(sys.argv[2])
    table = pandas.read_csv('/home/nick/Projects/dplus/data/cytoBand.txt', 
        sep='\t', names=['chrom', 'start', 'end', 'name', 'stain'])
    df = table[table['chrom'] == chrom]
    start = 0
    end = np.array(df[df['stain'] == 'acen']['start'])[0]
    ps = (np.linspace(start, end, (end-start) // scale) // 100000) * 100000 + 1
    start = np.array(df[df['stain'] == 'acen']['end'])[-1]
    end = np.array(df['end'])[-1]
    qs = (np.linspace(start, end, (end-start) // scale) // 100000) * 100000 + 1
    with open(f'windows.{chrom}.txt', 'w') as fout:
        for start, end in zip(ps[:-1], ps[1:]):
            fout.write(f'{int(start)}\t{int(end)}\t{int(ps[-1])}\n')
        for start, end in zip(qs[:-1], qs[1:]):
            fout.write(f'{int(start)}\t{int(end)}\t{int(qs[-1])}\n')


main()
