import os
import numpy as np
import pandas as pd
import math
import argparse


ref_genome_folder = '../Data/ref_genome/hg19'
exclude_chr = ['chrM', 'chrY', 'chrMT']
resolution = 50000


# load ref genome length info
chr_length = dict()
with open(ref_genome_folder + '/hg19.len', 'r') as f:
    lines = f.read().splitlines()
    for line in lines:
        chr_name , len_ = line.split('\t')[:2]
        if chr_name in exclude_chr:
            continue
        chr_length[chr_name] = int(len_)

# calculate bin number in each chr
chr_bin_number = dict()
for key in chr_length.keys():
    chr_bin_number[key] = math.floor(chr_length[key] / resolution) + 1



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin_dir", help = "folder of bin file", type = str)
    parser.add_argument("--output", help = "output .bw file path", type = str)
    args = parser.parse_args()

    bin_dict = dict()
    for key in chr_bin_number.keys():
        bin_dict[key] = np.zeros(chr_bin_number[key])

    bin_dir = args.bin_dir
    for file in sorted(os.listdir(bin_dir)):
        if (file[0] == '.') or not(file.split('.')[-1] == 'bin'):
            continue
        chr_ = file.split('.')[0].split('_')[-1]
        with open(f'{bin_dir}/{file}', 'r') as f:
            for line in f.readlines()[1:]:
                start_, end_, obs, exp ,var = line.strip().split()
                residual = float(obs) - float(exp)
                bin_ = math.floor(int(start_)/resolution)
                bin_dict[chr_][bin_] = residual
    
    # Create a new directory because it does not exist 
    dir_ = '/'.join(args.output.split('/')[:-1])
    if not os.path.exists(dir_):
        os.makedirs(dir_)

    with open(args.output, 'w+') as f:
        for key in bin_dict.keys():
            for bin_, residual in enumerate(bin_dict[key]):
                bin_start = bin_ * 50000
                bin_end = (bin_+1) * 50000 if (bin_+1) * 50000 <= chr_length[key] else chr_length[key]
                residual = round(residual,4)
                f.write(f'{key}\t{bin_start}\t{bin_end}\t{residual}\n')