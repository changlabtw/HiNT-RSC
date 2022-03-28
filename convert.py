import os
import numpy as np
import pandas as pd
import math
import argparse


ref_genome_folder = 'Data/ref_genome/hg19'
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


def convert_bin_to_bw(bin_dir, output):

    # read each bin file in bin_dir
    bin_dict = dict()
    for key in chr_bin_number.keys():
        bin_dict[key] = np.zeros(chr_bin_number[key])
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
    dir_ = '/'.join(output.split('/')[:-1])
    if not (os.path.exists(dir_) or dir_ == None):
        os.makedirs(dir_)

    # convert to bw line and store
    with open(output, 'w+') as f:
        for key in bin_dict.keys():
            for bin_, residual in enumerate(bin_dict[key]):
                bin_start = bin_ * 50000
                bin_end = (bin_+1) * 50000 if (bin_+1) * 50000 <= chr_length[key] else chr_length[key]
                residual = round(residual,4)
                f.write(f'{key}\t{bin_start}\t{bin_end}\t{residual}\n')


def convert_seg_to_bed(seg_file, 
                        output_bed_path, 
                        copy_ratio_threshold : float = 0.3, 
                        has_offset : bool = False):
    with open(seg_file, 'r') as f:
        # check if seg file has p-value column
        line_split = f.readline().strip().split()
        has_p_value = True if len(line_split) == 8 else False

        # read seg and save to bed
        output_bed = open(output_bed_path, 'w+')
        output_bed.write(f'chrom\tstart\tend\n')
        for line in f.readlines():
            if has_p_value:
                chr_, start, end, binNum, observed, expected, copy_ratio, p_value = line.strip().split()
            else:
                chr_, start, end, binNum, observed, expected, copy_ratio = line.strip().split()
            
            # if the seg file has offset (every chr start from 1 ), then offset back
            if has_offset:
                start = int(start) - 1 if not (int(start) - 1) < 0 else 0
                end = int(end)-1 if not (int(end)-1) > chr_length[chr_] else chr_length[chr_]
            else:
                start = int(start) if not int(start) < 0 else 0
                end = int(end) if not int(end) > chr_length[chr_] else chr_length[chr_]

            # keep only significant copy ratio
            if abs(float(copy_ratio)) >= copy_ratio_threshold:
                start = start if not int(start) < 0 else 0
                end = end if not int(end) > chr_length[chr_] else chr_length[chr_]
                output_bed.write(f'{chr_}\t{start}\t{end}\n')
            else:
                continue

        output_bed.close()