import os
import numpy as np
import math

ref_genome_folder = 'Data/ref_genome/hg19'
exclude_chr = ['chrM', 'chrY','chrMT']
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


def output_bin(
    dict_ : dict, 
    trg_output_dir : str, 
    trg_name : str = '') -> None:

    if not os.path.isdir(trg_output_dir):
        os.mkdir(trg_output_dir)
    config_file = open(f'{trg_output_dir}/{trg_name}_config.txt', 'w+')
    config_file.write('chr\tpath\n')

    for key in dict_.keys():
        bin_file_path = f'{trg_output_dir}/{trg_name}_{key}.bin'
        config_file.write(f'{key}\t{bin_file_path}\n')
        with open(bin_file_path, 'w+') as f:
            f.write(f'start\tend\tobs\texp\tvar\n')
            min_value =  min(dict_[key])
            for i, value in enumerate(dict_[key]):
                start_ = i * resolution
                end_ = ( i +1 ) * resolution - 1
                obs = value - min_value
                exp = 0 - min_value
                var = obs - exp
                # only store if var == 0
                if not var == 0:
                    f.write(f'{start_}\t{end_}\t{obs}\t{exp}\t{var}\n')
    config_file.close()


