import os
import numpy as np
import math


def ref_len_file():
    """Chromosome-length file of the reference genome.

    HINT_RSC_REF is either the reference root (default 'data/ref_genome', relative to paper/; the genome is then
    chosen with HINT_RSC_GENOME, default 'hg19') or a genome folder such as 'data/ref_genome/hg38'.
    The genome folder contains <genome>.len (tab-separated: chromosome, length).
    """
    folder = os.environ.get('HINT_RSC_REF', 'data/ref_genome')
    name = os.path.basename(os.path.normpath(folder))
    if not os.path.exists(os.path.join(folder, name + '.len')):
        name = os.environ.get('HINT_RSC_GENOME', 'hg19')
        folder = os.path.join(folder, name)
    return os.path.join(folder, name + '.len')


exclude_chr = ['chrM', 'chrY','chrMT']
resolution = 50000

# load ref genome length info
chr_length = dict()
with open(ref_len_file(), 'r') as f:
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


