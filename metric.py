import numpy as np
import math
from scipy import stats
from load import load_seg_to_bin, load_seg_pair

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



def calculate_correlation_with_groundtruth(trg_seg_file, 
                                            groundtruth_file = '../../CNV/K562_WGS_CNV/K562_WGS_CNV_lambda4.txt'):
    gt_seg = load_seg_to_bin(groundtruth_file)
    trg_seq = load_seg_to_bin(trg_seg_file)

    gt_seg_overall, trg_seq_overall = [], []
    for key in chr_bin_number.keys():
        gt_seg_overall = gt_seg_overall + list(gt_seg[key])
        trg_seq_overall = trg_seq_overall + list(trg_seq[key])
    
    pearson = round(stats.pearsonr(gt_seg_overall, trg_seq_overall)[0],3)
    soearman = round(stats.spearmanr(gt_seg_overall, trg_seq_overall)[0],3)
    print(f'Pearson corr. : ', pearson)
    print(f'Spearman corr. : ', soearman)
    return pearson, soearman



def calculate_precision_racall(trg_seg_file, 
                                groundtruth_file = 'CNV/K562_WGS_CNV/K562_WGS_CNV_lambda4.txt', 
                                copy_ratio_threshold : float = 0.3, 
                                iou_threshold : float = 0.50):
    gt_seg_pair_dict = load_seg_pair(groundtruth_file, threshold = copy_ratio_threshold, has_offset = True)
    trg_seg_pair_dicr = load_seg_pair(trg_seg_file, threshold = copy_ratio_threshold)

    # report the number of ground truth and targer
    num_gt, num_trg = 0, 0
    for chr_ in gt_seg_pair_dict.keys():
        num_gt = num_gt + len(gt_seg_pair_dict[chr_])
        num_trg = num_trg + len(trg_seg_pair_dicr[chr_])
    print(f'Number of ground truth : {num_gt}')
    print(f'Number of target : {num_trg}')

    # generate seg that overlap over threshold
    def getIntercept(a, b):
        return max(0, min(a[1], b[1]) - max(a[0], b[0]))
    def getUnion(a, b):
        return max(a[1], b[1]) - min(a[0], b[0])
    def IOU(a, b):
        return getIntercept(a,b) / getUnion(a,b)
    hit_pair_dict = dict()
    overlap_frac = list()
    for key in chr_bin_number.keys():
        hit_pair_dict[key] = list()
    for chr_ in trg_seg_pair_dicr.keys():
        for pair in trg_seg_pair_dicr[chr_]:
            for gt_pair in gt_seg_pair_dict[chr_]:
                overlap = IOU(pair, gt_pair)
                if overlap >= iou_threshold:
                    hit_pair_dict[chr_].append(pair)
                    overlap_frac.append(overlap)
    print(f'Number of hits : {len(overlap_frac)}')

    # precision and recall
    num_gt_pair = 0
    for chr_ in gt_seg_pair_dict.keys():
        num_gt_pair = num_gt_pair + len(gt_seg_pair_dict[chr_]) 
    num_trg_pair = 0
    for chr_ in trg_seg_pair_dicr.keys():
        num_trg_pair = num_trg_pair + len(trg_seg_pair_dicr[chr_])
    num_hit_pair = 0 
    for chr_ in hit_pair_dict.keys():
        num_hit_pair = num_hit_pair + len(hit_pair_dict[chr_])
    precision = round(num_hit_pair / num_trg_pair,3)
    recall = round(num_hit_pair / num_gt_pair,3)
    print(f'precision : {precision}\nrecall : {recall}')

    # avg. hit size
    avg_overlap_frac = round(np.mean(overlap_frac),3)
    print(f'average hit fraction : {avg_overlap_frac}')

