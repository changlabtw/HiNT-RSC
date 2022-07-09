from logging import raiseExceptions
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

def calculate_dict_correlation(dict1, dict2, pearson = True, spearman = True):
    if not len(dict1.keys()) == len(dict2.keys()):
        raise Exception('dict size inconsistent') 
    
    dict_list_1, dict_list_2 = [], []
    for key in chr_bin_number.keys():
        dict_list_1 = dict_list_1 + list(dict1[key])
        dict_list_2 = dict_list_2 + list(dict2[key])
    if pearson:
        pearson = round(stats.pearsonr(dict_list_1, dict_list_2)[0],3)
        print(f'Pearson corr. : ', pearson)
    if spearman:
        spearman = round(stats.spearmanr(dict_list_1, dict_list_2)[0],3)
        print(f'Spearman corr. : ', spearman)
    return pearson, spearman

def calculate_correlation_with_groundtruth(trg_seg_file, 
                                            groundtruth_file = '../../CNV/WGS/K562/BICseq2/K562_WGS_CNV_lambda4.txt'):
    gt_seg = load_seg_to_bin(groundtruth_file)
    trg_seq = load_seg_to_bin(trg_seg_file)

    gt_seg_overall, trg_seq_overall = [], []
    for key in chr_bin_number.keys():
        gt_seg_overall = gt_seg_overall + list(gt_seg[key])
        trg_seq_overall = trg_seq_overall + list(trg_seq[key])
    
    pearson = round(stats.pearsonr(gt_seg_overall, trg_seq_overall)[0],3)
    spearman = round(stats.spearmanr(gt_seg_overall, trg_seq_overall)[0],3)
    print(f'Pearson corr. : ', pearson)
    print(f'Spearman corr. : ', spearman)
    return pearson, spearman



def calculate_precision_racall(trg_seg_file, 
                                groundtruth_file = '../../CNV/WGS/K562/BICseq2/K562_WGS_CNV_lambda4.txt', 
                                copy_ratio_threshold : float = 0.3, 
                                iou_threshold : float = 0.50):
    
    # get pair dict of chr, each chr is a list of [start, end, state] (state = 'gain' or 'loss')
    gt_seg_pair_dict = load_seg_pair(groundtruth_file, threshold = copy_ratio_threshold, has_offset = True)
    trg_seg_pair_dict = load_seg_pair(trg_seg_file, threshold = copy_ratio_threshold)

    # report the number of ground truth and targer
    num_gt_pair, num_trg_pair = 0, 0
    all_gt_sizes = list()
    for chr_ in gt_seg_pair_dict.keys():
        num_gt_pair = num_gt_pair + len(gt_seg_pair_dict[chr_])
        num_trg_pair = num_trg_pair + len(trg_seg_pair_dict[chr_])
        for gt_pair in gt_seg_pair_dict[chr_]: all_gt_sizes.append(gt_pair[1]-gt_pair[0])
    print(f'Number of ground truth : {num_gt_pair}')
    print(f'Number of target : {num_trg_pair}')

    # generate seg that overlap over threshold
    def getIntercept(a, b):
        return max(0, min(a[1], b[1]) - max(a[0], b[0]))
    def getUnion(a, b):
        return max(a[1], b[1]) - min(a[0], b[0])
    def IOU(a, b):
        overlap_size = getIntercept(a,b)
        IOU_frac = getIntercept(a,b) / getUnion(a,b)
        return IOU_frac, overlap_size
    
    hit_pair_dict = dict()
    hit_pair_sizes, overlap_sizes = list(), list()
    for chr_ in trg_seg_pair_dict.keys():
        hit_pair_dict[chr_] = list()
        for pair in trg_seg_pair_dict[chr_]:
            for gt_pair in gt_seg_pair_dict[chr_]:
                overlap_frac, overlap_size = IOU(pair, gt_pair)
                if (overlap_frac > iou_threshold) and (pair[2] == gt_pair[2]):  # iou > threshold and there are in the same state (gain or loss)
                    hit_pair_dict[chr_].append(pair)
                    hit_pair_sizes.append(gt_pair[1] - gt_pair[0])
                    overlap_sizes.append(overlap_size)
    print(f'Number of hits : {len(hit_pair_sizes)}')

    # precision and recall
    num_hit_pair = 0 
    for chr_ in hit_pair_dict.keys():
        num_hit_pair = num_hit_pair + len(hit_pair_dict[chr_])
    precision = round(num_hit_pair / num_trg_pair,3)
    recall = round(num_hit_pair / num_gt_pair,3)
    print(f'precision : {precision}\nrecall : {recall}\n')
    
    # hit size by length
    print("hits by size : ")
    print("\t < 1Mb : ", len([i for i in hit_pair_sizes if i < 1000000]), ' / ', len([i for i in all_gt_sizes if i < 1000000]))
    print("\t1Mb - 2Mb : ", len([i for i in hit_pair_sizes if 1000000 <= i < 2000000]), ' / ', len([i for i in all_gt_sizes if 1000000 <= i < 2000000]))
    print("\t2Mb - 5Mb : ", len([i for i in hit_pair_sizes if 2000000 <= i < 5000000]), ' / ', len([i for i in all_gt_sizes if 2000000 <= i < 5000000]))
    print("\t5Mb - 10Mb : ", len([i for i in hit_pair_sizes if 5000000 <= i < 10000000]), ' / ', len([i for i in all_gt_sizes if 5000000 <= i < 10000000]))
    print("\t10Mb - 20Mb : ", len([i for i in hit_pair_sizes if 10000000 <= i < 20000000]), ' / ', len([i for i in all_gt_sizes if 10000000 <= i < 20000000]))
    print("\t20Mb - 50Mb : ", len([i for i in hit_pair_sizes if 20000000 <= i < 50000000]), ' / ', len([i for i in all_gt_sizes if 20000000 <= i < 50000000]))
    print("\t>= 50Mb : ", len([i for i in hit_pair_sizes if i >= 50000000]), ' / ', len([i for i in all_gt_sizes if i >= 50000000]))

    # avg. hit size and hit frac
    sum_overlap_size = round(np.sum(overlap_sizes),3)
    print(f'sum of overlap size : {sum_overlap_size}')
