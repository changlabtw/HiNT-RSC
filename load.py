import pandas as pd
import hicstraw
import numpy as np
import math
import os
import tabix

ref_genome_folder = 'Data/ref_genome/hg19'
exclude_chr = ['chrM', 'chrY']
contact_matrix = 'ContactMatrix/GSE63525_K562_combined.hic'
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

def generate_coverage_profile_intra_inter(
    contact_matrix_uri, 
    trg_chr : str = 'chr1', 
    resolution : int = 50000) -> dict:

    chr_name_list = list(chr_length.keys())
    bin_num = math.floor(chr_length[trg_chr]/resolution) + 1
    coverage_profile = np.zeros(bin_num)

    for chr_ in chr_length.keys():
        
        # get contacts
        result = hicstraw.straw("observed", 'NONE', contact_matrix_uri, trg_chr.split('chr')[1], chr_.split('chr')[1], 'BP', resolution)

        # build 1D coverage profile
        for i in range(len(result)):
            # binX binY will be lower index -> higher index of chr (no matter how you place in straw.straw)
            if chr_name_list.index(trg_chr) <= chr_name_list.index(chr_):
                trg_bin = math.floor(result[i].binX / resolution)
                another_bin = math.floor(result[i].binY / resolution)
            else:
                another_bin = math.floor(result[i].binX / resolution)
                trg_bin = math.floor(result[i].binY / resolution)

            count = 0 if np.isnan(result[i].counts) else result[i].counts

            # straw.straw only gives upper triangular matrix
            if (trg_chr == chr_) & (trg_bin != another_bin):
                coverage_profile[trg_bin] += count
                coverage_profile[another_bin] += count
            else:
                coverage_profile[trg_bin] += count
    return coverage_profile


def generate_coverage_profile_intra(
    contact_matrix_uri, 
    trg_chr : str = 'chr1', 
    resolution : int = 50000) -> dict:

    bin_num = math.floor(chr_length[trg_chr]/resolution) + 1
    coverage_profile = np.zeros(bin_num)
        
    # get contacts
    result = hicstraw.straw("observed", 'NONE', contact_matrix_uri, trg_chr.split('chr')[1], trg_chr.split('chr')[1], 'BP', resolution)

    # build 1D coverage profile
    for i in range(len(result)):
        # binX binY will be lower index -> higher index of chr (no matter how you place in straw.straw)
        trg_bin = math.floor(result[i].binX / resolution)
        another_bin = math.floor(result[i].binY / resolution)
        count = 0 if np.isnan(result[i].counts) else result[i].counts

        # straw.straw only gives upper triangular matrix
        if (trg_bin != another_bin):
            coverage_profile[trg_bin] += count
            coverage_profile[another_bin] += count
        else:
            coverage_profile[trg_bin] += count
    
    return coverage_profile


def load_GC_comtent(GC_uri) -> dict:
    GC_dict, counts = dict(), dict()
    for key in chr_length.keys():
        GC_dict[key] = np.zeros(chr_bin_number[key])
        counts[key] = np.zeros(chr_bin_number[key])

    with open(GC_uri, 'r') as f:
        for line in f.readlines()[1:]:
            line_split = line.strip().split()
            
            if line_split[0] in exclude_chr:
                continue
            else:
                start = math.floor(int(line_split[1]) / resolution)
                GC_dict[line_split[0]][start] += float(line_split[3])
                counts[line_split[0]][start] += 1

    for key in GC_dict.keys():
        GC_dict[key] = GC_dict[key] / counts[key]
    
    return GC_dict


# function from HiNT source code
def getmappability(mappablity_track, matrixchrom, matrixstart, matrixend, resolution):
	tb = tabix.open(mappablity_track)
	records = tb.query(matrixchrom, int(matrixstart), int(matrixend))
	Allstart = int(matrixstart)
	Allend = int(matrixend)
	values = []
	binvalues = []
	values += [0]*(Allstart%resolution)
	for r in records:
		chrom,start,end,value = r
		length = int(end) - int(start)
		if int(start) == Allstart:
			temp = [float(value)]*length
		else:
			temp = [0]*(int(start)-Allstart) + [float(value)]*length
		values += temp
		Allstart = int(end)
	if Allstart == Allend:
		pass
	else:
		values += [0]*(Allend-Allstart)
	values += [0]*(resolution-Allend%resolution)
	tempcount = 0
	for i,v in enumerate(values):
		if i%resolution == (resolution-1):
			binvalues.append(float(tempcount)/float(resolution))
			tempcount = 0
		else:
			tempcount += v
	return binvalues

# function from HiNT source code
def getGCpercent(gcfile,chrom,resolution):
	GCdata = pd.read_csv(gcfile, header = 0,sep="\t",index_col=0)
	GCper = GCdata.loc[chrom:chrom,'5_pct_gc']
	end = GCdata.loc[chrom:chrom,'3_usercol']
	length = end[-1:][0]
	steps = int(length/resolution)
	GC = []
	GC += [0]*(steps+1)
	for i in range(steps):
		binstart = i*(int(resolution/1000))
		binend = (i+1)*(int(resolution/1000))
		averageGC = np.average(GCper[binstart:binend])
		GC[i] = averageGC
	lastbinStart = steps*(int(resolution/1000))
	lastbinEnd = len(GCper)
	GC[-1] = np.average(GCper[lastbinStart:lastbinEnd])
	#print GCper
	GC = np.asarray(GC)

	return GC

# function from HiNT source code
def getFragmentsNumber(sitesInfo, matrixchrom, resolution,chromlength):
	sites = sitesInfo[matrixchrom]
	bincounts = [0.0] * (int(int(chromlength)/int(resolution)) + 1)
	for site in sites:
		idx = int(int(site) / int(resolution))
		if idx >= len(bincounts):
			pass
		else:
			bincounts[idx] += 1.0
	return bincounts


def load_num_res_site(num_res_site_uri) -> dict:
    num_res_sites_dict = dict()
    for key in chr_length.keys() :
        num_res_sites_dict[key] = np.zeros(chr_bin_number[key])

    with open(num_res_site_uri, 'r') as f:
        for line in f.readlines():
            line_split = line.split()
            if line_split[0] in exclude_chr:
                continue
            else:
                for site in line_split[1:]:
                    bin_ = math.floor( int(site) / resolution )
                    num_res_sites_dict[line_split[0]][bin_] += 1
    
    return num_res_sites_dict


def load_seg_to_bin(seg_file) -> dict:
    cnv_dict = dict()
    for key in chr_bin_number.keys():
        cnv_dict[key] = np.zeros(chr_bin_number[key])

    with open(seg_file, 'r') as f:
        has_p_value = True if len(f.readline().strip().split()) == 8 else False
        for line in f.readlines():
            if has_p_value:
                chr_, start, end, binNum, observed, expected, copy_ratio, p_value = line.strip().split()
            else:
                chr_, start, end, binNum, observed, expected, copy_ratio = line.strip().split()
            start_bin = math.floor(int(start) / resolution)
            end_bin = math.floor(int(end) / resolution)
            for bin_num in range(start_bin, end_bin+1):
                if bin_num >= chr_bin_number[chr_]:
                    break
                cnv_dict[chr_][bin_num] = float(copy_ratio)
    return cnv_dict


def load_seg_pair(
    seg_file, 
    threshold :float = 0.3, 
    has_offset : bool = False) -> dict:
    
    cnv_pair_dict = dict()
    for key in chr_bin_number.keys():
        cnv_pair_dict[key] = list()
    
    with open(seg_file, 'r') as f:
        has_p_value = True if len(f.readline().strip().split()) == 8 else False
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
            if abs(float(copy_ratio)) > threshold :
                cnv_pair_dict[chr_].append([start, end])
    return cnv_pair_dict


def load_bin_dir_to_ratio(bin_dir) -> dict :
   
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
                bin_ = math.floor(int(start_)/resolution)
                if  not ( (float(exp) == 0) or (float(obs) == 0) ):
                    log2_ratio = math.log2(float(obs)/float(exp))
                else:
                    log2_ratio = 0
                bin_dict[chr_][bin_] = log2_ratio
    return bin_dict