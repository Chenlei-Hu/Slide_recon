"""
analysis of fiducial diffusion sequencing result without bead barcode matching
input fastq file
output collasped barcode information
"""

import os
import time
import gzip
import argparse
import numpy as np
import mappy as mp
import pandas as pd
import editdistance
import matplotlib.pyplot as plt
from collections import Counter
from multiprocessing import Pool
from umi_tools import UMIClusterer
from scipy.signal import find_peaks
from bead_matching import barcode_matching
from scipy.ndimage import gaussian_filter1d
from matplotlib.backends.backend_pdf import PdfPages


def barcode_extract(fq1_file, fq2_file, r2type):
    """
    input: fastq file
    output: dict of barcode without matching
    """
    aln_dict = {}
    R1_bc_list = []
    R2_bc_list = []
    alignment_stat = Counter()
    for fq1, fq2 in zip(mp.fastx_read(fq1_file, read_comment=False), mp.fastx_read(fq2_file, read_comment=False)):
        alignment_stat["total_reads"] += 1

        if alignment_stat["total_reads"] % 1000000 == 0:
            print(alignment_stat["total_reads"])
    
        # check read length
        if len(fq1[1])<46 or len(fq2[1])<46:  
            alignment_stat["Read_too_short"] += 1
            continue
        R1_bc = fq1[1][:8]+fq1[1][26:32]
        R1_bumi = fq1[1][33:41] # if R1: 32:40, if V8 or V10: 33:41
        R1_UP = fq1[1][8:26]
        R1_poly = fq1[1][42:]

        # different Read 2 bead type
        if r2type == 'V9':
            R2_bc = fq2[1][:8]+fq2[1][26:32]
            R2_bumi = fq2[1][33:41]
            R2_UP = fq2[1][8:26]
            R2_poly = fq2[1][42:]
             # check UP site with mismatch < 3bp
            if editdistance.eval(R1_UP,'TCTTCAGCGTTCCCGAGA')>3 or editdistance.eval(R2_UP,'TCTTCAGCGTTCCCGAGA')>3:
                alignment_stat["UP_not_matched"] += 1
                continue 
        elif r2type == 'V15':
            R2_bc = fq2[1][:14]
            R2_bumi = fq2[1][25:33]
            R2_UP = fq2[1][15:25]
            R2_poly = fq2[1][34:]
            # check UP site with mismatch < 3bp
            if editdistance.eval(R1_UP,'TCTTCAGCGTTCCCGAGA')>3 or editdistance.eval(R2_UP,'CTGTTTCCTG')>2:
                alignment_stat["UP_not_matched"] += 1
                continue 

       
        # check if end up with polyA or polyT
        # if R1_poly.count('T')/len(R1_poly) < 0.75 or R2_poly.count('A')/len(R1_poly) < 0.75:
        #     alignment_stat["lack_of_polyAT"] += 1
        #     continue

        alignment_stat['effective_read'] += 1
        aln_dict.setdefault(R1_bc,[]).append((R2_bc,R1_bumi,R2_bumi)) 
        R1_bc_list.append(R1_bc)
        R2_bc_list.append(R2_bc)
    return aln_dict, alignment_stat, R1_bc_list, R2_bc_list


def bc_rankplot(bc_list,sample,position,qc_pdfs,max_expected_barcodes=100000):

    bc_dict = Counter(bc_list).most_common()
    sub = bc_dict[100:max_expected_barcodes].copy()
    x = np.histogram([np.log10(bc[1]) for bc in sub], 100)
    smooth = gaussian_filter1d(x[0], 3)
    peak_idx,_ = find_peaks(-smooth)
    mean_hist = (x[1][1:][peak_idx]+x[1][:-1][peak_idx])/2
    mean_hist = mean_hist[-1]

    bc_wl = [bc for bc in bc_dict if bc[1]>=10**mean_hist].copy()
    white_list_size=len(bc_wl)
    
    plt.figure(figsize=(4,3))
    log10_ranks=np.log10(np.arange(1,len(bc_dict)+1))
    log10_reads=[np.log10(bc[1]) for bc in bc_dict]
    plt.plot(log10_ranks,log10_reads)#,label='Rank Plot of Reads')
    plt.xlabel('Log10 Ranks')
    plt.ylabel('Log10 Reads')
    plt.title(f'{sample} {position} {white_list_size}')
    plt.plot([0, log10_ranks[-1]], [mean_hist, mean_hist], linewidth=1,label='log10 threshold',c='tab:green')
    log10_wl=np.log10(white_list_size)
    plt.plot([log10_wl, log10_wl], [np.min(log10_reads), np.max(log10_reads)], linewidth=1,label='log10 size',c='tab:orange')
    plt.legend(loc="best");
    
    qc_pdfs.savefig(bbox_inches='tight')
    
    plt.figure(figsize=(4,3))
    plt.plot(x[1][:-1],x[0], label='Raw Histogram')
    plt.plot(x[1][:-1],smooth, label='Gaussian Smoothed')
    plt.xlabel('Log10 UMI Counts')
    plt.ylabel('Bin Height')
    plt.title(f'{sample} {position}')
    plt.plot([mean_hist, mean_hist], [0, np.max(x[0])], linewidth=2,label='Whitelist Threshold')
    plt.legend(loc="best");
    
    qc_pdfs.savefig(bbox_inches='tight')
    return 10**mean_hist


def umi_collapsing(cnt_dict, max_dist=1):
    """
    input: dict of barcode without matching
    output: list of barcode after collapsing
    """
    start_time = time.time()
    clusterer = UMIClusterer(cluster_method="directional")
    clustered_bc = clusterer(cnt_dict, threshold=max_dist)
    clustering_time = time.time()
    cluster_bc = [bc_group[0].decode('utf-8') for bc_group in clustered_bc]
    end_time = time.time()
    print("Clustering time: {}s".format(clustering_time-start_time))
    print("Dict creation time is: {}s".format(end_time-clustering_time))
    print("Total time is: {}s".format(end_time-start_time))
    return cluster_bc


def bc_collapsing(aln_dict, R1_bc_list, R2_bc_list, min_reads_R1, min_reads_R2, alignment_stat):
    """ 
    input: dict of barcode without matching
    output: dict of barcode after filtering and collapsing
    """
    # filter for reads and collapse to whitelist
    R1_list = [s.encode('utf-8') for s in R1_bc_list]
    R1_dict = dict(Counter(R1_list))
    R1_dict_top = {k: v for k, v in R1_dict.items() if v > min_reads_R1}
    R1_whitelist = umi_collapsing(R1_dict_top)
    print("R1 total {}, after filter {}, whitelist {}".format(len(R1_dict),len(R1_dict_top),len(R1_whitelist)))
    print("read percentage: {}".format(np.sum(list(R1_dict_top.values()))/np.sum(list(R1_dict.values()))))
    R2_list = [s.encode('utf-8') for s in R2_bc_list]
    R2_dict = dict(Counter(R2_list))
    R2_dict_top = {k: v for k, v in R2_dict.items() if v > min_reads_R2}
    R2_whitelist = umi_collapsing(R2_dict_top)
    print("R2 total {}, after filter {}, whitelist {}".format(len(R2_dict),len(R2_dict_top),len(R2_whitelist)))
    print("read percentage: {}".format(np.sum(list(R2_dict_top.values()))/np.sum(list(R2_dict.values()))))

    # match to whitelist
    R1_bc_matching_dict,_,_ = barcode_matching(Counter(R1_whitelist), list(set(R1_bc_list)), max_dist=1)
    R2_bc_matching_dict,_,_ = barcode_matching(Counter(R2_whitelist), list(set(R2_bc_list)), max_dist=1)

    # generate dict with matched bc
    aln_dict_new = {}
    for bc_R1 in aln_dict:
        if bc_R1 in R1_bc_matching_dict:
            for R2 in range(len(aln_dict[bc_R1])):
                bc_R2 = aln_dict[bc_R1][R2][0]
                if bc_R2 in R2_bc_matching_dict:
                    alignment_stat["after_filter_reads"] += 1
                    aln_dict_new.setdefault(R1_bc_matching_dict[bc_R1],[]).append(
                        (R2_bc_matching_dict[bc_R2],aln_dict[bc_R1][R2][1],aln_dict[bc_R1][R2][2])) 
    print(len(aln_dict_new))
                    
    return aln_dict_new, alignment_stat


def write_blind(aln_dict_new, alignment_stat, sample, out_dir):
    # collapse for reads
    for bc in aln_dict_new:
        tmp = Counter(aln_dict_new[bc])
        aln_dict_new[bc] = tmp

    # write result to csv
    raw_f = gzip.open(os.path.join(out_dir,(sample+"_blind_raw_reads_filtered.csv.gz")),"wb")
    raw_f.write(b'R1_bc,R2_bc,R1_bumi,R2_bumi,reads\n')
    for bc_R1 in aln_dict_new:
        raw_f.write(bytes('\n'.join(['{},{},{},{},{}'.format(
            bc_R1, it[0], it[1], it[2], aln_dict_new[bc_R1][it]) for it in aln_dict_new[bc_R1]])+'\n',"UTF-8"))
    raw_f.close()
    print("Write matched data to {}".format("blind_raw_reads_filtered.csv.gz"))

    with open(os.path.join(out_dir,(sample+"_blind_statistics_filtered.csv")),"w") as f:
        f.write("alignment_status,counts\n")
        for aln_stat in alignment_stat:
            f.write("{},{}\n".format(aln_stat, alignment_stat[aln_stat]) )


def get_args():
    parser = argparse.ArgumentParser(description='Process recon seq data.')
    parser.add_argument("-d", "--date",
        help="input experiment data.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-s", "--sample",
        help="input sample id.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-r2", "--read2type",
        help="input bead type of read2.",
        type=str,
        default='V9',
    )
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = get_args()

    fq_dir = os.path.join("/broad/thechenlab/Chenlei/spotmapping/fiducial/data",args.date,"fastq")
    # fq_dir = os.path.join("/Volumes/broad_thechenlab/Chenlei/spotmapping/fiducial/data",args.date,"fastq")
    # fq_dir = os.path.join("/mnt/thechenlab/Chenlei/spotmapping/fiducial/data",args.date,"fastq")
    fq_files = [f for f in os.listdir(fq_dir) if not f.startswith('.')]
    fq1_file = os.path.join(fq_dir,[it for it in fq_files if args.sample in it and "R1" in it][0])
    fq2_file = os.path.join(fq_dir,[it for it in fq_files if args.sample in it and "R2" in it][0])

    # make output dir'
    out_dir = os.path.join("/broad/thechenlab/Chenlei/spotmapping/fiducial/data",args.date,args.sample)
    # out_dir = os.path.join("/Volumes/broad_thechenlab/Chenlei/spotmapping/fiducial/data",args.date,args.sample)
    # out_dir = os.path.join("/mnt/thechenlab/Chenlei/spotmapping/fiducial/data",args.date,args.sample)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    # bacode blind
    aln_dict, stat, R1_bc_list, R2_bc_list = barcode_extract(fq1_file, fq2_file, args.read2type)
    print("barcode extracted")

    qc_pdf_file = os.path.join(out_dir, args.sample+'_QC.pdf')
    qc_pdfs = PdfPages(qc_pdf_file)
    R1_threshold = bc_rankplot(R1_bc_list,args.sample,'R1',qc_pdfs,max_expected_barcodes=100000)
    R2_threshold = bc_rankplot(R2_bc_list,args.sample,'R2',qc_pdfs,max_expected_barcodes=100000)
    qc_pdfs.close()

    # R1_threshold = 50
    # R2_threshold = 50

    aln_dict_new, stat_new = bc_collapsing(aln_dict, R1_bc_list, R2_bc_list, min_reads_R1=R1_threshold, min_reads_R2=R2_threshold, alignment_stat = stat)
    write_blind(aln_dict_new, stat_new, args.sample, out_dir)
