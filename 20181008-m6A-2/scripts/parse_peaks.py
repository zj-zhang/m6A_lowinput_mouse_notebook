'''
read in peaks and compute t-statistic
ZZJ
8.27.2018
'''

import os
import pandas as pd
import numpy as np
import scipy.stats as ss
import pysam
from collections import defaultdict
import math
from tqdm import tqdm

def read_t2g(fn='/u/nobackup/yxing/NOBACKUP/frankwoe/CLAM_m6A_Snakemake/t2g_mm10.txt'):
	tmp = pd.read_table(fn, header=0)
	t2g = {}
	for i in range(tmp.shape[0]):
		t2g[tmp.iloc[i,0]] = tmp.iloc[i,1]
	return t2g

def read_tx_file(fn, t2g):
	tx = pd.read_table(fn, header=0)
	gene_dict = defaultdict(float)
	unmapped = 0
	for i in range(tx.shape[0]):
		t = tx.loc[i, 'target_id']
		t = t.split('.')[0]
		if t in t2g:
			g = t2g[t]
		else:
			unmapped += 1
			continue
		gene_dict[g] += tx.loc[i, 'tpm']
	#if unmapped:
	#	print('%i tx are not mappable in this t2g index.'%unmapped)
	return gene_dict


def read_project_kallisto(_project_dir, t2g):
	project_dir = os.path.join(_project_dir, 'kallisto')
	tx_fn_list = [os.path.join(project_dir, x, 'abundance.tsv') for x in os.listdir(project_dir)]
	sample_gene_dict = {}
	for tx_fn in tx_fn_list:
		gene_dict = read_tx_file(tx_fn, t2g)
		sample_name = tx_fn.split('/')[-2]
		sample_gene_dict[sample_name] = gene_dict
	return sample_gene_dict


def read_peak_file(fn, peak_to_gene=None):
	peak_dict = {}
	if peak_to_gene is None:
		peak_to_gene = {}
	with open(fn, 'r') as f:
		for line in f:
			ele = line.strip().split()
			peak = ':'.join(ele[0:3])
			score = float(ele[-4])
			gene = ele[3].split('-')[0]
			peak_dict[peak] = score
			peak_to_gene[peak] =  gene
	return peak_dict, peak_to_gene


def read_project_clam(_project_dir):
	project_dir = os.path.join(_project_dir, 'clam')
	project_dict = {}
	peak_fn_list = [os.path.join(project_dir, x, 'narrow_peak.unique.bed.all') for x in os.listdir(project_dir) if x.startswith('peaks-')]
	for peak_fn in peak_fn_list:
		peak_name = peak_fn.split('/')[-2].split('-')[1].rstrip('_IP')
		project_dict[peak_name] = read_peak_file(peak_fn)[0]

	df = pd.DataFrame.from_dict(project_dict)

	sig_peaks = set()
	peak_to_gene = dict()
	peak_dict = dict()
	peak_fn_list2 = [os.path.join(project_dir, x, 'narrow_peak.unique.bed') for x in os.listdir(project_dir) if x.startswith('peaks-')]
	for peak_fn in peak_fn_list2:
		peak_name = peak_fn.split('/')[-2].split('-')[1].rstrip('_IP')
		tmp, peak_to_gene = read_peak_file(peak_fn, peak_to_gene)
		_ = map(sig_peaks.add, tmp.keys())
		peak_dict[peak_name] = tmp

	df_with_sig = df.loc[sig_peaks]
	peak_df = pd.DataFrame.from_dict(peak_dict)
	return df_with_sig, peak_to_gene, peak_df


def wilcoxon_test(df_with_sig):
	group_comparisons = {
		'KI_gender': [
			['F_KI_06', 'F_KI_07', 'F_KI_08', 'F_KI_09', 'F_KI_10',],
			['M_KI_01', 'M_KI_02', 'M_KI_03', 'M_KI_04', 'M_KI_5',]
		],
		'V_position': [
			['LV_01', 'LV_02', 'LV_03', 'LV_04', 'LV_05', 'LV_06', 'LV_07', 'LV_08','LV_09','LV_10',],
			['RV_01', 'RV_02', 'RV_03', 'RV_04', 'RV_05', 'RV_06', 'RV_07', 'RV_08','RV_09','RV_10',],
		],
	}

	res = pd.DataFrame(1, index=df_with_sig.index, columns=list(group_comparisons.keys()))
	for peak in df_with_sig.index:
		for comparison in group_comparisons:
			a = df_with_sig.loc[peak, group_comparisons[comparison][0]]
			b = df_with_sig.loc[peak, group_comparisons[comparison][1]]
			a = a[~np.isnan(a)]
			b = b[~np.isnan(b)]
			fc = np.exp(abs(np.mean(a) - np.mean(b)))
			if fc>2:
				pv = ss.mannwhitneyu(a,b).pvalue
				res.loc[peak, comparison] = pv


def RPKM_test_new(ratio, input_readcounts, peak_to_gene, gene_dict, peak_df):
	'''
	1) Input gene FPKM >= 1 in all 4 samples; 
	2) Input window RPKM >= 10 in all 4 samples; 
	3) At least 1.5 fold (or 2 fold) change of peak 
		intensities in both replicates in the same direction; 
	4) The maximum peak intensity of all samples >= 2; 
	5) In each replicate, the sample with higher peak 
		intensity must be called as having peak.
	'''
	group_comparisons = {
		'KI_age': [
			['P1_NKI_1', 'P1_NKI_2', 'P1_NKI_3', 'P1_NKI_4','P1_NKI_5','P1_NKI_6','P1_NKI_7','P1_NKI_8','P1_NKI_9','P1_NKI_10','P1_NKI_11','P1_NKI_12',],
			['M_KI_01', 'M_KI_02', 'M_KI_03', 'M_KI_04', 'M_KI_05','F_KI_06', 'F_KI_07', 'F_KI_08', 'F_KI_09', 'F_KI_10',]
		],
		'V_age': [
			['NH_1','NH_2','NH_3','NH_4','NH_5','NH_6','NH_7','NH_8','NH_9','NH_10','NH_11','NH_12',],
			['LV_01', 'LV_02', 'LV_03', 'LV_04', 'LV_05', 'LV_06', 'LV_07', 'LV_08','LV_09','LV_10', 'RV_01', 'RV_02', 'RV_03', 'RV_04', 'RV_05', 'RV_06', 'RV_07', 'RV_08','RV_09','RV_10',]
		],
		'adult_KI_V': [
			['M_KI_01', 'M_KI_02', 'M_KI_03', 'M_KI_04', 'M_KI_05','F_KI_06', 'F_KI_07', 'F_KI_08', 'F_KI_09', 'F_KI_10',],
			['LV_01', 'LV_02', 'LV_03', 'LV_04', 'LV_05', 'LV_06', 'LV_07', 'LV_08','LV_09','LV_10', 'RV_01', 'RV_02', 'RV_03', 'RV_04', 'RV_05', 'RV_06', 'RV_07', 'RV_08','RV_09','RV_10',]
		],
		'neonatal_KI_V': [
			['P1_NKI_1', 'P1_NKI_2', 'P1_NKI_3', 'P1_NKI_4','P1_NKI_5','P1_NKI_6','P1_NKI_7','P1_NKI_8','P1_NKI_9','P1_NKI_10','P1_NKI_11','P1_NKI_12',],
			['NH_1','NH_2','NH_3','NH_4','NH_5','NH_6','NH_7','NH_8','NH_9','NH_10','NH_11','NH_12',],
		],

	}

	res = pd.DataFrame(0., index=ratio.index, columns=list(group_comparisons.keys()))
	for comparison in group_comparisons:
		print(comparison)
		this_res = pd.DataFrame(columns=['peak', 'group1', 'group2', 'log1.5_fc', 't_pval', 'wilcox_pval'])
		for peak in tqdm(ratio.index):
			# unpack comparison index
			a = group_comparisons[comparison][0]
			b = group_comparisons[comparison][1]
			# peak intensity
			a_ratio = ratio.loc[peak, a]
			b_ratio = ratio.loc[peak, b]
			## 4) The maximum peak intensity of all samples >= 2; 
			if not (max(np.concatenate([a_ratio.values, b_ratio.values]))>=2):
				continue

			# fold change between groups
			log_a_fc = np.log(a_ratio) - np.mean(np.log(b_ratio))
			log_b_fc = np.log(b_ratio) - np.mean(np.log(a_ratio))
			#abs_avg_fc = np.exp(abs(np.mean(np.log(a_ratio)) - np.mean(np.log(b_ratio))))
			abs_avg_fc = max(np.mean(a_ratio)/np.mean(b_ratio), np.mean(b_ratio)/np.mean(a_ratio))
			avg_fc = np.exp(np.mean(np.log(a_ratio)) - np.mean(np.log(b_ratio)))
			## 	3) At least 1.5 fold (or 2 fold) change of peak 
			##	intensities in both replicates in the same direction; 
			#if min(log_a_fc) < 0 < max(log_a_fc) or min(log_b_fc) < 0 < max(log_b_fc):  # check if all the same sign
			#	continue												# by whether list stradles zero
			if not abs_avg_fc > 1.5:
				continue

			## 2) Input window RPKM >= 10 in all 4 samples; 
			a_input_win = input_readcounts.loc[peak, a]
			b_input_win = input_readcounts.loc[peak, b]
			if not ( all(a_input_win>=10) and all(b_input_win>=10) ):
				continue

			# gene expression
			target_gene = peak_to_gene[peak]
			a_geneexp = np.array([ gene_dict[sam+'_Input'][target_gene] for sam in a ])
			b_geneexp = np.array([ gene_dict[sam+'_Input'][target_gene] for sam in b ])
			## 1) Input gene FPKM >= 1 in all 4 samples;
			if not all(a_geneexp)>=1 and not all(b_geneexp)>=1:
				continue
			## 	5) In each replicate, the sample with higher peak 
			## intensity must be called as having peak.
			#if all(np.isnan(peak_df.loc[peak, a])) and all(np.isnan(peak_df.loc[peak, b])):
			if avg_fc > 1. and sum(np.isnan(peak_df.loc[peak, a])) > len(a) - 2:
				continue
			if avg_fc < 1. and sum(np.isnan(peak_df.loc[peak, b])) > len(b) - 2:
				continue

			## 6) wilcoxon & t-test must be significant
			pv1 = ss.mannwhitneyu(a_ratio, b_ratio).pvalue
			if pv1>0.05:
				continue
			pv2 = ss.ttest_ind(a_ratio, b_ratio).pvalue
			if pv2>0.05:
				continue
			foldchange = math.log(np.mean(a_ratio)/np.mean(b_ratio), 1.5)
			chrom, start, end = peak.split(':')
			new_peak = '{}:{}-{}'.format(chrom, start, end)
			this_res = this_res.append(
				{
					'peak': new_peak, 
					'group1':','.join([str(round(x,2)) for x in a_ratio]),
					'group2':','.join([str(round(x,2)) for x in b_ratio]),
					'log1.5_fc': foldchange,
					't_pval': pv2,
					'wilcox_pval': pv1 
				}, 
				ignore_index=True)
			res.loc[peak, comparison] = foldchange
		this_res.to_csv("{}.DiffPeak.tsv".format(comparison), sep="\t")
	new_index=[]
	url=[]
	url_template='=HYPERLINK("https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&position={0}","{0}")'
	for peak in res.index:
		chrom, start, end = peak.split(':')
		new_peak = '{}:{}-{}'.format(chrom, start, end)
		new_index.append(new_peak)
		url.append(url_template.format(new_peak))
	res.index = new_index
	res['url'] = pd.Series(url, index=new_index)
	return res


def get_total_reads(bam_filename):
	idxstats  = pysam.idxstats(bam_filename).split('\n')
	tot = 0
	for l in idxstats:
		if not l:
			continue
		ele = l.split('\t')
		tot += int(ele[-2])
	return tot


def get_reads_window(bam, chrom, start, end):
	reads  = [ 1 for x in bam.fetch(chrom ,start, end) ]
	return len(reads)


def make_file_handlers(_project_dir):
	project_dir = os.path.join(_project_dir, 'clam')
	bam_dict = {}
	bam_fn_list = [os.path.join(project_dir, x, 'unique.sorted.collapsed.bam') for x in os.listdir(project_dir) \
		if not x.startswith('peaks-') and os.path.isdir(project_dir+'/'+x)]
	#print(bam_fn_list)
	for bam_fn in bam_fn_list:
		sample_name = bam_fn.split('/')[-2]
		bam_dict[sample_name] = pysam.Samfile(bam_fn, 'rb')
	return bam_dict


def count(df_with_sig, bam_dict, project_dir):
	'''
	Given a list of peaks and samples, count the RPKM,
	then save to file.
	'''
	input_readcounts = pd.DataFrame(0, index=df_with_sig.index, columns=df_with_sig.columns)
	ip_readcounts = pd.DataFrame(0, index=df_with_sig.index, columns=df_with_sig.columns)
	
	input_total_counts = {
		x+'_Input' : get_total_reads(project_dir+'/clam/'+x+'_Input/unique.sorted.bam')/float(10**6)
		 for x in df_with_sig.columns
		 }
	ip_total_counts = {
		x+'_IP':get_total_reads(project_dir+'/clam/'+x+'_IP/unique.sorted.bam')/float(10**6)
	 	for x in df_with_sig.columns
	 	}

	for peak in df_with_sig.index:
		chrom, start, end = peak.split(':')
		start = int(start)
		end = int(end)
		for sam in df_with_sig.columns:
			input_name = sam+'_Input'
			this_bam = bam_dict[input_name]
			this_input = get_reads_window(this_bam, chrom, start, end)
			input_readcounts.loc[peak, sam] = \
				this_input / float(input_total_counts[input_name]) / ((end-start)/1000.)

			ip_name = sam+'_IP'
			this_bam = bam_dict[ip_name]
			this_ip = get_reads_window(this_bam, chrom, start, end)
			ip_readcounts.loc[peak, sam] = \
				this_ip / float(ip_total_counts[ip_name]) / ((end-start)/1000.)
	return input_readcounts, ip_readcounts


def get_peak_intensity(ip_readcounts, input_readcounts):
	pseudo_count = 1.
	ratio = pd.DataFrame(0, index=ip_readcounts.index, columns=ip_readcounts.columns)
	for peak in ratio.index:
		input_rpkm = input_readcounts.loc[peak]
		ip_rpkm = ip_readcounts.loc[peak]
		ratio.loc[peak] = (ip_rpkm + pseudo_count) / (input_rpkm + pseudo_count)
	return ratio


def run_and_save():
	PROJECT_DIR = '../_data/m6A/'
	T2G_FILE = '../_data/t2g/t2g_mm10.txt'

	# read in kallisto gene level estimates
	print('read t2g')
	t2g = read_t2g(T2G_FILE)
	print('read gene quant')
	gene_dict = read_project_kallisto(PROJECT_DIR, t2g)

	# read in clam peaks
	print('read peak')
	df_with_sig, peak_to_gene, peak_df = read_project_clam(PROJECT_DIR)

	# get peak intensities
	#bam_dict = make_file_handlers(PROJECT_DIR)
	#input_readcounts, ip_readcounts = count(df_with_sig, bam_dict, PROJECT_DIR)

	# save the RPKM values
	#input_readcounts.to_csv('input_peak.RPKM.csv', sep='\t')
	#ip_readcounts.to_csv('ip_peak.RPKM.csv', sep='\t')

	# load previously counts
	input_readcounts = pd.read_table('../_data/m6A/parsed_peaks/input_peak.RPKM.csv', index_col=0)
	ip_readcounts = pd.read_table('../_data/m6A/parsed_peaks/ip_peak.RPKM.csv', index_col=0)

	# compute the peak intensity
	#ratio = get_peak_intensity(ip_readcounts, input_readcounts)
	#ratio.to_csv('peak_intensity.csv', sep='\t')

	ratio = pd.read_table('../_data/m6A/parsed_peaks/peak_intensity.csv', index_col=0)

	# test for differential
	print('compute test')
	differential_sites = RPKM_test_new(ratio, input_readcounts, peak_to_gene, gene_dict, peak_df)
	differential_sites.to_csv('../_data/m6A/parsed_peaks/differential_sites_log1.5.csv', sep='\t')


if __name__ == '__main__':
	run_and_save()
