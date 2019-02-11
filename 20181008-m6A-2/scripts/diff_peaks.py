'''
read in peaks parsed by `parse_peaks.py` and
perform a variety of ad-hoc filtering 
ZZJ
2.8.2019
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


def read_deseq_file(fn, fdr_cutoff=0.01, logfc_cutoff=2):
	df = pd.read_table(fn)
	log2fc_cutoff = np.log2(logfc_cutoff)
	df = df.dropna()
	sig_genes = df[ (df.padj<fdr_cutoff) & ( (df.log2FoldChange>log2fc_cutoff) | (df.log2FoldChange < -log2fc_cutoff))].index
	return set(sig_genes)

def read_project_deseq(_project_dir):
	project_dir = os.path.join(_project_dir, 'deseq')
	project_dict = {}
	deseq_fn_list = [os.path.join(project_dir, x) for x in os.listdir(project_dir) if x.endswith('.txt')]
	for deseq_fn in deseq_fn_list:
		comparison_name = deseq_fn.split('/')[-1].split('.')[0]
		project_dict[comparison_name] = read_deseq_file(deseq_fn)
	return project_dict



def RPKM_test_new(ratio, input_readcounts, peak_to_gene, gene_dict, peak_df, sig_diff_gene_dict):
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
		'Heart_Age': [
			['NH_1','NH_2','NH_3','NH_4','NH_5','NH_6','NH_7','NH_8','NH_9','NH_10','NH_11','NH_12',],
			['LV_01', 'LV_02', 'LV_03', 'LV_04', 'LV_05', 'LV_06', 'LV_07', 'LV_08','LV_09','LV_10', 'RV_01', 'RV_02', 'RV_03', 'RV_04', 'RV_05', 'RV_06', 'RV_07', 'RV_08','RV_09','RV_10',]
		],
		'Kidney_Age': [
			['P1_NKI_1', 'P1_NKI_2', 'P1_NKI_3', 'P1_NKI_4','P1_NKI_5','P1_NKI_6','P1_NKI_7','P1_NKI_8','P1_NKI_9','P1_NKI_10','P1_NKI_11','P1_NKI_12',],
			['M_KI_01', 'M_KI_02', 'M_KI_03', 'M_KI_04', 'M_KI_05','F_KI_06', 'F_KI_07', 'F_KI_08', 'F_KI_09', 'F_KI_10',]
		],

		'AHK': [
			['M_KI_01', 'M_KI_02', 'M_KI_03', 'M_KI_04', 'M_KI_05','F_KI_06', 'F_KI_07', 'F_KI_08', 'F_KI_09', 'F_KI_10',],
			['LV_01', 'LV_02', 'LV_03', 'LV_04', 'LV_05', 'LV_06', 'LV_07', 'LV_08','LV_09','LV_10', 'RV_01', 'RV_02', 'RV_03', 'RV_04', 'RV_05', 'RV_06', 'RV_07', 'RV_08','RV_09','RV_10',]
		],
		'NHK': [
			['P1_NKI_1', 'P1_NKI_2', 'P1_NKI_3', 'P1_NKI_4','P1_NKI_5','P1_NKI_6','P1_NKI_7','P1_NKI_8','P1_NKI_9','P1_NKI_10','P1_NKI_11','P1_NKI_12',],
			['NH_1','NH_2','NH_3','NH_4','NH_5','NH_6','NH_7','NH_8','NH_9','NH_10','NH_11','NH_12',],
		],

	}

	res = pd.DataFrame(0., index=ratio.index, columns=list(group_comparisons.keys()))
	res_specific = pd.DataFrame(0., index=ratio.index, columns=list(group_comparisons.keys()))
	for comparison in group_comparisons:
		print(comparison)
		this_res = pd.DataFrame(columns=['peak', 'group1_peakIntensity', 'group2_peakIntensity', 'log1.5_fc', 't_pval', 'wilcox_pval'])
		sig_num = 0
		for peak in tqdm(ratio.index):
			# unpack comparison index
			a = group_comparisons[comparison][0]
			b = group_comparisons[comparison][1]
			# peak intensity
			a_ratio = ratio.loc[peak, a]
			b_ratio = ratio.loc[peak, b]

			##--- DEPRECATED ---##
			## because requiring num. of peaks already achieve this	
			## 4) The maximum peak intensity of all samples >= 2; 
			#if not (max(np.concatenate([a_ratio.values, b_ratio.values]))>=2):
			#	continue
			##--- DONE DEPRECATED ---##


			# fold change between groups
			log_a_fc = np.log(a_ratio) - np.mean(np.log(b_ratio))
			log_b_fc = np.log(b_ratio) - np.mean(np.log(a_ratio))
			abs_avg_fc = max(np.mean(a_ratio)/np.mean(b_ratio), np.mean(b_ratio)/np.mean(a_ratio))
			avg_fc = np.exp(np.log(np.mean(a_ratio)) - np.log(np.mean(b_ratio)))

			## 	3) At least 1.5 fold (or 2 fold) change of peak 
			##	intensities in both replicates in the same direction; 
			##--- DEPRECATED ---##
			## let t-test and wilcox-test take care of the distributional things
			#if min(log_a_fc) < 0 < max(log_a_fc) or min(log_b_fc) < 0 < max(log_b_fc):  # check if all the same sign
			#	continue												# by whether list stradles zero
			##--- DONE DEPRECATED ---##

			if not abs_avg_fc > 1.5:
			#if not abs_avg_fc > 2:
				continue

			## 2) Input window RPKM >= 10 in all 4 samples; 
			a_input_win = input_readcounts.loc[peak, a]
			b_input_win = input_readcounts.loc[peak, b]
			#if not ( all(a_input_win>=10) and all(b_input_win>=10) ):
			if np.mean(a_input_win)<=1 or np.mean(b_input_win)<=1:
				continue

			# gene expression
			target_gene = peak_to_gene[peak]
			a_geneexp = np.array([ gene_dict[sam+'_Input'][target_gene] for sam in a ])
			b_geneexp = np.array([ gene_dict[sam+'_Input'][target_gene] for sam in b ])
			
			## 1) Input gene TPM >= 1 in all samples;
			#if all(a_geneexp<=1) and all(b_geneexp<=1):
			if np.mean(a_geneexp)<=1 and np.mean(b_geneexp)<=1:
				continue

			if target_gene in sig_diff_gene_dict[comparison]:
				is_tissue_specific = True
			else:
				is_tissue_specific = False
			
			## 	5) In each replicate, the sample with higher peak 
			## intensity must be called as having peak in half of the replicates.			
			if avg_fc > 1. and sum(np.isnan(peak_df.loc[peak, a])) > len(a) - len(a)*0.5:
				continue
			if avg_fc < 1. and sum(np.isnan(peak_df.loc[peak, b])) > len(b) - len(b)*0.5:
				continue

			## 6) wilcoxon & t-test must be significant
			try:
				pv1 = ss.mannwhitneyu(a_ratio, b_ratio).pvalue
			except:
				pv1 = 1.
			if pv1>0.05:
				continue
			try:
				pv2 = ss.ttest_ind(a_ratio, b_ratio).pvalue
			except:
				pv2 = 1.
			if pv2>0.05:
				continue

			foldchange = math.log(np.mean(a_ratio)/np.mean(b_ratio), 1.5)
			chrom, start, end = peak.split(':')
			new_peak = '{}:{}-{}'.format(chrom, start, end)
			this_res = this_res.append(
				{
					'peak': new_peak, 
					'target_gene': target_gene,
					'is_tissue_specific': is_tissue_specific,
					'group1_avgPeakInts': np.mean(a_ratio),
					'group2_avgPeakInts': np.mean(b_ratio),
					'group1_peakIntensity':','.join([str(round(x,2)) for x in a_ratio]),
					'group2_peakIntensity':','.join([str(round(x,2)) for x in b_ratio]),
					'log1.5_fc': foldchange,
					't_pval': pv2,
					'wilcox_pval': pv1,
					#'group1_inputWinRPKM': ','.join([str(round(x,2)) for x in a_input_win]),
					#'group2_inputWinRPKM': ','.join([str(round(x,2)) for x in b_input_win]),
					'group1_inputWinRPKM': np.mean(a_input_win),
					'group2_inputWinRPKM': np.mean(b_input_win),
					#'group1_geneExp': ','.join([str(round(x,2)) for x in a_geneexp]), 
					#'group2_geneExp': ','.join([str(round(x,2)) for x in b_geneexp]), 
					'group1_geneExp': np.mean(a_geneexp), 
					'group2_geneExp': np.mean(b_geneexp), 
					'group1_numPeaks': sum(~np.isnan(peak_df.loc[peak, a])),
					'group2_numPeaks': sum(~np.isnan(peak_df.loc[peak, b])),

				}, 
				ignore_index=True)
			sig_num += 1
			if is_tissue_specific:
				res_specific.loc[peak, comparison] = foldchange
			else:
				res.loc[peak, comparison] = foldchange
		this_res.to_csv("_data/DiffPeaks_tsv/{}.DiffPeak.tsv".format(comparison), sep="\t")
	res = reindex_res_df(res)
	res_specific = reindex_res_df(res_specific)
	return res, res_specific


def reindex_res_df(res):
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


def run_and_save():
	PROJECT_DIR = '../_data/m6A/'
	T2G_FILE = '../_data/t2g/t2g_mm10.txt'

	# read in kallisto gene level estimates
	print('read t2g')
	t2g = read_t2g(T2G_FILE)
	print('read gene quant')
	gene_dict = read_project_kallisto(PROJECT_DIR, t2g)
	sig_diff_gene_dict = read_project_deseq(PROJECT_DIR)

	# read in clam peaks
	print('read peak')
	df_with_sig, peak_to_gene, peak_df = read_project_clam(PROJECT_DIR)

	# load previously counts
	input_readcounts = pd.read_table('../_data/m6A/parsed_peaks/input_peak.RPKM.csv', index_col=0)
	ip_readcounts = pd.read_table('../_data/m6A/parsed_peaks/ip_peak.RPKM.csv', index_col=0)

	# compute the peak intensity
	ratio = pd.read_table('../_data/m6A/parsed_peaks/peak_intensity.csv', index_col=0)

	# test for differential
	print('compute test')
	differential_sites, differential_sites_specific = RPKM_test_new(ratio, input_readcounts, peak_to_gene, gene_dict, peak_df, sig_diff_gene_dict)
	differential_sites.to_csv('../_data/m6A/parsed_peaks/differential_sites_log1.5.csv', sep='\t')
	differential_sites_specific.to_csv('../_data/m6A/parsed_peaks/differential_sites_log1.5-tissue_specific.csv', sep='\t')


if __name__ == '__main__':
	run_and_save()
