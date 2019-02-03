import os
import pandas as pd
import numpy as np
import scipy.stats as ss
from collections import defaultdict


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


def read_t2g_symbol(fn='/u/nobackup/yxing/NOBACKUP/frankwoe/CLAM_m6A_Snakemake/t2g_mm10.txt'):
	tmp = pd.read_table(fn, header=0)
	ens2symbol = {}
	for i in range(tmp.shape[0]):
		ens2symbol[tmp.iloc[i,1]] = tmp.iloc[i,2]
	return ens2symbol


def write_gene_list():
	PROJECT_DIR = '/u/home/f/frankwoe/nobackup/CLAM_m6A_Snakemake/projects/Mouse-all/projects/Mouse-all'
	T2G_FILE = '/u/nobackup/yxing/NOBACKUP/frankwoe/CLAM_m6A_Snakemake/t2g_mm10.txt'

	# read in clam peaks
	df_with_sig, peak_to_gene, peak_df = read_project_clam(PROJECT_DIR)

	# read in kallisto gene level estimates
	ens2symbol = read_t2g_symbol(T2G_FILE)

	# read in diff sites
	diff_sites = pd.read_table('differential_sites_log1.5.csv', index_col=0)

	# write_gene_list
	for compr in diff_sites.columns:
		with open('genes_%s.txt'%compr, 'w') as f:
			ens_genes = [peak_to_gene[x] for x in diff_sites.index[diff_sites[compr]!=0].values]
			gene_symbol = set([ens2symbol[x] for x in ens_genes if x in ens2symbol])
			f.write('\n'.join(gene_symbol))

