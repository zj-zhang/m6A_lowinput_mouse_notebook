'''
Assign peaks to genomic regions
Zijun Zhang
8.1.2018
10.25.2018: wrapped to a function with document
'''

def intersect_gtf_regions(peak_fp, outfn, gtf_dir, do_plot=False):
	'''function: intersect_gtf_regions(peak_fp, outfn, gtf_dir, [do_plot=False])
	Intersect a peak BED file with a list of genomic region annotations (e.g. start/stop codon, UTR, intron), 
	output the intersection counts and peak-region annotations.
	:param peak_fp: filepath to a BED-format peak
	:param outfn: filepath to output count file, has to end with ".txt"; annotation will be "NNN.annot.txt"
	:param do_plot: plot the bar-chart, optional
	:rtype: pd.DataFrame, the count for each region annotation.
	.. example::
			category	count	norm_count
			0	start	593	0.246556
			1	stop	1606	0.720184
			2	3UTR	6737	0.006684
			3	5UTR	2375	0.011469
	'''

	import sys
	import os
	import pybedtools
	import matplotlib
	matplotlib.use('agg')
	import matplotlib.pyplot as plt
	import pandas as pd
	from collections import defaultdict


	### input arguments
	#peak_fp, genome, outfn = sys.argv[1], sys.argv[2], sys.argv[3]
	#par_dir = sys.argv[4]


	### make pybedtools objects
	peaks = pybedtools.BedTool(peak_fp)
	ref_dict = {
	    'start': pybedtools.BedTool(os.path.join(gtf_dir, 'start.bed')),
	    'stop': pybedtools.BedTool(os.path.join(gtf_dir, 'stop.bed')),
		'exon': pybedtools.BedTool(os.path.join(gtf_dir, 'exons.bed')),
		'3UTR': pybedtools.BedTool(os.path.join(gtf_dir, '3UTRs.bed')),
		'5UTR': pybedtools.BedTool(os.path.join(gtf_dir, '5UTRs.bed')),
		'cds': pybedtools.BedTool(os.path.join(gtf_dir, 'cds.bed')),
		'intron': pybedtools.BedTool(os.path.join(gtf_dir, 'introns.bed')),
		'proximal200': pybedtools.BedTool(os.path.join(gtf_dir, 'proximal200_intron.bed')),
		'proximal500': pybedtools.BedTool(os.path.join(gtf_dir, 'proximal500_intron.bed'))
	}

	get_total_length = lambda featureset: sum([x.end-x.start-1 for x in featureset])

	### process reference for use
	target = {
		"start": ref_dict['start'],
		"stop": ref_dict["stop"],
		"3UTR": ref_dict['3UTR'],
		"5UTR": ref_dict['5UTR'],
		"CDS": ref_dict['cds'],
		"other_exon": ref_dict['exon']-ref_dict['3UTR']-ref_dict['5UTR']-ref_dict['cds'],
		"px200_intron": ref_dict['proximal200'],
		"px500_intron": ref_dict['proximal500'].subtract(ref_dict['proximal200']),
		"distal_intron": ref_dict['intron'].subtract(ref_dict['exon']).subtract(ref_dict['proximal500'])
	}

	category_list = ['start', 'stop', '3UTR', '5UTR', 'CDS', 'other_exon', "px200_intron", "px500_intron", "distal_intron"]

	### compute counts and length-normalized counts
	len_dict = {x:get_total_length(target[x]) for x in category_list}

	#count_dict = {x: (peaks+target[x]).count() for x in category_list}
	intersect_handler_dict = {x: (peaks+target[x]) for x in category_list}
	count_dict = {x: intersect_handler_dict[x].count() for x in category_list}


	count_df = pd.DataFrame({
		'category':category_list, 
		'count':[count_dict[x] for x in category_list],
		'norm_count':[count_dict[x]/float(len_dict[x]) for x in category_list],
		})
	count_df.norm_count /= sum(count_df.norm_count)

	count_df.to_csv(outfn, sep='\t')
	

	### plot
	if do_plot:
		fig = plt.figure()
		ax  =fig.add_subplot(111)
		#ax.tick_params(axis='x',rotation=45)
		ax2 = ax.twinx()
		width = 0.25

		count_df.plot(x='category', y='count', kind='bar', color='red', ax=ax, width=width, position=1)
		count_df.plot(x='category', y='norm_count', kind='bar', color='blue', ax=ax2, width=width, position=0)

		ax.set_ylabel('count')
		ax2.set_ylabel('length-normalized count')
		fig.tight_layout()
		plt.xticks(rotation=45)
		fig.savefig(outfn.rstrip('txt')+'png')
		plt.close()

	### output the annotations
	annot_fn = outfn.rstrip('txt')+'annot.txt'
	annot_dict = defaultdict(list)
	for x in intersect_handler_dict:
		for peak in intersect_handler_dict[x]:
			# chrom, start, end, gene, strand
			annot_dict[(peak[0], peak[1], peak[2], peak[3], peak[5])].append(x)
	

	### fine-tune the counts: if start/stop > UTR > CDS > exon > px_intron > distal_intron
	new_count_dict = defaultdict(int)
	for peak in annot_dict:
		if 'start' in annot_dict[peak]:
			annot_dict[peak] = ['start']
			new_count_dict['start'] += 1
		elif 'stop' in annot_dict[peak]:
			annot_dict[peak] = ['stop']
			new_count_dict['stop'] += 1
		elif '3UTR' in annot_dict[peak]:
			annot_dict[peak] = ['3UTR']
			new_count_dict['3UTR'] += 1
		elif '5UTR' in annot_dict[peak]:
			annot_dict[peak] = ['5UTR']
			new_count_dict['5UTR'] += 1
		elif 'CDS' in annot_dict[peak]:
			annot_dict[peak] = ['CDS']
			new_count_dict['CDS'] += 1
		elif 'other_exon' in annot_dict[peak]:
			annot_dict[peak] = ['other_exon']
			new_count_dict['other_exon'] += 1
		elif 'px200_intron' in annot_dict[peak] or 'px500_intron' in annot_dict[peak]:
			annot_dict[peak] = ['px_intron']
			new_count_dict['px_intron'] += 1
		elif 'distal_intron' in annot_dict[peak]:
			annot_dict[peak] = ['distal_intron']
			new_count_dict['distal_intron'] += 1

	new_cat_list = ['start', 'stop', '3UTR', '5UTR', 'CDS', 'other_exon', 'px_intron', 'distal_intron']
	len_dict['px_intron'] = len_dict['px200_intron'] + len_dict['px500_intron']
	new_count_df = pd.DataFrame({
		'category':new_cat_list, 
		'count':[new_count_dict[x] for x in new_cat_list],
		'norm_count':[new_count_dict[x]/float(len_dict[x]) for x in new_cat_list],
		})
	new_count_df.norm_count /= sum(new_count_df.norm_count)

	new_count_df.to_csv(outfn.rstrip('txt')+'fineTune.txt', sep='\t')


	with open(annot_fn, 'w') as fo:
		for peak in annot_dict:
			fo.write("\t".join(peak)+'\t'+','.join(annot_dict[peak])+'\n')

	return new_count_df