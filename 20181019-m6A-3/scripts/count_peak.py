

def intersect_peaks(replicate_list, peak_dir):
    '''
    function:: intersect_peaks(replicate_list, peak_dir)
    Intersect peaks for a group of replicates in the same condition.
    :param replicate_list: a list of filepaths replicates
    :param peak_dir: str, directory to all CLAM peaks, usually the 
      "project/clam" folder in `CLAM_eCLIP_Snakemake`
    :rtype: a Pandas DataFrame, with each row being peaks
    .. example::
        An example peak(row) will be like this:
        gene                   ENSMUSG00000019966
        strand                                  +
        replicate         F_KI_10,F_KI_09,F_KI_06
        ratio                      2.75,2.73,2.80
        num_replicates                          3
        Name: chr10:100096030-100096130:+, dtype: object
    '''
    import os
    import pandas as pd
    from collections import defaultdict
    this_peaks = defaultdict(dict)
    # read through to collect info from CLAM output
    for replicate in replicate_list:
        rep_name = '_'.join(replicate.split('-')[1].split('_')[:-1])
        fn = os.path.join(peak_dir, replicate, 'narrow_peak.unique.bed')
        with open(fn, 'r') as f:
            for line in f:
                ele = line.strip().split()
                chrom = ele[0]
                start = ele[1]
                end = ele[2]
                gene = ele[3].split('-')[0]
                strand = ele[5]
                ratio = ele[6]
                peak_id = '{}:{}-{}:{}'.format(chrom, start, end, strand)
                this_peaks[peak_id]['gene'] = gene
                this_peaks[peak_id]['strand'] = strand
                if 'replicate' in this_peaks[peak_id]:                    
                    this_peaks[peak_id]['replicate'].append(rep_name)
                    this_peaks[peak_id]['ratio'].append(ratio)
                else:
                    this_peaks[peak_id]['replicate'] = [rep_name]
                    this_peaks[peak_id]['ratio'] = [ratio]
    # clean up the list objects and add counts
    for peak_id in this_peaks:
        this_peaks[peak_id]['num_replicates'] = len(set(this_peaks[peak_id]['replicate']))
        this_peaks[peak_id]['replicate'] = ','.join(this_peaks[peak_id]['replicate'])
        this_peaks[peak_id]['ratio'] = ','.join(this_peaks[peak_id]['ratio'])
    # convert to pd.DataFrame
    this_peaks_df = pd.DataFrame.from_dict(this_peaks, orient="index")
    return this_peaks_df