#!/usr/bin/env python3
#%%
import pandas as pd
import numpy as np
import csv
import re
import argparse
import gtfparse
import multiprocessing

#%%
REQUIRED_COLUMNS = [
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
]
#%%
def make_attribute_column(row):
    attribute = ''
    for i in row.index:
        if i not in REQUIRED_COLUMNS and row[i] != '':
            tmp = f'{i} "{row[i]}"; '
            attribute = attribute + tmp
    attribute = attribute.strip()
    return attribute
#%%
def set_transcript_ranges(group):
    exons = group.query('feature == "exon"')
    if len(exons) == 0:
        return group
    exon_min = min( min(exons['start']), min(exons['end']))
    exon_max = max( max(exons['start']), max(exons['end']))
    group.loc[group.feature == 'transcript', 'start'] = exon_min
    group.loc[group.feature == 'transcript', 'end'] = exon_max
    return group


def transform_cds(gtf):
    # only keep transcriipt and CDS rows
    gtf_cds = (
        gtf
            .query('feature in ["transcript","CDS"]')
    )
    # rename CDS to exon and update transcript range
    gtf_cds['feature'] = np.where(gtf_cds['feature'] =='CDS','exon','transcript')
    gtf_cds = gtf_cds.groupby('transcript_id').apply(set_transcript_ranges)
    # only keep transcripts that have a CDS
    sizes = gtf_cds.groupby('transcript_id').size()
    transcripts_with_cds = list(sizes[sizes > 1].index)
    gtf_cds = gtf_cds.query(f'transcript_id in {transcripts_with_cds}')
    return gtf_cds   


#%%

def process_gtf(sample, name):
    # transform cds info
    sample_cds = transform_cds(sample)
    # make attribute column for cds df and only keep gtf and attribute columns
    sample_cds['attribute'] = sample_cds.apply(make_attribute_column, axis=1)
    sample_cds = sample_cds.filter(REQUIRED_COLUMNS + ['attribute'])
    sample_cds.to_csv(f'{name}.cds_renamed_exon.gtf', sep='\t',index=False, header=False, quoting=csv.QUOTE_NONE)
    #make attribute column for exon df
    sample['attribute'] = sample.apply(make_attribute_column, axis=1)
    sample = sample.filter(REQUIRED_COLUMNS + ['attribute'])
    sample_exon = (
        sample
            .query('feature in ["transcript","exon"]')
            .copy()  
    )
    sample_exon.to_csv(f'{name}.transcript_exons_only.gtf', sep='\t',index=False, header=False, quoting=csv.QUOTE_NONE)

def process_gtf_single(sample):
    transformed = transform_cds(sample)
    transformed['attribute'] = transformed.apply(make_attribute_column, axis=1)
    transformed = transformed.filter(REQUIRED_COLUMNS + ['attribute'])
    return transformed
    
    
def process_gtf_multiprocess(sample,name, num_cores):
    # transform cds info
    chromosomes = sample['seqname'].unique()
    sample_split = [sample[sample['seqname'] == csome] for csome in chromosomes]
    pool = multiprocessing.Pool(processes = num_cores)
    sample_cds_split = pool.map(process_gtf_single, sample_split)
    sample_cds = pd.concat(sample_cds_split)
    sample_cds.to_csv(f'{name}.cds_renamed_exon.gtf', sep='\t',index=False, header=False, quoting=csv.QUOTE_NONE)
    #make attribute column for exon df
    sample['attribute'] = sample.apply(make_attribute_column, axis=1)
    sample = sample.filter(REQUIRED_COLUMNS + ['attribute'])
    sample_exon = (
        sample
            .query('feature in ["transcript","exon"]')
    )
    sample_exon.to_csv(f'{name}.transcript_exons_only.gtf', sep='\t',index=False, header=False, quoting=csv.QUOTE_NONE)
    

def process_sample_rename(sample_file, name, num_cores):

    sample = gtfparse.read_gtf(sample_file)
    sample['transcript_id'] = sample['transcript_id'].apply(lambda x: x.split('|')[1])
    process_gtf_multiprocess(sample, name, num_cores)

def process_reference_rename(reference_file, name, num_cores):
    
    ref = gtfparse.read_gtf(reference_file)
    process_gtf_multiprocess(ref, name, num_cores)

    

    # gtf_colnames = ['seqname','source','feature','start','end','score','strand','frame','attribute']
    # ref = pd.read_table(reference_file, names=gtf_colnames, skiprows=skiprows, )
    # ref_cds = (
    #     ref
    #         .query('feature in ["transcript","CDS"]')
    #         .copy()  
    # )
    # ref_cds['feature'] = np.where(ref_cds['feature'] =='CDS','exon','transcript')

    # savefile = f'{name}.cds_renamed_exon.gtf'
    
    # ref_cds.to_csv(savefile, sep='\t',index=False, header=False, quoting=csv.QUOTE_NONE)

    # ref_exons = (
    #     ref
    #         .query('feature in ["transcript","exon"]')
    # )
    # savefile = f'{name}.transcript_exon_only.gtf'
    # # with open(reference_file, 'r') as iref, open(savefile, 'w') as out:
    # #     for i in range(5):
    # #         out.write(iref.readline())
    # ref_exons.to_csv(savefile, sep='\t',index=False, header=False, quoting=csv.QUOTE_NONE)


def main():
    parser = argparse.ArgumentParser(description='rename cds to exon for sqanti protein module')
    parser.add_argument('--sample_gtf', action='store', dest= 'sample_gtf',help='sample gtf file')
    parser.add_argument('--sample_name', action='store', dest= 'sample_name',help='sample name')
    parser.add_argument('--reference_gtf', action='store', dest= 'reference_gtf',help='sample gtf file')
    parser.add_argument('--reference_name', action='store', dest= 'reference_name',help='sample name')
    parser.add_argument('--num_cores', action='store',dest='num_cores',help='number of cores to use in multiprocessing',default=8,type=int)
    results = parser.parse_args()
    
    process_sample_rename(results.sample_gtf, results.sample_name, results.num_cores)
    process_reference_rename(results.reference_gtf, results.reference_name, results.num_cores)

#%%
    
if __name__ == "__main__":
    main()

# %%
