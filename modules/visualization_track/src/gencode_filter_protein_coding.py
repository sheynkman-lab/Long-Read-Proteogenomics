#!/usr/bin/env python3
import argparse
import re

def filter_protein_coding_gtf(gencode_gtf_file, isonames_with_cds):
    gene_name_pattern = 'gene_name \"(.*?)\";'
    isoname_pattern = 'transcript_name \"(.*?)\";'
    with open(gencode_gtf_file) as ifile, open("gencode.filtered.gtf", "w") as ofile:
        for line in ifile:
            if line.startswith('#'): continue
            wds = line.split('\t')
            if wds[2] in ('transcript', 'exon', 'CDS'):
                isoname_result = re.search(isoname_pattern, wds[8])
                if isoname_result is not None:
                    isoname = isoname_result.group(1)
                    if isoname in isonames_with_cds:
                        gene = re.search(gene_name_pattern, wds[8]).group(1)
                        if wds[0] in ('chrX','chrY'):
                            gene=f'{gene}.{wds[0]}'
                            isoname=f'{isoname}.{wds[0]}'
                        odata = wds[0:8]
                        odata[1] = 'gencode'
                        odata.append(f'gene_id "{gene}"; transcript_id "{isoname}";')
                        ofile.write('\t'.join(odata) + "\n")

def get_transcripts_with_cds(gencode_gtf_file):
    isoname_pattern = 'transcript_name \"(.*?)\";'
    isonames = set()
    with open(gencode_gtf_file) as ifile:
        for line in ifile:
            if line.startswith('#'): continue
            wds = line.split('\t')
            if wds[2] == 'CDS':
                isoname_result = re.search(isoname_pattern, line)
                if isoname_result is not None:
                    isoname = isoname_result.group(1)
                    isonames.add(isoname)
    return isonames

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference_gtf", action="store", dest="reference_gtf",)
    args = parser.parse_args()

    isonames_with_cds = get_transcripts_with_cds(args.reference_gtf)
    filter_protein_coding_gtf(args.reference_gtf,isonames_with_cds)


if __name__ == "__main__":
    main()