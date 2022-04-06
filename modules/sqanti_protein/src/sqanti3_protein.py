#!/usr/bin/env python
# SQANTI: Structural and Quality Annotation of Novel Transcript Isoforms
# Authors: Lorena de la Fuente, Hector del Risco, Cecile Pereira and Manuel Tardaguila
# Modified by Liz (etseng@pacb.com) as SQANTI2/3 versions
# Modified by Fran (francisco.pardo.palacios@gmail.com) currently as SQANTI3 version (05/15/2020)
# Modified by Gloria (gs9yr@virginia.edu) for protein classification

# %%
__author__  = "etseng@pacb.com"
__version__ = '2.0.0'  # Python 3.7

# import pdb
import os, re, sys, subprocess, timeit, glob, copy
# import shutil
# import distutils.spawn
import itertools
import bisect
import argparse
import math
import numpy as np
from scipy import mean
from collections import defaultdict, Counter, namedtuple
from collections.abc import Iterable
from csv import DictWriter, DictReader

utilitiesPath = os.path.join(os.path.dirname(os.path.realpath(__file__)), "utilities")
sys.path.insert(0, utilitiesPath)
# from multiprocessing import Process

# utilitiesPath = os.path.join(os.path.dirname(os.path.realpath(__file__)), "utilities")
# sys.path.insert(0, utilitiesPath)
# from rt_switching import rts
# from indels_annot import calc_indels_from_sam


try:
    from Bio.Seq import Seq
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
except ImportError:
    print("Unable to import Biopython! Please make sure Biopython is installed.", file=sys.stderr)
    sys.exit(-1)

try:
    from bx.intervals import Interval, IntervalTree
except ImportError:
    print("Unable to import bx-python! Please make sure bx-python is installed.", file=sys.stderr)
    sys.exit(-1)


# so gloria can run cupcake and gtf tools on mac
# cupcake_dir = '/Users/gloriasheynkman/Documents/research_drive/bioinfo_tools/cDNA_Cupcake/'
# if cupcake_dir not in sys.path:
#     sys.path.append('/Users/gloriasheynkman/Documents/research_drive/bioinfo_tools/cDNA_Cupcake/')
from cupcake.tofu.compare_junctions import compare_junctions
from cupcake.io.GFF import collapseGFFReader, write_collapseGFF_format
GTF2GENEPRED_PROG = "gtfToGenePred"
utilitiesPath = os.path.join(os.path.dirname(os.path.realpath(__file__)), "utilities")
sys.path.insert(0, utilitiesPath)
GTF2GENEPRED_PROG = os.path.join(utilitiesPath,"gtfToGenePred")
GFFREAD_PROG = "gffread"

# importing functions from original sqanti3
from sqanti3_qc import myQueryTranscripts, genePredReader, myQueryProteins, associationOverlapping


# functions that were modified to run "sqanti protein" (gloria and liz)

def reference_parser(args, genome_chroms):
    """
    Read the reference GTF file
    :param args:
    :param genome_chroms: list of chromosome names from the genome fasta, used for sanity checking
    :return: (refs_1exon_by_chr, refs_exons_by_chr, junctions_by_chr, junctions_by_gene)
    """
    global referenceFiles

    referenceFiles = os.path.join(args.dir, "refAnnotation_" + args.output_prefix + ".genePred")
    print("**** Parsing Reference Transcriptome....", file=sys.stdout)

    # don't check for exisitng reference file, because need to create a new one for exon and cds centric runs
    # if os.path.exists(referenceFiles):
    #     print("{0} already exists. Using it.".format(referenceFiles), file=sys.stdout)
    # else:
    ## gtf to genePred
    if not args.genename:
        subprocess.call([GTF2GENEPRED_PROG, args.annotation, referenceFiles, '-genePredExt', '-allErrors', '-ignoreGroupsWithoutExons'])
    else:
        subprocess.call([GTF2GENEPRED_PROG, args.annotation, referenceFiles, '-genePredExt', '-allErrors', '-ignoreGroupsWithoutExons', '-geneNameAsName2'])

    ## parse reference annotation
    # 1. ignore all miRNAs (< 200 bp)
    # 2. separately store single exon and multi-exon references
    refs_1exon_by_chr = defaultdict(lambda: IntervalTree()) #
    refs_exons_by_chr = defaultdict(lambda: IntervalTree())
    # store donors as the exon end (1-based) and acceptor as the exon start (0-based)
    # will convert the sets to sorted list later
    junctions_by_chr = defaultdict(lambda: {'donors': set(), 'acceptors': set(), 'da_pairs': set()})
    # dict of gene name --> set of junctions (don't need to record chromosome)
    junctions_by_gene = defaultdict(lambda: set())
    # dict of gene name --> list of known begins and ends (begin always < end, regardless of strand)
    known_5_3_by_gene = defaultdict(lambda: {'begin':set(), 'end': set()})

    ## dictionary of record.id (e.g., gencode enst) to genePred record objects
    ## need this for later computation of 5' and 3' overhangs for protein classification
    refDict = {}

    for r in genePredReader(referenceFiles):
        refDict[r.id] = r
        if r.length < args.min_ref_len and not args.is_fusion: continue # ignore miRNAs
        if r.exonCount == 1:
            refs_1exon_by_chr[r.chrom].insert(r.txStart, r.txEnd, r)
            known_5_3_by_gene[r.gene]['begin'].add(r.txStart)
            known_5_3_by_gene[r.gene]['end'].add(r.txEnd)
        else:
            refs_exons_by_chr[r.chrom].insert(r.txStart, r.txEnd, r)
            # only store junctions for multi-exon transcripts
            for d, a in r.junctions:
                junctions_by_chr[r.chrom]['donors'].add(d)
                junctions_by_chr[r.chrom]['acceptors'].add(a)
                junctions_by_chr[r.chrom]['da_pairs'].add((d,a))
                junctions_by_gene[r.gene].add((d,a))
            known_5_3_by_gene[r.gene]['begin'].add(r.txStart)
            known_5_3_by_gene[r.gene]['end'].add(r.txEnd)

    # check that all genes' chromosomes are in the genome file
    ref_chroms = set(refs_1exon_by_chr.keys()).union(list(refs_exons_by_chr.keys()))
    if genome_chroms is not None:
        diff = ref_chroms.difference(genome_chroms)
        if len(diff) > 0:
            print("WARNING: ref annotation contains chromosomes not in genome: {0}\n".format(",".join(diff)), file=sys.stderr)

    # convert the content of junctions_by_chr to sorted list
    for k in junctions_by_chr:
        junctions_by_chr[k]['donors'] = list(junctions_by_chr[k]['donors'])
        junctions_by_chr[k]['donors'].sort()
        junctions_by_chr[k]['acceptors'] = list(junctions_by_chr[k]['acceptors'])
        junctions_by_chr[k]['acceptors'].sort()
        junctions_by_chr[k]['da_pairs'] = list(junctions_by_chr[k]['da_pairs'])
        junctions_by_chr[k]['da_pairs'].sort()

    return dict(refs_1exon_by_chr), dict(refs_exons_by_chr), dict(junctions_by_chr), dict(junctions_by_gene), dict(known_5_3_by_gene), refDict


class myProteinTranscripts(myQueryTranscripts):
    def __init__(self, id, tss_diff, tts_diff, num_exons, length, str_class, subtype=None,
                 genes=None, transcripts=None, chrom=None, strand=None, bite ="NA",
                 RT_switching ="????", canonical="NA", min_cov ="NA",
                 min_cov_pos ="NA", min_samp_cov="NA", sd ="NA", FL ="NA", FL_dict={},
                 nIndels ="NA", nIndelsJunc ="NA", proteinID=None,
                 ORFlen="NA", CDS_start="NA", CDS_end="NA",
                 CDS_genomic_start="NA", CDS_genomic_end="NA", 
                 ORFseq="NA",
                 is_NMD="NA",
                 isoExp ="NA", geneExp ="NA", coding ="non_coding",
                 refLen ="NA", refExons ="NA",
                 refStart = "NA", refEnd = "NA",
                 q_splicesite_hit = 0,
                 q_exon_overlap = 0,
                 FSM_class = None, percAdownTTS = None, seqAdownTTS=None,
                 dist_cage='NA', within_cage='NA',
                 dist_polya_site='NA', within_polya_site='NA',
                 polyA_motif='NA', polyA_dist='NA', ref_obj=None):

        super().__init__(id, tss_diff=tss_diff,
                         tts_diff=tts_diff,
                         num_exons=num_exons,
                         length=length,
                         str_class=str_class,
                         subtype=subtype,
                 genes=genes, transcripts=transcripts, chrom=chrom, strand=strand, bite=bite,
                 RT_switching=RT_switching, canonical=canonical, min_cov=min_cov,
                 min_cov_pos=min_cov_pos, min_samp_cov=min_samp_cov, sd=sd, FL=FL, FL_dict=FL_dict,
                 nIndels=nIndels, nIndelsJunc=nIndelsJunc , proteinID=proteinID,
                 ORFlen=ORFlen, CDS_start=CDS_start, CDS_end=CDS_end,
                 CDS_genomic_start=CDS_genomic_start, CDS_genomic_end=CDS_genomic_end,
                 ORFseq=ORFseq, is_NMD=is_NMD,
                 isoExp=isoExp, geneExp=geneExp, coding=coding,
                 refLen=refLen, refExons=refExons,
                 refStart=refStart, refEnd=refEnd,
                 q_splicesite_hit=q_splicesite_hit,
                 q_exon_overlap=q_exon_overlap,
                 FSM_class=FSM_class, percAdownTTS=percAdownTTS, seqAdownTTS=seqAdownTTS,
                 dist_cage=dist_cage, within_cage=within_cage,
                 dist_polya_site=dist_polya_site, within_polya_site=within_polya_site,
                 polyA_motif=polyA_motif, polyA_dist=polyA_dist)

        self.ref_obj = ref_obj
        self.num_junc_after_stop = None

    def modify(self, ref_transcript, ref_gene, tss_diff, tts_diff, refLen, refExons, ref):
        super().modify(ref_transcript=ref_transcript,
                         ref_gene=ref_gene,
                         tss_diff=tss_diff,
                         tts_diff=tts_diff,
                         refLen=refLen,
                         refExons=refExons)
        self.ref_obj = ref

    def __str__(self):
        return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.chrom, self.strand,
                                                                                                                                                           str(self.length), str(self.num_exons),
                                                                                                                                                           str(self.str_class), "_".join(set(self.genes)),
                                                                                                                                                           self.id, str(self.refLen), str(self.refExons),
                                                                                                                                                           str(self.tss_diff), str(self.tts_diff),
                                                                                                                                                           self.subtype, self.RT_switching,
                                                                                                                                                           self.canonical, str(self.min_samp_cov),
                                                                                                                                                           str(self.min_cov), str(self.min_cov_pos),
                                                                                                                                                           str(self.sd), str(self.FL), str(self.nIndels),
                                                                                                                                                           str(self.nIndelsJunc), self.bite, str(self.isoExp),
                                                                                                                                                           str(self.geneExp), str(self.ratioExp()),
                                                                                                                                                           self.FSM_class, self.coding, str(self.ORFlen),
                                                                                                                                                           str(self.CDSlen()), str(self.CDS_start), str(self.CDS_end),
                                                                                                                                                           str(self.CDS_genomic_start), str(self.CDS_genomic_end), str(self.is_NMD),
                                                                                                                                                           str(self.percAdownTTS),
                                                                                                                                                           str(self.seqAdownTTS),
                                                                                                                                                           str(self.dist_cage),
                                                                                                                                                           str(self.within_cage),str(self.pos_cage_peak),
                                                                                                                                                           str(self.dist_polya_site),
                                                                                                                                                           str(self.within_polya_site),
                                                                                                                                                           str(self.polyA_motif),
                                                                                                                                                           str(self.polyA_dist))




def protein_isoforms_parser(args):
    """
    Parse input isoforms (GTF) to dict (chr --> sorted list)
    """
    queryFile = os.path.splitext(args.corrGTF)[0] +".genePred"

    print("**** Parsing Isoforms....", file=sys.stderr)

    # gtf to genePred
    cmd = GTF2GENEPRED_PROG + " {0} {1} -genePredExt -allErrors -ignoreGroupsWithoutExons".format(\
        args.corrGTF, queryFile)
    if subprocess.check_call(cmd, shell=True)!=0:
        print("ERROR running cmd: {0}".format(cmd), file=sys.stderr)
        sys.exit(-1)

    isoforms_list = defaultdict(lambda: []) # chr --> list to be sorted later

    # FOR PROTEINS only
    # to compute 5' and 3' overhang for protein classification, need to retrieve the genePred object later
    queryDict = {}
    for r in genePredReader(queryFile):
        isoforms_list[r.chrom].append(r)
        queryDict[r.id] = r
    for k in isoforms_list:
        isoforms_list[k].sort(key=lambda r: r.txStart)

    return isoforms_list, queryDict


def transcriptsKnownSpliceSites(refs_1exon_by_chr, refs_exons_by_chr, start_ends_by_gene, trec, genome_dict, nPolyA):
    """
    :param refs_1exon_by_chr: dict of single exon references (chr -> IntervalTree)
    :param refs_exons_by_chr: dict of multi exon references (chr -> IntervalTree)
    :param trec: id record (genePredRecord) to be compared against reference
    :param genome_dict: dict of genome (chrom --> SeqRecord)
    :param nPolyA: window size to look for polyA
    :return: myQueryTranscripts object that indicates the best reference hit
    """
    def calc_overlap(s1, e1, s2, e2):
        if s1=='NA' or s2=='NA': return 0
        if s1 > s2:
            s1, e1, s2, e2 = s2, e2, s1, e1
        return max(0, min(e1,e2)-max(s1,s2))

    def calc_splicesite_agreement(query_exons, ref_exons):
        q_sites = {}
        for e in query_exons:
            q_sites[e.start] = 0
            q_sites[e.end] = 0
        for e in ref_exons:
            if e.start in q_sites: q_sites[e.start] = 1
            if e.end in q_sites: q_sites[e.end] = 1
        return sum(q_sites.values())

    def calc_exon_overlap(query_exons, ref_exons):
        q_bases = {}
        for e in query_exons:
            for b in range(e.start, e.end): q_bases[b] = 0

        for e in ref_exons:
            for b in range(e.start, e.end):
                if b in q_bases: q_bases[b] = 1
        return sum(q_bases.values())

    def get_diff_tss_tts(trec, ref):
        if trec.strand == '+':
            diff_tss = trec.txStart - ref.txStart
            diff_tts = ref.txEnd - trec.txEnd
        else:
            diff_tts = trec.txStart - ref.txStart
            diff_tss = ref.txEnd - trec.txEnd
        return diff_tss, diff_tts


    def get_gene_diff_tss_tts(isoform_hit):
        # now that we know the reference (isoform) it hits
        # add the nearest start/end site for that gene (all isoforms of the gene)
        nearest_start_diff, nearest_end_diff = float('inf'), float('inf')
        for ref_gene in isoform_hit.genes:
            for x in start_ends_by_gene[ref_gene]['begin']:
                d = trec.txStart - x
                if abs(d) < abs(nearest_start_diff):
                    nearest_start_diff = d
            for x in start_ends_by_gene[ref_gene]['end']:
                d = trec.txEnd - x
                if abs(d) < abs(nearest_end_diff):
                    nearest_end_diff = d

        if trec.strand == '+':
            isoform_hit.tss_gene_diff = nearest_start_diff if nearest_start_diff!=float('inf') else 'NA'
            isoform_hit.tts_gene_diff = nearest_end_diff if nearest_end_diff!=float('inf') else 'NA'
        else:
            isoform_hit.tss_gene_diff = -nearest_end_diff if nearest_start_diff!=float('inf') else 'NA'
            isoform_hit.tts_gene_diff = -nearest_start_diff if nearest_end_diff!=float('inf') else 'NA'

    def categorize_incomplete_matches(trec, ref):
        """
        intron_retention --- at least one trec exon covers at least two adjacent ref exons
        complete --- all junctions agree and is not IR
        5prime_fragment --- all junctions agree but trec has less 5' exons. The isoform is a 5' fragment of the reference transcript
        3prime_fragment --- all junctions agree but trec has less 3' exons. The isoform is a 3' fragment of the reference transcript
        internal_fragment --- all junctions agree but trec has less 5' and 3' exons
        """
        # check intron retention
        ref_exon_tree = IntervalTree()
        for i,e in enumerate(ref.exons): ref_exon_tree.insert(e.start, e.end, i)
        for e in trec.exons:
            if len(ref_exon_tree.find(e.start, e.end)) > 1: # multiple ref exons covered
                return "intron_retention"

        agree_front = trec.junctions[0]==ref.junctions[0]
        agree_end   = trec.junctions[-1]==ref.junctions[-1]
        if agree_front:
            if agree_end:
                return "complete"
            else: # front agrees, end does not
                return ("5prime_fragment" if trec.strand=='+' else '3prime_fragment')
        else:
            if agree_end: # front does not agree, end agrees
                return ("3prime_fragment" if trec.strand=='+' else '5prime_fragment')
            else:
                return "internal_fragment"


    # NOTE liz - not using polya annotations
    seq_downTTS = None
    percA = None

    isoform_hit = myProteinTranscripts(id=trec.id, tts_diff="NA", tss_diff="NA", \
                                       num_exons=trec.exonCount,
                                       length=trec.length,
                                       str_class="", \
                                       chrom=trec.chrom,
                                       strand=trec.strand, \
                                       subtype="no_subcategory", \
                                       percAdownTTS=str(percA), \
                                       seqAdownTTS=seq_downTTS)

    ##***************************************##
    ########### SPLICED TRANSCRIPTS ###########
    ##***************************************##

    cat_ranking = {'full-splice_match': 5, 'incomplete-splice_match': 4, 'anyKnownJunction': 3, 'anyKnownSpliceSite': 2,
                   'geneOverlap': 1, '': 0}

    #if trec.id.startswith('PB.1961.2'):
    #    pdb.set_trace()
    if trec.exonCount >= 2:

        hits_by_gene = defaultdict(lambda: [])  # gene --> list of hits
        best_by_gene = {}  # gene --> best isoform_hit

        if trec.chrom in refs_exons_by_chr:
            for ref in refs_exons_by_chr[trec.chrom].find(trec.txStart, trec.txEnd):
                hits_by_gene[ref.gene].append(ref)
        if trec.chrom in refs_1exon_by_chr:
            for ref in refs_1exon_by_chr[trec.chrom].find(trec.txStart, trec.txEnd):
                hits_by_gene[ref.gene].append(ref)

        if len(hits_by_gene) == 0: return isoform_hit

        for ref_gene in hits_by_gene:
            isoform_hit = myProteinTranscripts(id=trec.id, tts_diff="NA", tss_diff="NA", \
                                               num_exons=trec.exonCount,
                                               length=trec.length,
                                               str_class="", \
                                               chrom=trec.chrom,
                                               strand=trec.strand, \
                                               subtype="no_subcategory", \
                                               percAdownTTS=str(percA), \
                                               seqAdownTTS=seq_downTTS)

            for ref in hits_by_gene[ref_gene]:
                if trec.strand != ref.strand:
                    # opposite strand, just record it in AS_genes
                    isoform_hit.AS_genes.add(ref.gene)
                    continue

                if ref.exonCount == 1: # mono-exonic reference, handle specially here
                    if calc_exon_overlap(trec.exons, ref.exons) > 0 and cat_ranking[isoform_hit.str_class] < cat_ranking["geneOverlap"]:
                        isoform_hit = myProteinTranscripts(trec.id, "NA", "NA", trec.exonCount, trec.length,
                                                            "geneOverlap",
                                                           subtype="mono-exon",
                                                           chrom=trec.chrom,
                                                           strand=trec.strand,
                                                           genes=[ref.gene],
                                                           transcripts=[ref.id],
                                                           refLen=ref.length,
                                                           refExons=ref.exonCount,
                                                           refStart=ref.txStart,
                                                           refEnd=ref.txEnd,
                                                           q_splicesite_hit=0,
                                                           q_exon_overlap=calc_exon_overlap(trec.exons, ref.exons),
                                                           percAdownTTS=str(percA),
                                                           seqAdownTTS=seq_downTTS,
                                                           ref_obj=ref)

                else: # multi-exonic reference
                    match_type = compare_junctions(trec, ref, internal_fuzzy_max_dist=0, max_5_diff=999999, max_3_diff=999999)

                    if match_type not in ('exact', 'subset', 'partial', 'concordant', 'super', 'nomatch'):
                        raise Exception("Unknown match category {0}!".format(match_type))

                    diff_tss, diff_tts = get_diff_tss_tts(trec, ref)
                    #has_overlap = gene_overlap(isoform_hit.genes[-1], ref.gene) if len(isoform_hit.genes) >= 1 else Fals
                    # #############################
                    # SQANTI's full-splice_match
                    # #############################
                    if match_type == "exact":
                        subtype = "multi-exon"
                        # assign as a new hit if
                        # (1) no prev hits yet
                        # (2) this one is better (prev not FSM or is FSM but worse tss/tts)
                        if cat_ranking[isoform_hit.str_class] < cat_ranking["full-splice_match"] or \
                                abs(diff_tss)+abs(diff_tts) < isoform_hit.get_total_diff():
                            isoform_hit = myProteinTranscripts(trec.id, diff_tss, diff_tts, trec.exonCount, trec.length,
                                                               str_class="full-splice_match",
                                                               subtype=subtype,
                                                               chrom=trec.chrom,
                                                               strand=trec.strand,
                                                               genes=[ref.gene],
                                                               transcripts=[ref.id],
                                                               refLen = ref.length,
                                                               refExons= ref.exonCount,
                                                               refStart=ref.txStart,
                                                               refEnd=ref.txEnd,
                                                               q_splicesite_hit=calc_splicesite_agreement(trec.exons, ref.exons),
                                                               q_exon_overlap=calc_exon_overlap(trec.exons, ref.exons),
                                                               percAdownTTS=str(percA),
                                                               seqAdownTTS=seq_downTTS,
                                                               ref_obj=ref)
                    # #######################################################
                    # SQANTI's incomplete-splice_match
                    # (only check if don't already have a FSM match)
                    # #######################################################
                    elif match_type == "subset":
                        subtype = categorize_incomplete_matches(trec, ref)
                        # assign as a new (ISM) hit if
                        # (1) no prev hit
                        # (2) prev hit not as good (is ISM with worse tss/tts or anyKnownSpliceSite)
                        if cat_ranking[isoform_hit.str_class] < cat_ranking["incomplete-splice_match"] or \
                            (isoform_hit.str_class=='incomplete-splice_match' and abs(diff_tss)+abs(diff_tts) < isoform_hit.get_total_diff()):
                            isoform_hit = myProteinTranscripts(trec.id, diff_tss, diff_tts, trec.exonCount, trec.length,
                                                               str_class="incomplete-splice_match",
                                                               subtype=subtype,
                                                               chrom=trec.chrom,
                                                               strand=trec.strand,
                                                               genes=[ref.gene],
                                                               transcripts=[ref.id],
                                                               refLen = ref.length,
                                                               refExons= ref.exonCount,
                                                               refStart=ref.txStart,
                                                               refEnd=ref.txEnd,
                                                               q_splicesite_hit=calc_splicesite_agreement(trec.exons, ref.exons),
                                                               q_exon_overlap=calc_exon_overlap(trec.exons, ref.exons),
                                                               percAdownTTS=str(percA),
                                                               seqAdownTTS=seq_downTTS,
                                                               ref_obj=ref)
                    # #######################################################
                    # Some kind of junction match that isn't ISM/FSM
                    # #######################################################
                    elif match_type in ('partial', 'concordant', 'super'):
                        q_sp_hit = calc_splicesite_agreement(trec.exons, ref.exons)
                        q_ex_overlap = calc_exon_overlap(trec.exons, ref.exons)
                        q_exon_d = abs(trec.exonCount - ref.exonCount)
                        if cat_ranking[isoform_hit.str_class] < cat_ranking["anyKnownJunction"] or \
                                (isoform_hit.str_class=='anyKnownJunction' and q_sp_hit > isoform_hit.q_splicesite_hit) or \
                                (isoform_hit.str_class=='anyKnownJunction' and q_sp_hit==isoform_hit.q_splicesite_hit and q_ex_overlap > isoform_hit.q_exon_overlap) or \
                                (isoform_hit.str_class=='anyKnownJunction' and q_sp_hit==isoform_hit.q_splicesite_hit and q_exon_d < abs(trec.exonCount-isoform_hit.refExons)):
                            isoform_hit = myProteinTranscripts(trec.id, "NA", "NA", trec.exonCount, trec.length,
                                                               str_class="anyKnownJunction",
                                                               subtype="no_subcategory",
                                                               chrom=trec.chrom,
                                                               strand=trec.strand,
                                                               genes=[ref.gene],
                                                               transcripts=["novel"],
                                                               refLen=ref.length,
                                                               refExons=ref.exonCount,
                                                               refStart=ref.txStart,
                                                               refEnd=ref.txEnd,
                                                               q_splicesite_hit=calc_splicesite_agreement(trec.exons, ref.exons),
                                                               q_exon_overlap=calc_exon_overlap(trec.exons, ref.exons),
                                                               percAdownTTS=str(percA),
                                                               seqAdownTTS=seq_downTTS,
                                                               ref_obj=ref)
                    else: # must be nomatch
                        assert match_type == 'nomatch'
                        # at this point, no junction overlap, but may be a single splice site (donor or acceptor) match?
                        # also possibly just exonic (no splice site) overlap
                        if cat_ranking[isoform_hit.str_class] < cat_ranking["anyKnownSpliceSite"] and calc_splicesite_agreement(trec.exons, ref.exons) > 0:
                            isoform_hit = myProteinTranscripts(trec.id, "NA", "NA", trec.exonCount, trec.length,
                                                               str_class="anyKnownSpliceSite",
                                                               subtype="no_subcategory",
                                                               chrom=trec.chrom,
                                                               strand=trec.strand,
                                                               genes=[ref.gene],
                                                               transcripts=["novel"],
                                                               refLen=ref.length,
                                                               refExons=ref.exonCount,
                                                               refStart=ref.txStart,
                                                               refEnd=ref.txEnd,
                                                               q_splicesite_hit=calc_splicesite_agreement(trec.exons, ref.exons),
                                                               q_exon_overlap=calc_exon_overlap(trec.exons,
                                                                                              ref.exons),
                                                               percAdownTTS=str(percA),
                                                               seqAdownTTS=seq_downTTS)

                        if isoform_hit.str_class=="": # still not hit yet, check exonic overlap
                            if cat_ranking[isoform_hit.str_class] < cat_ranking["geneOverlap"] and calc_exon_overlap(trec.exons, ref.exons) > 0:
                                isoform_hit = myProteinTranscripts(trec.id, "NA", "NA", trec.exonCount, trec.length,
                                                                   str_class="geneOverlap",
                                                                   subtype="no_subcategory",
                                                                   chrom=trec.chrom,
                                                                   strand=trec.strand,
                                                                   genes=[ref.gene],
                                                                   transcripts=["novel"],
                                                                   refLen=ref.length,
                                                                   refExons=ref.exonCount,
                                                                   refStart=ref.txStart,
                                                                   refEnd=ref.txEnd,
                                                                   q_splicesite_hit=calc_splicesite_agreement(trec.exons, ref.exons),
                                                                   q_exon_overlap=calc_exon_overlap(trec.exons, ref.exons),
                                                                   percAdownTTS=str(percA),
                                                                   seqAdownTTS=seq_downTTS)

            best_by_gene[ref_gene] = isoform_hit
        # now we have best_by_gene:
        # start with the best scoring one (FSM is best) --> can add other genes if they don't overlap
        #if trec.id.startswith('PB.1252.'):
        #    pdb.set_trace()
        geneHitTuple = namedtuple('geneHitTuple', ['score', 'rStart', 'rEnd', 'rGene', 'iso_hit'])
        best_by_gene = [geneHitTuple(cat_ranking[iso_hit.str_class],iso_hit.refStart,iso_hit.refEnd,ref_gene,iso_hit) for ref_gene,iso_hit in best_by_gene.items()]
        best_by_gene = list(filter(lambda x: x.score > 0, best_by_gene))
        if len(best_by_gene) == 0: # no hit
            return isoform_hit

        # sort matching genes by ranking, allow for multi-gene match as long as they don't overlap
        # cat_ranking = {'full-splice_match': 5, 'incomplete-splice_match': 4, 'anyKnownJunction': 3, 'anyKnownSpliceSite': 2,
        #                    'geneOverlap': 1, '': 0}

        best_by_gene.sort(key=lambda x: (x.score,x.iso_hit.q_splicesite_hit+(x.iso_hit.q_exon_overlap)*1./sum(e.end-e.start for e in trec.exons)+calc_overlap(x.rStart,x.rEnd,trec.txStart,trec.txEnd)*1./(x.rEnd-x.rStart)-abs(trec.exonCount-x.iso_hit.refExons)), reverse=True)  # sort by (ranking score, overlap)
        isoform_hit = best_by_gene[0].iso_hit
        cur_start, cur_end = best_by_gene[0].rStart, best_by_gene[0].rEnd
        for t in best_by_gene[1:]:
            if t.score==0: break
            if calc_overlap(cur_start, cur_end, t.rStart, t.rEnd) <= 0:
                isoform_hit.genes.append(t.rGene)
                cur_start, cur_end = min(cur_start, t.rStart), max(cur_end, t.rEnd)

    ##***************************************####
    ########### UNSPLICED TRANSCRIPTS ###########
    ##***************************************####
    else: # single exon id
        if trec.chrom in refs_1exon_by_chr:
            for ref in refs_1exon_by_chr[trec.chrom].find(trec.txStart, trec.txEnd):
                if ref.strand != trec.strand:
                    # opposite strand, just record it in AS_genes
                    isoform_hit.AS_genes.add(ref.gene)
                    continue
                diff_tss, diff_tts = get_diff_tss_tts(trec, ref)

                # see if there's already an existing match AND if so, if this one is better
                if isoform_hit.str_class == "": # no match so far
                    isoform_hit = myProteinTranscripts(trec.id, diff_tss, diff_tts, trec.exonCount, trec.length, "full-splice_match",
                                                       subtype="mono-exon",
                                                       chrom=trec.chrom,
                                                       strand=trec.strand,
                                                       genes=[ref.gene],
                                                       transcripts=[ref.id],
                                                       refLen=ref.length,
                                                       refExons = ref.exonCount,
                                                       percAdownTTS=str(percA),
                                                       seqAdownTTS=seq_downTTS,
                                                       ref_obj=ref)
                elif abs(diff_tss)+abs(diff_tts) < isoform_hit.get_total_diff():
                    isoform_hit.modify(ref.id, ref.gene, diff_tss, diff_tts, ref.length, ref.exonCount, ref)

        if isoform_hit.str_class == "" and trec.chrom in refs_exons_by_chr:
            # no hits to single exon genes, let's see if it hits multi-exon genes
            # (1) if it overlaps with a ref exon and is contained in an exon, we call it ISM
            # (2) else, if it is completely within a ref gene start-end region, we call it NIC by intron retention
            for ref in refs_exons_by_chr[trec.chrom].find(trec.txStart, trec.txEnd):
                if calc_exon_overlap(trec.exons, ref.exons) == 0:   # no exonic overlap, skip!
                    continue
                if ref.strand != trec.strand:
                    # opposite strand, just record it in AS_genes
                    isoform_hit.AS_genes.add(ref.gene)
                    continue
                diff_tss, diff_tts = get_diff_tss_tts(trec, ref)

                for e in ref.exons:
                    if e.start <= trec.txStart < trec.txEnd <= e.end:
                        isoform_hit.str_class = "incomplete-splice_match"
                        isoform_hit.subtype = "mono-exon"
                        isoform_hit.modify(ref.id, ref.gene, diff_tss, diff_tts, ref.length, ref.exonCount, ref)
                        # this is as good a match as it gets, we can stop the search here
                        get_gene_diff_tss_tts(isoform_hit)
                        return isoform_hit

                # if we haven't exited here, then ISM hit is not found yet
                # instead check if it's NIC by intron retention
                # but we don't exit here since the next gene could be a ISM hit
                if ref.txStart <= trec.txStart < trec.txEnd <= ref.txEnd:
                    isoform_hit.str_class = "novel_in_catalog"
                    isoform_hit.subtype = "mono-exon"
                    # check for intron retention
                    if len(ref.junctions) > 0:
                        for (d,a) in ref.junctions:
                            if trec.txStart < d < a < trec.txEnd:
                                isoform_hit.subtype = "mono-exon_by_intron_retention"
                                break
                    isoform_hit.modify("novel", ref.gene, 'NA', 'NA', ref.length, ref.exonCount, ref)
                    get_gene_diff_tss_tts(isoform_hit)
                    return isoform_hit

                # if we get to here, means neither ISM nor NIC, so just add a ref gene and categorize further later
                isoform_hit.genes.append(ref.gene)

    get_gene_diff_tss_tts(isoform_hit)
    isoform_hit.genes.sort(key=lambda x: start_ends_by_gene[x]['begin'])
    return isoform_hit


def novelIsoformsKnownGenes(isoforms_hit, trec, junctions_by_chr, junctions_by_gene, start_ends_by_gene):
    """
    At this point: definitely not FSM or ISM, see if it is NIC, NNC, or fusion
    :return isoforms_hit: updated isoforms hit (myQueryTranscripts object)
    """
    def has_intron_retention():
        for e in trec.exons:
            m = bisect.bisect_left(junctions_by_chr[trec.chrom]['da_pairs'], (e.start, e.end))
            if m < len(junctions_by_chr[trec.chrom]['da_pairs']) and e.start <= junctions_by_chr[trec.chrom]['da_pairs'][m][0] < junctions_by_chr[trec.chrom]['da_pairs'][m][1] < e.end:
                return True
        return False

    ref_genes = list(set(isoforms_hit.genes))

    #if trec.id.startswith('PB.37872'):
    #pdb.set_trace()
    #
    # at this point, we have already found matching genes/transcripts
    # hence we do not need to update refLen or refExon
    # or tss_diff and tts_diff (always set to "NA" for non-FSM/ISM matches)
    #
    isoforms_hit.transcripts = ["novel"]
    if len(ref_genes) == 1:
        # hits exactly one gene, must be either NIC or NNC
        ref_gene_junctions = junctions_by_gene[ref_genes[0]]
        # 1. check if all donors/acceptor sites are known (regardless of which ref gene it came from)
        # 2. check if this query isoform uses a subset of the junctions from the single ref hit
        all_junctions_known = True
        all_junctions_in_hit_ref = True
        for d,a in trec.junctions:
            all_junctions_known = all_junctions_known and (d in junctions_by_chr[trec.chrom]['donors']) and (a in junctions_by_chr[trec.chrom]['acceptors'])
            all_junctions_in_hit_ref = all_junctions_in_hit_ref and ((d,a) in ref_gene_junctions)
        if all_junctions_known:
            isoforms_hit.str_class="novel_in_catalog"
            if all_junctions_in_hit_ref:
                isoforms_hit.subtype = "combination_of_known_junctions"
            else:
                isoforms_hit.subtype = "combination_of_known_splicesites"
        else:
            isoforms_hit.str_class="novel_not_in_catalog"
            isoforms_hit.subtype = "at_least_one_novel_splicesite"
    else: # see if it is fusion
        # list of a ref junctions from all genes, including potential shared junctions
        # NOTE: some ref genes could be mono-exonic so no junctions
        all_ref_junctions = list(itertools.chain(junctions_by_gene[ref_gene] for ref_gene in ref_genes if ref_gene in junctions_by_gene))

        # (junction index) --> number of refs that have this junction
        junction_ref_hit = dict((i, all_ref_junctions.count(junc)) for i,junc in enumerate(trec.junctions))

        # if the same query junction appears in more than one of the hit references, it is not a fusion
        if max(junction_ref_hit.values()) > 1:
            isoforms_hit.str_class = "moreJunctions"
        else:
            isoforms_hit.str_class = "fusion"
            isoforms_hit.subtype = "mono-exon" if trec.exonCount==1 else "multi-exon"

    if has_intron_retention():
        isoforms_hit.subtype = "intron_retention"

    return isoforms_hit


def isoformClassification(args, isoforms_by_chr, refs_1exon_by_chr, refs_exons_by_chr, junctions_by_chr, junctions_by_gene, start_ends_by_gene, genome_dict, indelsJunc, orfDict):
    # if args.is_fusion: # read GFF to get fusion components
    #     # ex: PBfusion.1.1 --> (1-based start, 1-based end) of where the fusion component is w.r.t to entire fusion
    #     fusion_components = get_fusion_component(args.isoforms)

    # ## read coverage files if provided
    # if args.coverage is not None:
    #     print("**** Reading Splice Junctions coverage files.", file=sys.stdout)
    #     SJcovNames, SJcovInfo = STARcov_parser(args.coverage)
    #     fields_junc_cur = FIELDS_JUNC + SJcovNames # add the samples to the header
    # else:
    #     SJcovNames, SJcovInfo = None, None
    #     print("Splice Junction Coverage files not provided.", file=sys.stdout)
    #     fields_junc_cur = FIELDS_JUNC

    # if args.cage_peak is not None:
    #     print("**** Reading CAGE Peak data.", file=sys.stdout)
    #     cage_peak_obj = CAGEPeak(args.cage_peak)
    # else:
    #     cage_peak_obj = None

    # if args.polyA_peak is not None:
    #     print("**** Reading polyA Peak data.", file=sys.stdout)
    #     polya_peak_obj = PolyAPeak(args.polyA_peak)
    # else:
    #     polya_peak_obj = None

    # if args.polyA_motif_list is not None:
    #     print("**** Reading PolyA motif list.", file=sys.stdout)
    #     polyA_motif_list = []
    #     for line in open(args.polyA_motif_list):
    #         x = line.strip().upper().replace('U', 'A')
    #         if any(s not in ('A','T','C','G') for s in x):
    #             print("PolyA motif must be A/T/C/G only! Saw: {0}. Abort!".format(x), file=sys.stderr)
    #             sys.exit(-1)
    #         polyA_motif_list.append(x)
    # else:
    #     polyA_motif_list = None


    # if args.phyloP_bed is not None:
    #     print("**** Reading PhyloP BED file.", file=sys.stdout)
    #     phyloP_reader = LazyBEDPointReader(args.phyloP_bed)
    # else:
    #     phyloP_reader = None

    # running classification
    print("**** Performing Classification of Isoforms....", file=sys.stdout)


    # accepted_canonical_sites = list(args.sites.split(","))

    # handle_class = open(outputClassPath+"_tmp", "w")
    # fout_class = DictWriter(handle_class, fieldnames=FIELDS_CLASS, delimiter='\t')
    # fout_class.writeheader()

    # #outputJuncPath = outputPathPrefix+"_junctions.txt"
    # handle_junc = open(outputJuncPath+"_tmp", "w")
    # fout_junc = DictWriter(handle_junc, fieldnames=fields_junc_cur, delimiter='\t')
    # fout_junc.writeheader()

    isoforms_info = {}
    novel_gene_index = 1

    for chrom,records in isoforms_by_chr.items():
        for rec in records:
            # Find best reference hit
            isoform_hit = transcriptsKnownSpliceSites(refs_1exon_by_chr, refs_exons_by_chr, start_ends_by_gene, rec, genome_dict, nPolyA=args.window)

            if isoform_hit.str_class in ("anyKnownJunction", "anyKnownSpliceSite"):
                # not FSM or ISM --> see if it is NIC, NNC, or fusion
                isoform_hit = novelIsoformsKnownGenes(isoform_hit, rec, junctions_by_chr, junctions_by_gene, start_ends_by_gene)
            elif isoform_hit.str_class in ("", "geneOverlap"):
                # possibly NNC, genic, genic intron, anti-sense, or intergenic
                isoform_hit = associationOverlapping(isoform_hit, rec, junctions_by_chr)

            # # write out junction information
            # write_junctionInfo(rec, junctions_by_chr, accepted_canonical_sites, indelsJunc, genome_dict, fout_junc, covInf=SJcovInfo, covNames=SJcovNames, phyloP_reader=phyloP_reader)

            if isoform_hit.str_class in ("intergenic", "genic_intron"):
                # Liz: I don't find it necessary to cluster these novel genes. They should already be always non-overlapping.
                if args.novel_gene_prefix is not None:  # used by splits to not have redundant novelGene IDs
                    isoform_hit.genes = ['novelGene_' + str(args.novel_gene_prefix) + '_' + str(novel_gene_index)]
                else:
                    isoform_hit.genes = ['novelGene_' + str(novel_gene_index)]
                isoform_hit.transcripts = ['novel']
                novel_gene_index += 1

            # # look at Cage Peak info (if available)
            # if cage_peak_obj is not None:
            #     if rec.strand == '+':
            #         within_cage, dist_cage , pos_cage_peak = cage_peak_obj.find(rec.chrom, rec.strand, rec.txStart)
            #     else:
            #         within_cage, dist_cage , pos_cage_peak = cage_peak_obj.find(rec.chrom, rec.strand, rec.txEnd)
            #     isoform_hit.within_cage = within_cage
            #     isoform_hit.dist_cage = dist_cage
            #     isoform_hit.pos_cage_peak = pos_cage_peak

            # # look at PolyA Peak info (if available)
            # if polya_peak_obj is not None:
            #     if rec.strand == '+':
            #         within_polya_site, dist_polya_site = polya_peak_obj.find(rec.chrom, rec.strand, rec.txStart)
            #     else:
            #         within_polya_site, dist_polya_site = polya_peak_obj.find(rec.chrom, rec.strand, rec.txEnd)
            #     isoform_hit.within_polya_site = within_polya_site
            #     isoform_hit.dist_polya_site = dist_polya_site

            # # polyA motif finding: look within 50 bp upstream of 3' end for the highest ranking polyA motif signal (user provided)
            # if polyA_motif_list is not None:
            #     if rec.strand == '+':
            #         polyA_motif, polyA_dist = find_polyA_motif(str(genome_dict[rec.chrom][rec.txEnd-50:rec.txEnd].seq), polyA_motif_list)
            #     else:
            #         polyA_motif, polyA_dist = find_polyA_motif(str(genome_dict[rec.chrom][rec.txStart:rec.txStart+50].reverse_complement().seq), polyA_motif_list)
            #     isoform_hit.polyA_motif = polyA_motif
            #     isoform_hit.polyA_dist = polyA_dist

            # Fill in ORF/coding info and NMD detection
            if orfDict:
                if args.is_fusion:
                    #pdb.set_trace()
                    # fusion - special case handling, need to see which part of the ORF this segment falls on
                    fusion_gene = 'PBfusion.' + str(seqid_fusion.match(rec.id).group(1))
                    rec_component_start, rec_component_end = fusion_components[rec.id]
                    rec_len = rec_component_end - rec_component_start + 1
                    if fusion_gene in orfDict:
                        orf_start, orf_end = orfDict[fusion_gene].cds_start, orfDict[fusion_gene].cds_end
                        if orf_start <= rec_component_start < orf_end:
                            isoform_hit.CDS_start = 1
                            isoform_hit.CDS_end = min(rec_len, orf_end - rec_component_start + 1)
                            isoform_hit.ORFlen = (isoform_hit.CDS_end - isoform_hit.CDS_start)/3
                            _s = (rec_component_start-orf_start)//3
                            _e = min(int(_s+isoform_hit.ORFlen), len(orfDict[fusion_gene].orf_seq))
                            isoform_hit.ORFseq = orfDict[fusion_gene].orf_seq[_s:_e]
                            isoform_hit.coding = "coding"
                        elif rec_component_start <= orf_start < rec_component_end:
                            isoform_hit.CDS_start = orf_start - rec_component_start
                            if orf_end >= rec_component_end:
                                isoform_hit.CDS_end = rec_component_end - rec_component_start + 1
                            else:
                                isoform_hit.CDS_end = orf_end - rec_component_start + 1
                            isoform_hit.ORFlen = (isoform_hit.CDS_end - isoform_hit.CDS_start) / 3
                            _e = min(int(isoform_hit.ORFlen), len(orfDict[fusion_gene].orf_seq))
                            isoform_hit.ORFseq = orfDict[fusion_gene].orf_seq[:_e]
                            isoform_hit.coding = "coding"
                elif rec.id in orfDict:  # this will never be true for fusion, so the above code seg runs instead
                    isoform_hit.coding = "coding"
                    isoform_hit.ORFlen = orfDict[rec.id].orf_length
                    isoform_hit.CDS_start = orfDict[rec.id].cds_start  # 1-based start
                    isoform_hit.CDS_end = orfDict[rec.id].cds_end      # 1-based end
                    isoform_hit.ORFseq  = orfDict[rec.id].orf_seq

            if isoform_hit.coding == "coding":
                m = {} # transcript coord (0-based) --> genomic coord (0-based)
                if rec.strand == '+':
                    i = 0
                    for exon in rec.exons:
                        for c in range(exon.start, exon.end):
                            m[i] = c
                            i += 1
                else: # - strand
                    i = 0
                    for exon in rec.exons:
                        for c in range(exon.start, exon.end):
                            m[rec.length-i-1] = c
                            i += 1

                isoform_hit.CDS_genomic_start = m[isoform_hit.CDS_start-1] + 1  # make it 1-based
                # NOTE: if using --orf_input, it is possible to see discrepancy between the exon structure
                # provided by GFF and the input ORF. For now, just shorten it
                isoform_hit.CDS_genomic_end = m[min(isoform_hit.CDS_end-1, max(m))] + 1    # make it 1-based
                #orfDict[rec.id].cds_genomic_start = m[orfDict[rec.id].cds_start-1] + 1  # make it 1-based
                #orfDict[rec.id].cds_genomic_end   = m[orfDict[rec.id].cds_end-1] + 1    # make it 1-based


            if isoform_hit.CDS_genomic_end!='NA':
                # NMD detection
                # if + strand, see if CDS stop is before the last junction
                if len(rec.junctions) > 0:
                    if rec.strand == '+':
                        dist_to_last_junc = isoform_hit.CDS_genomic_end - rec.junctions[-1][0]
                    else: # - strand
                        dist_to_last_junc = rec.junctions[0][1] - isoform_hit.CDS_genomic_end
                    isoform_hit.is_NMD = "TRUE" if dist_to_last_junc < 0 else "FALSE"
                    # can change dist_to_last_junct (above) to < 50, to match gencode nmd definition

            # find number of junctions downstream of stop codon
            # added as an ad hoc attribute of isoform_hit
            if isoform_hit.CDS_genomic_end != 'NA':
                num_junc_after_stop_codon = 0
                if rec.strand == '+':
                    for donor_coord, accept_coord in rec.junctions:
                        if donor_coord > isoform_hit.CDS_genomic_end:
                            num_junc_after_stop_codon += 1
                else: # - strand
                    for donor_coord, accept_coord in rec.junctions:
                        if accept_coord < isoform_hit.CDS_genomic_end:
                            num_junc_after_stop_codon += 1
                isoform_hit.num_junc_after_stop = num_junc_after_stop_codon

            isoforms_info[rec.id] = isoform_hit
            # fout_class.writerow(isoform_hit.as_dict())

    # handle_class.close()
    # handle_junc.close()
    return isoforms_info



def read_in_custom_orf_calls_into_orfDict(orf_file):
    # read in orfs called as part of lrp pipeline
    # stick with same format of orfDict as in expected in sqanti
    orfDict = {} # pb_acc -> myQueryProtein
    for line in open(orf_file):
        if line.startswith('pb_acc'): continue
        wds = line.split()
        pb_acc = wds[0]
        cds_start = int(wds[3])
        cds_end = int(wds[4])
        orf_length = int(wds[5])
        seq = ''
        transcript_len = int(wds[1])
        num_3utr_nt = transcript_len - orf_length
        orf_obj = myQueryProteins(cds_start, cds_end, orf_length, seq, pb_acc)
        orf_obj.num_3utr_nt = num_3utr_nt
        orfDict[pb_acc] = orf_obj 
    return orfDict


### code to get "perfect subset" information ###
### originally from cupcake compare_junctions.py ###

### begin ###

MatchIndexTuple = namedtuple('MatchIndexTuple', ['query_idx', 'ref_idx'])


def find_indices_for_exons_with_upstream_most_common_splicsite(r1, r2):
    for query_idx, qexon in enumerate(r1.segments):
        for ref_idx, rexon in enumerate(r2.segments):
            if rexon.end == qexon.end:
                return MatchIndexTuple(query_idx=query_idx, ref_idx=ref_idx)
    return MatchIndexTuple(query_idx=None, ref_idx=None)


def determine_extent_of_upstream_overhang(r1, r2):
    match = find_indices_for_exons_with_upstream_most_common_splicsite(r1, r2)
    if match is not None and match.ref_idx is not None:
        # pre: this exon in r1 (query) and r2 (ref) matches in the end
        # calculate overhang as the difference in the start
        overhang = r2.segments[match.ref_idx].start - r1.segments[match.query_idx].start
        return overhang
    else:
        return None


def determine_extent_of_downstream_overhang(r1, r2):
    # r2 is ref, and matched at the most 3' start site (if + strand)
    match = find_indices_for_exons_with_downstream_most_common_splicsite(r1, r2)
    if match is not None and match.ref_idx is not None:
        # pre: this exon in r1 (query) and r2 (ref) matches in the start
        # calculate overhang as the difference in the end
        overhang = r1.segments[match.query_idx].end - r2.segments[match.ref_idx].end
        return overhang
    else:
        return None


def find_indices_for_exons_with_downstream_most_common_splicsite(r1, r2):
    # r2 is the reference so we prioritze by search for the most 3' of r2 that can be matched by a r1 3' end
    for ref_idx in range(len(r2.segments) - 1, -1, -1):
        for query_idx in range(len(r1.segments) - 1, -1, -1):
            # the acceptor site (if + strand) matches
            # the donor site (if - strand) matces
            if r2.segments[ref_idx].start == r1.segments[query_idx].start:
                return MatchIndexTuple(query_idx=query_idx, ref_idx=ref_idx)


def get_perfect_subset_status(r1, r2):
    """
    Pre-req: r1, r2 overlap by some splice junction

    Check if r1 (query) is either a perfect or subset of r2 (query)

    :param r1: query isoform (r1.segments are the exons)
    :param r2: ref isoform (r2.segments are the exons)
    :return: (5_overhang_diff, 3_overhang_diff)
    """
    assert r1.strand in ('+', '-')
    if r1.strand != r2.strand:
        raise Exception("ERROR! query is {0} strand but ref is {1}!".format(r1.strand, r2.strand))

    if r1.strand == '+':
        five_prime_overhang = determine_extent_of_upstream_overhang(r1, r2)
        three_prime_overhang = determine_extent_of_downstream_overhang(r1, r2)
    else: # - strand
        three_prime_overhang = determine_extent_of_upstream_overhang(r1, r2)
        five_prime_overhang = determine_extent_of_downstream_overhang(r1, r2)

    return five_prime_overhang, three_prime_overhang

### end ###






###############################
##### start of sqanti run #####
###############################

if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument("isoform_gff", help="Input isoform GFF3, for transcript")
    parser.add_argument("cds_isoform_gff", help="Input isoform GFF3, for CDS")
    parser.add_argument("orf_tsv", help="Predicted ORF tsv")
    parser.add_argument("annotation_gtf", help="Annotation GTF, for transcript")
    parser.add_argument("cds_annotation_gtf", help="Annotation GTF, for CDS")
    parser.add_argument("-d", "--output_dir", default="output", help="Output directory (default: output)")
    parser.add_argument("-p", "--output_prefix", default="out", help="Output prefix (default: out")


    args = parser.parse_args()

    output_dir = args.output_dir
    output_prefix = args.output_prefix
    output_filename = os.path.join(output_dir, output_prefix+'.sqanti_protein_classification.tsv')

    if os.path.exists(output_dir):
        if not os.path.isdir(output_dir):
            print("ERROR: {output_dir} is not a directory! Abort!", file=sys.stderr)
            sys.exit(-1)
    else:
        os.makedirs(output_dir)

    SQANTIArgs = namedtuple('SQANTIArgs', 'isoform annotation dir output_prefix genename min_ref_len is_fusion corrGTF orf_tsv coverage window novel_gene_prefix')
    # NOTE - liz - i need to stick with these names since they are originally in sqanti input
    # for now, not changing into *_filename
    cds_isoform_gff = args.isoform_gff
    cds_annotation_gtf = args.annotation_gtf

    sqanti_args = SQANTIArgs(cds_isoform_gff,
                      cds_annotation_gtf,
                      dir=output_dir,
                      output_prefix=output_prefix,
                      corrGTF=args.isoform_gff,
                      orf_tsv=args.orf_tsv,
                      genename=None,
                      min_ref_len=0,
                      is_fusion=False,
                      coverage=None,
                      window=None,
                      novel_gene_prefix=None)


    #### process exon-based comparisons ####

    ## parse reference transcripts(GTF) to dicts
    refs_1exon_by_chr, refs_exons_by_chr, junctions_by_chr, junctions_by_gene, start_ends_by_gene, refDict = \
        reference_parser(sqanti_args, genome_chroms=None)

    ## parse query isoforms
    isoforms_by_chr, queryDict = protein_isoforms_parser(sqanti_args)

    ## read in orf calls from cpat (from lrp pipeline) into sqanti orfDict format
    orfDict = read_in_custom_orf_calls_into_orfDict(sqanti_args.orf_tsv)

    ## transcript isoform classification
    isoforms_info = isoformClassification(sqanti_args,
                                          isoforms_by_chr,
                                          refs_1exon_by_chr,
                                          refs_exons_by_chr,
                                          junctions_by_chr,
                                          junctions_by_gene,
                                          start_ends_by_gene,
                                          orfDict=orfDict,
                                          genome_dict=None,
                                          indelsJunc=None)


    #### process cds-based comparisons ####


    # updated named tuple to point to cds files (genocode, pacbio)
    ProteinArgs = namedtuple('ProteinArgs', 'isoform annotation dir output_prefix genename min_ref_len is_fusion corrGTF orf_tsv coverage window novel_gene_prefix')
    # NOTE - liz - i need to stick with these names since they are originally in sqanti input
    cds_isoform_gff = args.cds_isoform_gff #os.path.abspath(ddir + 'jurkat_cds_chr22.gff')
    cds_annotation_gtf = args.cds_annotation_gtf #os.path.abspath(ddir + 'gencode_cds_chr22.gtf')
    protein_args = ProteinArgs(cds_isoform_gff,
                               cds_annotation_gtf,
                               dir=output_dir,
                               output_prefix=output_prefix,
                               corrGTF=cds_isoform_gff,
                               orf_tsv=args.orf_tsv,
                               genename=None,
                               min_ref_len=0,
                               is_fusion=False,
                               coverage=None,
                               window=None,
                               novel_gene_prefix=None)

    ## parse reference transcripts(GTF) to dicts
    protein_refs_1exon_by_chr, protein_refs_exons_by_chr, protein_junctions_by_chr, \
    protein_junctions_by_gene, protein_start_ends_by_gene, protein_refDict_cds = \
        reference_parser(protein_args, genome_chroms=None)

    ## parse query isoforms
    protein_isoforms_by_chr, queryDict_cds = protein_isoforms_parser(protein_args)

    # isoform classification
    # note - orfDict is input, but results not in use for cds compare
    isoforms_info_cds = isoformClassification(protein_args,
                                              protein_isoforms_by_chr,
                                              protein_refs_1exon_by_chr,
                                              protein_refs_exons_by_chr,
                                              protein_junctions_by_chr,
                                              protein_junctions_by_gene,
                                              protein_start_ends_by_gene,
                                              genome_dict=None,
                                              indelsJunc=None,
                                              orfDict=None)

    #### write out results
    #### note - need to get "perfect subset" data while write-out
    f = open(output_filename, 'w')
    FIELDNAMES = ['pb', 'tx_cat', 'pr_splice_cat', 'tx_subcat', 'pr_splice_subcat',
                  'tx_tss_diff', 'tx_tts_diff', 'tx_tss_gene_diff', 'tx_tts_gene_diff',
                  'pr_nterm_diff', 'pr_cterm_diff', 'pr_nterm_gene_diff', 'pr_cterm_gene_diff',
                  'tx_transcripts', 'pr_transcripts',
                  'tx_gene', 'pr_gene',
                  'tx_num_exons', 'pr_num_exons',
                  'is_nmd',
                  'num_junc_after_stop_codon', 'num_nt_after_stop_codon',
                  'tx_5hang', 'tx_3hang',
                  'pr_nhang', 'pr_chang']
    writer = DictWriter(f, delimiter='\t', fieldnames=FIELDNAMES)
    writer.writeheader()
    for pb, pr in isoforms_info_cds.items():
        # pr is the protein object
        tx = isoforms_info[pb] # get transcript object

        # get perfect subset info for transcript (exon)
        transcript_match_id = tx.transcripts[0] #only get first transcript for overhang calc
        if transcript_match_id in refDict and pb in queryDict:
            ref = refDict[transcript_match_id]
            query = queryDict[pb]
            tx_5hang, tx_3hang = get_perfect_subset_status(query, ref)
        else:
            tx_5hang, tx_3hang = None, None

        # get perfect subset info for protein (cds)
        protein_match_id = pr.transcripts[0] #only get first transcript for overhang calc
        if protein_match_id in protein_refDict_cds and pb in queryDict_cds:
            pr_ref = protein_refDict_cds[protein_match_id]
            pr_query = queryDict_cds[pb]
            pr_5hang, pr_3hang = get_perfect_subset_status(pr_query, pr_ref)
        else:
            pr_5hang, pr_3hang = None, None

        info = {'pb': pb,
                'tx_cat': tx.str_class,
                'pr_splice_cat': pr.str_class,
                'tx_subcat': tx.subtype,
                'pr_splice_subcat': pr.subtype,
                'tx_tss_diff': tx.tss_diff,
                'tx_tts_diff': tx.tts_diff,
                'tx_tss_gene_diff': tx.tss_gene_diff,
                'tx_tts_gene_diff': tx.tts_gene_diff,
                'pr_nterm_diff': pr.tss_diff,
                'pr_cterm_diff': pr.tts_diff,
                'pr_nterm_gene_diff': pr.tss_gene_diff,
                'pr_cterm_gene_diff': pr.tts_gene_diff,
                'tx_transcripts': ','.join(tx.transcripts),
                'pr_transcripts': ','.join(pr.transcripts),
                'tx_gene': ','.join(tx.genes),
                'pr_gene': ','.join(pr.genes),
                'tx_num_exons': tx.num_exons,
                'pr_num_exons': pr.num_exons,
                'is_nmd': tx.is_NMD * 1,
                'num_junc_after_stop_codon': tx.num_junc_after_stop,
                'num_nt_after_stop_codon': orfDict[pb].num_3utr_nt,
                'tx_5hang': tx_5hang,
                'tx_3hang': tx_3hang,
                'pr_nhang': pr_5hang,
                'pr_chang': pr_3hang
                }
        writer.writerow(info)
    f.close()
    print(f"Output written to: {output_filename}")


