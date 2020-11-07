def orf_mapping(orf_coord, gencode, sample_gtf, seq_seq):
    exons = sample_gtf[sample_gtf['feature'] == 'exon']
    exons['exon_length'] = abs(exons['end'] - exons['start']) + 1

    exons = exons[['transcript_id', 'exon_id', 'exon_number', 'start', 'end', 'strand']]
    exons.rename(columns = {'start' : 'exon_start', 'end': 'exon_end'}, inplace = True)

    start_codons = gencode[gencode['feature'] == 'start_codon']
    start_codons = start_codons[['seqname','transcript_id','strand',  'start', 'end']]
    start_codons.rename(columns = {'start' : 'start_codon_start', 'end': 'start_codon_end'}, inplace = True)

    plus = plus_mapping(eons, orf_coord)
    minus = minus_mapping(exons, orf_coord)
    all_cds = pd.concat([plus, minus])

    def get_num_upstream_atgs(row):
        orf_start = int(row['orf_start'])
        acc = row['pb_acc']
        seq = orfs[acc] # get orf seq
        upstream_seq = seq[0:orf_start-1] # sequence up to the predicted ATG
        num_atgs = upstream_seq.count('ATG')
        return num_atgs

    all_cds['upstream_atgs'] = all_cds.apply(get_num_upstream_atgs, axis=1)


def plus_mapping(exons,orf_coord, start_codons):
    """
    Map minus strand
    """
    def compare_start_plus(row, start_codons):
        start = int(row['cds_start'])
        match = start_codons[start_codons['start_codon_start'] == start]
        match = match[match['seqname'].str.strip() == row['seqname'].strip()]
        if len(match) > 0:
            return list(match['transcript_id'])
        return None

    plus_exons = exons[exons['strand'] == '+']
    plus_exons['current_size'] = plus_exons.sort_values(by = ['transcript_id', 'exon_start']).groupby('transcript_id')['exon_length'].cumsum()
    plus_exons['prior_size'] = plus_exons['current_size'] - plus_exons['exon_length']

    plus_comb = pd.merge(orf_coord, plus_exons, left_on = 'pb_acc', right_on = 'transcript_id', how = 'inner')
    plus_comb = plus_comb[(plus_comb['prior_size'] <= plus_comb['orf_start']) &  (plus_comb['orf_start'] <= plus_comb['current_size'])]

    plus_comb['start_diff'] = plus_comb['orf_start'] - plus_comb['prior_size']
    plus_comb['cds_start'] = plus_comb['exon_start'] + plus_comb['start_diff'] - 1

    plus_comb['gencode_atg'] = plus_comb.apply(lambda row : compare_start_plus(row, start_codons), axis = 1)
    plus_comb.drop(columns=['exon_length', 'current_size', 'prior_size', 'start_diff'], inplace = True)
    return plus_comb
    

def minus_mapping(exons, orf_coord, start_codons):
    def compare_start_minus(row, start_codons):
        start = int(row['cds_start'])
        match = start_codons[(start_codons['start_codon_end'] == start) ]
        match = match[match['seqname'].str.strip() == row['seqname'].strip()]
        if len(match) > 0:
            return list(match['transcript_id'])
        return None
    
    minus_exons = exons[exons['strand'] == '-']
    minus_exons['current_size'] = minus_exons.sort_values(by = ['transcript_id', 'exon_start'], ascending=[True,False]).groupby('transcript_id')['exon_length'].cumsum()
    minus_exons['prior_size'] = minus_exons['current_size'] - minus_exons['exon_length']

    minus_comb = pd.merge(orf_coord, minus_exons, left_on = 'pb_acc', right_on = 'transcript_id', how = 'inner')

    minus_comb = minus_comb[minus_comb['orf_start'].between(minus_comb['prior_size'],minus_comb['current_size'])]

    minus_comb['start_diff'] = minus_comb['orf_start'] - minus_comb['prior_size']
    minus_comb['cds_start'] = minus_comb['exon_end'] - minus_comb['start_diff'] + 1
    minus_comb.drop(columns=['exon_length', 'current_size', 'prior_size', 'start_diff'], inplace = True)

    minus_comb['gencode_atg'] = minus_comb.apply(lambda row : compare_start_minus(row, start_codons), axis = 1)

def read_orf(filename):
    """
    Reads the ORF file and formats the column names
    Keep only predictions with score higher than protein-coding threshold

    Parameters
    ---------
    filename : str
        location of orf coordinate file

    Returns
    --------
    orf: pandas DataFrame
    """
    orf = pd.read_csv(filename, sep = '\t')
    orf[['pb_acc', 'misc', 'orf']] = orf['ID'].str.split('_', expand=True)
    orf = orf.drop(labels=['ID', 'misc'], axis=1)
    orf.columns = ['len', 'orf_strand', 'orf_frame', 'orf_start', 'orf_end', 'orf_len', 'fickett', 'hexamer', 'coding_score', 'pb_acc', 'orf_rank',]
    return orf



def orf_calling(orf):
    def call_orf(group):

        score_threshold = 0.364

        good_score = group[group['coding_score'] >= score_threshold]
        if(len(good_score) >= 1):
            group = good_score

            
        highscore = 0.9
        score_threshold = 0.364

        group['orf_calling_confidence'] = 'None'
        group['atg_rank'] = group['upstream_atgs'].rank(ascending=True)
        group['score_rank'] = group['coding_score'].rank(ascending=False)

        group['atg_score'] = group['upstream_atgs'].apply(lambda x : 1/x  if x > 1 else 0.99)
        group['gencode_score'] = group['gencode_atg'].apply(lambda x : 0 if x == '' else 0.8)
        group['orf_score'] = group.apply(lambda row: 1 - (1-row['coding_score']*0.99)*(1-row['atg_score'])*(1-row['gencode_score']), axis = 1)

        def calling_confidence(row):
            if row['atg_rank'] == 1 and row['score_rank'] == 1:
                return 'Clear Best ORF'
            elif row['coding_score'] > highscore and row['atg_rank'] == 1:
                return 'Nonsense Mediated Decay'
            elif row['coding_score'] < score_threshold:
                return 'Low Quality ORF'
            return 'None'

        group['orf_calling_confidence'] = group.apply(lambda row : calling_confidence(row), axis = 1)
        group = group.sort_values(by='atg_rank').reset_index(drop=True)
        if group.loc[0,'atg_rank'] == group.loc[0,'score_rank']:
            return group.head(1)
            

        
        group = group.sort_values(by='orf_score', ascending=False).reset_index(drop=True)
        return group.head(1)
        
    called_orf = orf.groupby('pb_acc').apply(newbest_filter).reset_index(drop=True)
    return called_orf
    
    