/*
 * Run pairwise blast of isoforms between protein databases
 */

process mapisoforms{
    tag "mapisoforms: $xxx"
    
    input:
    set val(name), file(db) from proteindbs 
    
    output:
    set val(name), file("${name}.tsv) into blast_result
    
    script:
    """
    blastp -db ${db} -query ${querydb} -outfmt "6 qseqid qlen qstart qend sseqid slen sstart send evalue bitscore length nident gaps" -out ${blast_result}
    """
}

