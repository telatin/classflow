#!/usr/bin/env python
"""
Parse the TSV output from Kraken2 and generate its taxonomy report
"""
# support lru_cache
from functools import lru_cache
import argparse
import sys
import os
import subprocess
import pandas
import functools
def kraken2_tax_confidence(string, min_confidence=0):
    """
    Return the confidence of a Kraken2 classification
    2755405:47 0:2 2755405:2 0:10 2755405:2 0:14 2755405:31 816:1 2755405:7 |:| 2755405:2 816:11 2755405:3 816:2 2755405:17 0:32 816:5 2755405:5 816:3 2755405:36
    """
    total = 0
    taxa_counts = {}
    for part in string.split():
        
        try:
            taxon, count = part.split(':')
        except ValueError:
            continue
        
        if taxon == '|':
            continue
        total += int(count)
        if taxon not in taxa_counts:
            taxa_counts[taxon] = int(count)
        else:
            taxa_counts[taxon] += int(count)
    
    top_taxon = max(taxa_counts, key=taxa_counts.get)
    top_count = taxa_counts[top_taxon]
    confidence = top_count / total

    if confidence < min_confidence:
        return 0, confidence
    return top_taxon, confidence
def kraken2_confidence(string, min_confidence=0):
    """
    Return the confidence of a Kraken2 classification
    2755405:47 0:2 2755405:2 0:10 2755405:2 0:14 2755405:31 816:1 2755405:7 |:| 2755405:2 816:11 2755405:3 816:2 2755405:17 0:32 816:5 2755405:5 816:3 2755405:36
    """
    total = 0
    taxa_counts = {}
    for part in string.split():
        
        try:
            taxon, count = part.split(':')
        except ValueError:
            continue
        
        if taxon == '|':
            continue
        total += int(count)
        if taxon not in taxa_counts:
            taxa_counts[taxon] = int(count)
        else:
            taxa_counts[taxon] += int(count)
    
    top_taxon = max(taxa_counts, key=taxa_counts.get)
    top_count = taxa_counts[top_taxon]
    confidence = top_count / total

    return confidence

 
def taxonkit_lineages(ids):
    db = {"0": [(0, 'no rank', 'unclassified')]}
    stdin = str(ids)
    command = ['taxonkit', 'lineage', '-Rtn']
    proc = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    stdout, stderr = proc.communicate(stdin.encode())
    for line in stdout.decode().splitlines():
        input_id, raw_taxa, raw_ids, taxon_label, raw_ranks = line.strip().split('\t')
        taxa = raw_taxa.split(';')
        ids = raw_ids.split(';')
        ranks = raw_ranks.split(';')
        db[str(input_id)] = list(zip(ids, ranks, taxa))
    
    return db

def taxon_counts(zipped_taxa):
    counts = {}
    for taxon, rank, label in zipped_taxa:
        if taxon not in counts:
            counts[taxon] = 1
        else:
            counts[taxon] += 1
    return counts


def label(query, db):

    for taxon, rank, label in db[query]:
        return label
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse the TSV output from Kraken2 and generate its taxonomy report')
    parser.add_argument('TSV', help='TSV file from Kraken2')
    parser.add_argument('-o', '--output', help='Output report file [default: stdout]')
    parser.add_argument('-m', '--min-confidence', type=float, default=0.0, help='Minimum confidence value to report')
    parser.add_argument('--min-support', type=int, default=0, help='Minimum number of reads to report')
    parser.add_argument('-x', '--max-reads', type=int, default=0, help='Stop after this many reads')
    
    parser.add_argument('--db', help='Kraken2 database directory')
    args = parser.parse_args()


    
    if not os.path.exists(args.TSV):
        print('Error: TSV file does not exist', file=sys.stderr)
        sys.exit(1)
    
    if args.output is None:
        output = sys.stdout
    else:
        output = open(args.output, 'w')
    
    # 1: Read the TSV file
    counts = {}
    c = 0

    if args.max_reads > 0:
        tsv_file = pandas.read_csv(args.TSV, sep='\t', header=None, names=['status', 'read_name', 'taxon', 'read_len', 'classification'], nrows=args.max_reads)
    else:
        tsv_file = pandas.read_csv(args.TSV, sep='\t', header=None, names=['status', 'read_name', 'taxon', 'read_len', 'classification'])

    total_records = len(tsv_file)
    total_unclassified = len(tsv_file[tsv_file['taxon'] == 0])

    # Delete records where 'taxon' is 0 in place
    tsv_file = tsv_file.loc[tsv_file["taxon"] != 0]
    classified_records = len(tsv_file)

    # Add column "confidence" defined as kraken2_confidence(classification)
    tsv_file['confidence'] = tsv_file['classification'].apply(kraken2_confidence)


    # Delete records where 'confidence' is less than args.min_confidence in place
    tsv_file = tsv_file.loc[tsv_file["confidence"] >= args.min_confidence]
    confident_records = len(tsv_file)

    # Delete column "confidence" in place
    tsv_file.drop(columns=['confidence'], inplace=True)

    print(f"% classified: {classified_records/total_records*100:.2f}%", file=sys.stderr)
   

    
    # Build taxonomy cache: join all tsv_file "taxon" columns except for "0"
    all_taxa_string = '\n'.join(tsv_file[tsv_file['taxon'] != 0]['taxon'].astype(str))
    taxonomy_db = taxonkit_lineages(all_taxa_string)

    # Add column "counts" defined as 1
    tsv_file['counts'] = 1
   

    # Collapse rows by taxon, summing counts in "counts" column
    tsv_file = tsv_file.groupby('taxon').sum()
    

    for row in tsv_file.itertuples():
        parent_taxa = taxonomy_db[str(row.Index)]
        for taxon, rank, label in parent_taxa:
            # Add row to tsv_file if taxon not in tsv_file, else add counts to existing row
            if taxon not in tsv_file.index:
                tsv_file.loc[taxon] = [0]
            else:
                # Increase counts by row.counts
                tsv_file.loc[taxon, 'counts'] += row.counts
        
       
    # print tsv in descending order of counts
    tsv_file.sort_values(by=['counts'], ascending=False, inplace=True)
    
    print(tsv_file.head(10))