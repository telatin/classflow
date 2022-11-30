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
import networkx as nx

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

def getId(taxon):
    return taxon[0]
def getRank(id, db):
    if int(id) == 0:
        return "U"
    elif int(id) == 1:
        return "R"
    for taxid in db:
        for tuple in db[taxid]:
            if int(tuple[0]) == int(id):
                if tuple[1] == "no rank":
                    return "-"
                elif tuple[1] == "superkingdom":
                    return "D"
                else:
                    return tuple[1][0].upper()

def getName(id, db):
    if int(id) == 0:
        return "unclassified"
    elif int(id) == 1:
        return "root"
    for taxid in db:
        for tuple in db[taxid]:
            if int(tuple[0]) == int(id):
                return '  ' * len(tuple) + tuple[2]

def getUpperRanks(id, db):
    upperRanksIds = []
    for taxon, rank, label in db[id]:
        upperRanksIds.append(taxon)
    return upperRanksIds

def label(query, db):

    for taxon, rank, label in db[query]:
        return label
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse the TSV output from Kraken2 and generate its taxonomy report')
    parser.add_argument('TSV', help='TSV file from Kraken2')
    parser.add_argument('-o', '--output', help='Output report file [default: stdout]')
    parser.add_argument('-b', '--basename', type=str, default="temp_kraken_", help='Prefix [default: %(default)s]')
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


    ## PANDAS
    ## Load TSV file or part of it

    colNames =['status', 'read_name', 'taxon', 'read_len', 'classification']

    if args.max_reads > 0:
        tsv_file = pandas.read_csv(args.TSV, sep='\t', header=None, names=colNames, nrows=args.max_reads)
    else:
        tsv_file = pandas.read_csv(args.TSV, sep='\t', header=None, names=colNames)

    total_records = len(tsv_file)
    total_unclassified = len(tsv_file[tsv_file['taxon'] == 0])
    classified_records = len(tsv_file[tsv_file['taxon'] != 0])

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
    original_taxa_ids = set(all_taxa_string.split())
    taxonomy_db = taxonkit_lineages(all_taxa_string)

    with open(args.basename + "taxonomy_db.txt", "w") as outfile: 
        for taxid, tuples in taxonomy_db.items():
            outfile.write(f"{taxid}\t{tuples}\n")

    
    # Add column "counts" defined as 1
    tsv_file['counts'] = 1
    tsv_file['up_counts'] = 1
    

    # Collapse rows by taxon, summing counts in "counts" column
    tsv_file = tsv_file.groupby('taxon').sum(numeric_only = True)
    
    
    tsv_file['up_counts'] = 0
    
    tsv_file.to_csv(args.basename + "1_start.tsv", sep='\t', header=True, index=True)
    
    # Drop column read_len
    tsv_file.drop(columns=['read_len'], inplace=True)
    tsv_file.sort_values(by=['counts'], ascending=False, inplace=True)

    
    # Extract all taxids from taxonomy_db
    upper_ranks_ids = list(taxonomy_db.keys())
    for id, tuples in taxonomy_db.items():
        for t in tuples:
            if t[0] not in upper_ranks_ids:
                upper_ranks_ids.append(t[0])

    
    # Add a row "1" with counts = classified
    #tsv_file.loc[1] = [classified_records, 0]
    # add to tsv_file all taxids from taxonomy_db that are not in tsv_file
    for taxid in upper_ranks_ids:
        if int(taxid) not in tsv_file.index:
            tsv_file.loc[int(taxid)] = [0, 0]
    
    tsv_file.to_csv(args.basename + "2_expand.tsv", sep='\t', header=True, index=True)
    # Now propagate each counts to all upper ranks in tsv_file

    for id in original_taxa_ids:
        intid = int(id)
        # Gather the value of "counts" for the row with id = id
        local_counts = tsv_file.loc[intid]['counts']

        # Gather its upper ranks
        local_upper_ranks = getUpperRanks(id, taxonomy_db) 
        local_upper_ranks = [int(x) for x in local_upper_ranks]

        # Add the "counts" to all upper ranks
        for taxid in local_upper_ranks:
            # increment the up_counts of the row with id = taxid
            if taxid in tsv_file.index and taxid not in original_taxa_ids:
                tsv_file.loc[taxid]['up_counts'] += int(local_counts)
            else:
                print(f"ERROR: taxid  {taxid} not in tsv_file.index {str(taxid) in original_taxa_ids} {str(taxid) in upper_ranks_ids}")
                pass

                #tsv_file.loc[taxid]['up_counts'] += int(local_counts)
    

    # Set up_counts = counts if taxon == 0
    tsv_file.loc[0]['up_counts'] = tsv_file.loc[0]['counts']
    # print tsv in descending order of counts
    tsv_file.sort_values(by=['up_counts'], ascending=False, inplace=True)
    
    # Add column "percentage" with the percentage of up_counts over the total of up_counts
    tsv_file['percentage'] = tsv_file['up_counts']/total_records*100

    # Add colun "rank" with the rank of the taxid
    tsv_file['rank'] = tsv_file.index.map(lambda x: getRank(x, taxonomy_db))
    tsv_file['description'] = tsv_file.index.map(lambda x: getName(x, taxonomy_db))
    
 
    # Make the index an ordinary columnt "taxid"
    tsv_file.reset_index(level=0, inplace=True)
    
    # change order of columns PERC	up_counts	count	rank	taxid	string
    tsv_file = tsv_file[['percentage', 'up_counts', 'counts', 'rank', 'taxon', 'description']]
    # Round percentage to 2 decimals exactly (not rounding)
    tsv_file['percentage'] = tsv_file['percentage'].apply(lambda x: f"{x:.2f}")
    # Convert percentage to string with padding spaces 
    tsv_file['percentage'] = tsv_file['percentage'].astype(str).str.pad(6, side='left')
    
    tsv_file.to_csv(args.basename + "3_end.tsv", sep='\t', header=False, index=False)
   
    # Print the sum of counts and up_counts
    print(f"Total records: {total_records}", file=sys.stderr)
    print(f"Classified records: {classified_records}", file=sys.stderr)
    print(f"Confident records: {confident_records}", file=sys.stderr)
    
    # Print the sum of counts and up_counts from tsv_file
    print(f"Sum of counts: {tsv_file['counts'].sum()}", file=sys.stderr)
    print(f"Sum of up_counts: {tsv_file['up_counts'].sum()}", file=sys.stderr)