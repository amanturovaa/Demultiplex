#!/usr/bin/env python3

import argparse
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="demultiplexer for paired-end reads")
    parser.add_argument("-r1", "--read1", help="r1 fastq.gz file (read 1)", required=True)
    parser.add_argument("-r2", "--read2", help="r2 fastq.gz file (index 1)", required=True)
    parser.add_argument("-r3", "--read3", help="r3 fastq.gz file (index 2)", required=True)
    parser.add_argument("-r4", "--read4", help="r4 fastq.gz file (read 2)", required=True)
    parser.add_argument("-x", "--indexes", help="index file with sequences", required=True)
    parser.add_argument("-q", "--quality", help="quality cutoff", type=int)
    parser.add_argument("-o", "--output", help="output directory prefix", default="")
    return parser.parse_args()

def reverse_complement(seq):
    """reverse complement"""
    comp_table = str.maketrans('ATGCN', 'TACGN')
    return seq.translate(comp_table)[::-1]

def avg_quality(qual_string):
    """average quality score"""
    total = 0
    for char in qual_string:
        total += ord(char) - 33
    return total / len(qual_string)

def main():
    args = get_args()
    
    #known indexes
    known_indexes = set()
    with open(args.indexes, 'r') as f:
        next(f)  #skip header
        for line in f:
            index_seq = line.strip().split()[-1]  ##last column
            known_indexes.add(index_seq)
    
    #starting open output files
    output_files = {}
    prefix = args.output + "/" if args.output else ""
    
    for index in known_indexes:
        r1_file = open(f"{prefix}{index}_R1.fastq", 'w')
        r2_file = open(f"{prefix}{index}_R2.fastq", 'w')
        output_files[index] = (r1_file, r2_file)
    
    #output files for unknown and hopped for both reads
    unknown_r1 = open(f"{prefix}unknown_R1.fastq", 'w')
    unknown_r2 = open(f"{prefix}unknown_R2.fastq", 'w')
    hopped_r1 = open(f"{prefix}hopped_R1.fastq", 'w')
    hopped_r2 = open(f"{prefix}hopped_R2.fastq", 'w')
    
    #initialize counters
    matched_count = 0
    hopped_count = 0
    unknown_count = 0
    total = 0
    
    #start reading files
    with gzip.open(args.read1, 'rt') as r1, \
         gzip.open(args.read2, 'rt') as i1, \
         gzip.open(args.read3, 'rt') as i2, \
         gzip.open(args.read4, 'rt') as r2:
        
        while True:
            # Read records
            r1_header = r1.readline().strip()
            if not r1_header:  #for the end of file
                break
                
            r1_seq = r1.readline().strip()
            r1_plus = r1.readline().strip()
            r1_qual = r1.readline().strip()
            
            i1_header = i1.readline().strip()
            i1_seq = i1.readline().strip()
            i1_plus = i1.readline().strip()
            i1_qual = i1.readline().strip()
            
            i2_header = i2.readline().strip()
            i2_seq = i2.readline().strip()
            i2_plus = i2.readline().strip()
            i2_qual = i2.readline().strip()
            
            r2_header = r2.readline().strip()
            r2_seq = r2.readline().strip()
            r2_plus = r2.readline().strip()
            r2_qual = r2.readline().strip()
            
            total += 1
            
            #get progress every 1 million reads to make sure it's working
            if total % 1000000 == 0:
                print(f"{total:,} reads")
            
            #reverse complement of index 2
            i2_rc = reverse_complement(i2_seq)
            
            #process quality
            i1_avg_qual = avg_quality(i1_qual)
            i2_avg_qual = avg_quality(i2_qual)
            
            #index pair string
            index_pair = f"{i1_seq}-{i2_rc}"
            
            #add indexes to headers
            r1_header_new = f"{r1_header} {index_pair}"
            r2_header_new = f"{r2_header} {index_pair}"
            
            # Create record strings once
            r1_record = f"{r1_header_new}\n{r1_seq}\n{r1_plus}\n{r1_qual}\n"
            r2_record = f"{r2_header_new}\n{r2_seq}\n{r2_plus}\n{r2_qual}\n"
            
            #where to write
            if (i1_seq not in known_indexes or 
                i2_rc not in known_indexes or 
                i1_avg_qual < args.quality or 
                i2_avg_qual < args.quality):
                #unknowns
                unknown_r1.write(r1_record)
                unknown_r2.write(r2_record)
                unknown_count += 1
                
            elif i1_seq == i2_rc:
                #match
                output_files[i1_seq][0].write(r1_record)
                output_files[i1_seq][1].write(r2_record)
                matched_count += 1
                
            else:
                #index hopped
                hopped_r1.write(r1_record)
                hopped_r2.write(r2_record)
                hopped_count += 1
    
    #close all files
    for r1_file, r2_file in output_files.values():
        r1_file.close()
        r2_file.close()
    
    unknown_r1.close()
    unknown_r2.close()
    hopped_r1.close()
    hopped_r2.close()
    
    #print out summary results
    total = matched_count + hopped_count + unknown_count
    print(f"Total reads processed: {total:,}")
    print(f"Matched pairs: {matched_count:,} ({matched_count/total*100:.2f}%)")
    print(f"Index hopped: {hopped_count:,} ({hopped_count/total*100:.2f}%)")
    print(f"Unknown/low quality: {unknown_count:,} ({unknown_count/total*100:.2f}%)")

if __name__ == "__main__":
    main()