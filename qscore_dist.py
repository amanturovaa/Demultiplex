#!/usr/bin/env python

import argparse #import everything to be able to argparse, plot and read .gz file
import matplotlib.pyplot as plt
import gzip
#from bioinfo import convert_phred  
def convert_phred(letter): #made the convert-phred function
    """Converts ASCII character into Phred score"""
    return ord(letter) - 33


def process_fastq(file): 
    qsums = [0] * 101 #initialize my running total of the quality scores at each position (101 total inclusive, exclusive)
    qcounts = [0] * 101  #initilaize counts for how many quality scores are at each base position 
    with gzip.open(file, "rt") as fh: #add gzip for .gz 
        for i, line in enumerate(fh): #iterate over each line 
            if i % 4 == 3:  #grab the quality score
                line = line.strip() #remove new lines
                for pos, letter in enumerate(line[:101]): #loop over all the letters in the qscore to position 101
                    q = convert_phred(letter) #covert the letters to their phred score and store in q
                    qsums[pos] += q #adds the score to the running qsums total
                    qcounts[pos] += 1 #increments by 1
    means = [sum / count if count != 0 else 0 for sum, count in zip(qsums, qcounts)] #take sum and count and divide to get average and avoid zero division if that is possible 
    return means

def plot_quality(means, label):
    plt.figure()
    plt.bar(range(len(means)), means)
    plt.xlabel("Base position")
    plt.ylabel("Mean Quality Score")
    plt.title(f"Mean Quality Scores - {label}")
    plt.xlim(0, 101)  # scales plot to be in correct range
    plt.ylim(0, 42)   # scales plt to be within phred score
    plt.savefig(f"{label}_quality_distribution.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate mean quality score histograms per base.")
    parser.add_argument("-r1", "--read1", required=True)
    parser.add_argument("-r2", "--read2", required=True)
    parser.add_argument("-r3", "--read3", required=True)
    parser.add_argument("-r4", "--read4", required=True)
    args = parser.parse_args()

    for fname, label in zip([args.read1, args.read2, args.read3, args.read4],
                            ["read1", "read2", "read3", "read4"]):
        means = process_fastq(fname)
        plot_quality(means, label)
