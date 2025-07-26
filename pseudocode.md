pseudocode

the problem... 
we have multiplexed reads, or many samples pooled together in one run with unique barcodes. we need to track how many have matched barcodes, index-hopped or unknown into separate files

informative output... 
counts of reads that have matched, index-hopped and unknown barcodes




initialize three variables, one for all matched pairs, one for index-hopped and one for unknown reads
initiate a dictionary to keep track of index 1 and 2 matched counts
initialize unknown counts and index hopped variables
Open all four zipped fastq files to read
    start the iteraions through the files and exctact the lines with the indexes and quility scores
        check the index1 and 2 records
            if the indexes have an N or do not meet the quality score requiremnt, 
                add the sequence of the indexes to header
                write the records from the files into their respective unknown files
                add one to the count of unknown counts
            else if the index 1 and 2 are a match... reverse complement
                add the sequence of the indexes to header
                write the records from the files into their respective matched files
                add one to the count of matched counts
            else the remaining files
                add the sequence of the indexes to header
                write the records from the files into their respective index hopped files
                add one to the count of index hopped counts
print the sums of the matched, index hopped and unknown pairs



high level functions...
def reverse_complement(seq: str) -> str:
    '''takes a string and returns reverse complement'''
    return reverse_complement
input AAACCC
output GGGTTT

def covert_phred(letter: str) -> int:
    '''takes ASCII character for Phred+33 and returns as an int'''
    return letter
input ?
output 30























