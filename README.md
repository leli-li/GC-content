# Introduction

### This program is a tool that uses sliding window to calculate GC content and plot graphs.
### The software first reads the sequence from the fasta file, then merges the sequences, and calculates the GC content of each window according to the sliding window and step size set by the user.

### Calculating GC content: (G+C)/(A+T+G+C)

# Input file
### The input file should be fasta format. The file can contain one or more sequences.


# Program function
## Single file processing :
### combine sequences in fasta file and plot GC content through this large sequence. The output is a line chart.
### plot distribution histogram based on each sequence in fasta file. The output is a histogram.
### filter sequence using a range of GC content. The output is a new fasta file.

## pairwise comparison : 
### Input two files and plot two lines of GC conteng in one graph. The output is one line chart. 

## File batch processing :
### Input more than one fasta files and do single file processing on each file. The number of output line charts is the same as the number of input files.
### The input files and this script should be put in the same directory without other files.

## Batch pairwise comparison :
### Input more than two fasta files and pairwise comparison on these files. If the number of input files is n, the number of output files should be [n*(n-1)]/2.
### The input files and this script should be put in the same directory without other files.

## all lines plotted in one graph :
### Plot all lines in one line chart. The output is one graph.
### The input files and this script should be put in the same directory without other files.


# Usage
python gc_content.py [-h] [-f FILE [FILE ...]] [-w WINDOW] [-s STEP]
                      [-r Lower limit upper limit]
                      {s,c,b,bc,a}

## required arguments:
### This parameter can only choose one from {s,c,b,bc,a}
### single(s)
### pairwise comparison(c)
### batch single file processing(b)
### batch pairwise comparison(bc) 
### all lines plotted in one graph(a)

## optional arguments:
###   -h, --help show this help message and exit
###   -f FILE [FILE ...] --file FILE [FILE ...] : input file
### 			file name should be the whole name (with suffix)
###   -w WINDOW, --window WINDOW : set the size of sliding window
###   -s STEP, --step STEP : step size
###   -r Lower limit upper limit, --range Lower limit upper limit : set the range of GC content to filter sequences. between 0 and 100

# example
### single file:
python gc_content.py s -f input_file -w 100000 -s 1000
### fileter sequence:
python gc_content.py s -f input_file -r 10 30


### pairwise comparison:
python gc_content.py c -f input_file1 input_file2 -w 10000 -s 100

### batch single file processing:
python gc_content.py b -w 100000 -s 1000

### batch pairwise comparison:
python gc_content.py bc -w 100000 -s 1000

### plot all lines in one graph:
python gc_content.py a -w 100000 -s 1000



