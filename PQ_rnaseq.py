### Python script to interact with os for RNA seq trimming, alignment, and counting
### Author: Pengfei Qiao (pq26@cornell.edu)

import os
import sys

#Get fastq files
os.system("ls ./rawfastq/*.fastq > fastqfiles.txt") # Change path to fastq files
fastq = [l.strip() for l in open("fastqfiles.txt","rU")]
fastq = sorted(fastq, key=lambda x: x.split("/")[2].split("_")[8]) #Change if fastq file names are differnet
os.system("rm fastqfiles.txt")
os.system("mkdir rawcounts")
os.system("mkdir mergedcounts")

#Index genome
os.system("hisat2-build %s genome" %(sys.argv[1])) #sys.argv[1] is the genome file

# Define function to Trim, Align, Sort and Count (TASC)
# Change parameters if desired
def tasc(fastqinput):
  trimmedoutput = fastqinput[53:58] + ".fastq"
  samoutput = fastqinput[53:58]+".sam" # Index - Need to look at file names and fill, same for the rest in this script
  bamoutput = fastqinput[53:58]+".bam"
  sortedbamoutput = fastqinput[53:58]+"_sorted.bam"
  txtoutput = "./rawcounts/"+fastqinput[53:58]+".txt"
  os.system("java -jar trimmomatic-0.36.jar SE -threads 12 -phred33 %s %s HEADCROP:12 ILLUMINACLIP:./trial/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36" %(fastqinput,trimmedoutput))
  os.system("hisat2 -p 6 -t -x genome -U %s -S %s" %(trimmedoutput,samoutput))
  os.system("samtools view -bS %s > %s" %(samoutput, bamoutput))
  os.system("rm %s" %(trimmedoutput))
  os.system("rm %s" %(samoutput))
  os.system("samtools sort %s -o %s" %(bamoutput, sortedbamoutput))
  os.system("htseq-count -i Name -f bam -r name --stranded=no %s %s > %s" %(sortedbamoutput, sys.argv[2], txtoutput)) #sys.argv[2] is the gff file
  return txtoutput

for i in fastq:
  tasc(i)
  
# Merge counts into one file
countfiles = ["./rawcounts/%s" %(l.strip()) for l in os.listdir("./rawcounts") if "PQ" in l]
mergedcounts = dict([(l.strip().split("\t")[0],[l.strip().split("\t")[1]]) for l in open(countfiles[0],"rU")])
for i in countfiles[1:]:
  infile = open(i,"rU")
  for l in infile:
    mergedcounts[l.strip().split("\t")[0]].append(l.strip().split("\t")[1])
  infile.close()

# Write out merged file:
outfile = open("rawcounts.txt","w")
outfile.write("Gene\t")
outfile.write("\t".join(countfiles))
outfile.write("\n")
for i in mergedcounts.keys():
  counts = mergedcounts[i]
  outfile.write("%s\t" %(i))
  outfile.write("\t".join(counts))
  outfile.write("\n")
outfile.close()

