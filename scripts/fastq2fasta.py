import os
import sys
from Bio import SeqIO

def fastq2fasta(input_fastq, result_dir):
    print('convert '+str(input_fastq)+' to fasta')
    SeqIO.convert(input_fastq, "fastq", result_dir+"/input.fasta", "fasta")

# def fastq2fasta():
#     argv = sys.argv
#     file = argv[1]
#     SeqIO.convert(file+".fastq", "fastq", file+".fasta", "fasta")

# if __name__ == "__main__":
#     fastq2fasta()

# from bioconvert import fasta2fastq

# def fasta_to_fastq(input_fasta, result_dir):
#     print('begin fasta2fastq')
#     # argv = sys.argv
#     # file_name = argv[1]
#     file_name = input_fasta
#     print(file_name)
#     if ".fa" == file_name[-3:]:
#         # transformer = fasta2fastq.Fasta2Fastq(file_name, file_name[:-3]+".fastq")
#         transformer = fasta2fastq.Fasta2Fastq(file_name, result_dir+"/resequencing/assembly_result.fastq")
#         transformer._method_pysam()
#     elif ".fasta" == file_name[-6:]:
#         # transformer = fasta2fastq.Fasta2Fastq(file_name, file_name[:-6]+".fastq")
#         transformer = fasta2fastq.Fasta2Fastq(file_name, result_dir+"/resequencing/assembly_result.fastq")
#         transformer._method_pysam()
#     print('end fasta2fastq')