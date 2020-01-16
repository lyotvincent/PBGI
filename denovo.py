from Bio.Blast.Applications import NcbiblastnCommandline
import subprocess
import os

import assembly
import scripts.handle_denovo_blast_result
import scripts.prokka_annotation
import scripts.prodigal_prediction

class Denovo:

    result_dir = None
    conf = None
    input_file = None
    input_file_1 = None
    input_file_2 = None
    input_file_12 = None

    def __init__(self, result_dir, conf, input_file=None, input_file_1=None, input_file_2=None, input_file_12=None):
        self.result_dir = result_dir
        self.conf = conf
        self.input_file = input_file
        self.input_file_1 = input_file_1
        self.input_file_2 = input_file_2
        self.input_file_12 = input_file_12

    def run_single(self):
        print('denovo single')
        print('begin assembly')
        assembly_conf = self.conf['denovo']['assembly']
        os.mkdir(os.path.abspath('.')+'/'+self.result_dir+'/denovo/assembly')
        assembly_obj = assembly.Assembly(self.result_dir+'/denovo/assembly', self.conf, input_file=self.input_file)
        if not assembly_conf['megahit']['enable'] == 0:
            assembly_obj.megahit_single()
            assembly_result = self.result_dir+'/denovo/assembly/megahit_out/final.contigs.fa'
        elif not assembly_conf['spades']['enable'] == 0:
            assembly_obj.spades_single()
            assembly_result = self.result_dir+'/denovo/assembly/spades_out/contigs.fasta'
        elif not assembly_conf['velvet']['enable'] == 0:
            assembly_obj.velvet_single()
            assembly_result = self.result_dir+'/denovo/assembly/velvet_output/contigs.fa'
        else:
            assembly_obj.megahit_single()
            assembly_result = self.result_dir+'/denovo/assembly/megahit_out/final.contigs.fa'
        print('end assembly')

        print("Begin QUAST")
        subprocess.run('quast.py '+assembly_result+' --min-contig 50 -o '+self.result_dir+'/denovo/quast_out', shell=True, check=True)
        print("End QUAST")

        self.blastn_and_annotation(assembly_result)

        print('end single')

    def run_paired(self):
        print('denovo paired')
        print('begin assembly')
        assembly_conf = self.conf['denovo']['assembly']
        os.mkdir(os.path.abspath('.')+'/'+self.result_dir+'/denovo/assembly')
        assembly_obj = assembly.Assembly(self.result_dir+'/denovo/assembly', self.conf, input_file_1=self.input_file_1, input_file_2=self.input_file_2)
        if not assembly_conf['megahit']['enable'] == 0:
            assembly_obj.megahit_paired()
            assembly_result = self.result_dir+'/denovo/assembly/megahit_out/final.contigs.fa'
        elif not assembly_conf['spades']['enable'] == 0:
            assembly_obj.spades_paired()
            assembly_result = self.result_dir+'/denovo/assembly/spades_out/contigs.fasta'
        elif not assembly_conf['velvet']['enable'] == 0:
            assembly_obj.velvet_paired()
            assembly_result = self.result_dir+'/denovo/assembly/velvet_output/contigs.fa'
        else:
            assembly_obj.megahit_paired()
            assembly_result = self.result_dir+'/denovo/assembly/megahit_out/final.contigs.fa'
        print('end assembly')

        print("Begin QUAST")
        subprocess.run('quast.py '+assembly_result+' --min-contig 50 -o '+self.result_dir+'/denovo/quast_out', shell=True, check=True)
        print("End QUAST")

        self.blastn_and_annotation(assembly_result)

        print('end paired')

    def run_interleaved(self):
        print('denovo interleaved')
        assembly_conf = self.conf['denovo']['assembly']
        os.mkdir(os.path.abspath('.')+'/'+self.result_dir+'/denovo/assembly')
        assembly_obj = assembly.Assembly(self.result_dir+'/denovo/assembly', self.conf, input_file_12=self.input_file_12)
        if not assembly_conf['megahit']['enable'] == 0:
            assembly_obj.megahit_interleaved()
            assembly_result = self.result_dir+'/denovo/assembly/megahit_out/final.contigs.fa'
        elif not assembly_conf['spades']['enable'] == 0:
            assembly_obj.spades_interlaced()
            assembly_result = self.result_dir+'/denovo/assembly/spades_out/contigs.fasta'
        elif not assembly_conf['velvet']['enable'] == 0:
            assembly_obj.velvet_interleaved()
            assembly_result = self.result_dir+'/denovo/assembly/velvet_output/contigs.fa'
        else:
            assembly_obj.megahit_interleaved()
            assembly_result = self.result_dir+'/denovo/assembly/megahit_out/final.contigs.fa'

        print("Begin QUAST")
        subprocess.run('quast.py '+assembly_result+' --min-contig 50 -o '+self.result_dir+'/denovo/quast_out', shell=True, check=True)
        print("End QUAST")

        self.blastn_and_annotation(assembly_result)

        print('end interleaved')

    def run_3gs(self):
        print('begin denovo 3gs')

        assembly_conf = self.conf['denovo']['assembly']
        if not assembly_conf['enable'] == 0:
            os.mkdir(os.path.abspath('.')+'/'+self.result_dir+'/denovo/assembly')
            assembly_obj = assembly.Assembly(self.result_dir+'/denovo/assembly', self.conf, input_file=self.input_file)
            if not assembly_conf['canu']['enable'] == 0:
                assembly_obj.canu()
                assembly_result = self.result_dir+'/denovo/assembly/canu_output/canu_assembly_result.contigs.fasta'
            elif not assembly_conf['spades']['enable'] == 0:
                assembly_obj.spades_3gs()
                assembly_result = self.result_dir+'/denovo/assembly/spades_out/contigs.fasta'
            else:
                assembly_obj.canu()
                assembly_result = self.result_dir+'/denovo/assembly/canu_output/canu_assembly_result.contigs.fasta'

            print("Begin QUAST")
            subprocess.run('quast.py '+assembly_result+' --min-contig 50 -o '+self.result_dir+'/denovo/quast_out', shell=True, check=True)
            print("End QUAST")
        else:
            scripts.fastq2fasta.fastq2fasta(self.input_file, self.result_dir+"/denovo")
            assembly_result = self.result_dir+"/denovo/input.fasta"

        self.blastn_and_annotation(assembly_result)

        print('end denovo 3gs')

    def blastn_and_annotation(self, assembly_result):
        if not self.conf['denovo']['ncbi_bacteria_blastn']['enable'] == 0:
            print('begin blast ncbi_bacteria_db')
            blastn_conf = self.conf['denovo']['ncbi_bacteria_blastn']
            blastn_cline = NcbiblastnCommandline(query=assembly_result, db="ncbi_bacteria_db", outfmt=5, out=self.result_dir+"/denovo/ncbi_bacteria_blast_out.xml", num_threads=4, evalue=1e-5, num_alignments=10)
            if not blastn_conf['num_threads'] == 0:
                blastn_cline.set_parameter('num_threads', blastn_conf['num_threads'])
            if not blastn_conf['evalue'] == 0:
                blastn_cline.set_parameter('evalue', blastn_conf['evalue'])
            if not blastn_conf['num_alignments'] == 0:
                blastn_cline.set_parameter('num_alignments', blastn_conf['num_alignments'])
            if not blastn_conf['task'] == 0:
                blastn_cline.set_parameter('task', blastn_conf['task'])
            if not blastn_conf['penalty'] == 0:
                blastn_cline.set_parameter('penalty', blastn_conf['penalty'])
            if not blastn_conf['reward'] == 0:
                blastn_cline.set_parameter('reward', blastn_conf['reward'])
            if not blastn_conf['dust'] == 0:
                blastn_cline.set_parameter('dust', blastn_conf['dust'])
            if not blastn_conf['filtering_db'] == 0:
                blastn_cline.set_parameter('filtering_db', blastn_conf['filtering_db'])
            if not blastn_conf['window_masker_taxid'] == 0:
                blastn_cline.set_parameter('window_masker_taxid', blastn_conf['window_masker_taxid'])
            if not blastn_conf['no_greedy'] == 0:
                blastn_cline.set_parameter('no_greedy', blastn_conf['no_greedy'])
            if not blastn_conf['min_raw_gapped_score'] == 0:
                blastn_cline.set_parameter('min_raw_gapped_score', blastn_conf['min_raw_gapped_score'])
            if not blastn_conf['ungapped'] == 0:
                blastn_cline.set_parameter('ungapped', blastn_conf['ungapped'])
            if not blastn_conf['off_diagonal_range'] == 0:
                blastn_cline.set_parameter('off_diagonal_range', blastn_conf['off_diagonal_range'])
            print(blastn_cline)
            stdout, stderr = blastn_cline()
            print(stdout)
            print(stderr)
            print('end blast ncbi_bacteria_db')

            print('begin handle_denovo_blast_result')
            scripts.handle_denovo_blast_result.handle(self.result_dir, self.result_dir+"/denovo/ncbi_bacteria_blast_out.xml", "ncbi_bacteria_annotation")
            print('end handle_denovo_blast_result')

        # blast patric_amr_db
        if not self.conf['denovo']['patric_blastn']['enable'] == 0:
            print('begin blast patric_amr_db')
            patric_blastn_conf = self.conf['denovo']['patric_blastn']
            blastn_cline = NcbiblastnCommandline(query=assembly_result, db="patric_amr_db", outfmt=5, out=self.result_dir+"/denovo/patric_amr_blast_out.xml", num_threads=4, evalue=1e-5, num_alignments=10)
            if not patric_blastn_conf['num_threads'] == 0:
                blastn_cline.set_parameter('num_threads', patric_blastn_conf['num_threads'])
            if not patric_blastn_conf['evalue'] == 0:
                blastn_cline.set_parameter('evalue', patric_blastn_conf['evalue'])
            if not patric_blastn_conf['num_alignments'] == 0:
                blastn_cline.set_parameter('num_alignments', patric_blastn_conf['num_alignments'])
            if not patric_blastn_conf['task'] == 0:
                blastn_cline.set_parameter('task', patric_blastn_conf['task'])
            if not patric_blastn_conf['penalty'] == 0:
                blastn_cline.set_parameter('penalty', patric_blastn_conf['penalty'])
            if not patric_blastn_conf['reward'] == 0:
                blastn_cline.set_parameter('reward', patric_blastn_conf['reward'])
            if not patric_blastn_conf['dust'] == 0:
                blastn_cline.set_parameter('dust', patric_blastn_conf['dust'])
            if not patric_blastn_conf['filtering_db'] == 0:
                blastn_cline.set_parameter('filtering_db', patric_blastn_conf['filtering_db'])
            if not patric_blastn_conf['window_masker_taxid'] == 0:
                blastn_cline.set_parameter('window_masker_taxid', patric_blastn_conf['window_masker_taxid'])
            if not patric_blastn_conf['no_greedy'] == 0:
                blastn_cline.set_parameter('no_greedy', patric_blastn_conf['no_greedy'])
            if not patric_blastn_conf['min_raw_gapped_score'] == 0:
                blastn_cline.set_parameter('min_raw_gapped_score', patric_blastn_conf['min_raw_gapped_score'])
            if not patric_blastn_conf['ungapped'] == 0:
                blastn_cline.set_parameter('ungapped', patric_blastn_conf['ungapped'])
            if not patric_blastn_conf['off_diagonal_range'] == 0:
                blastn_cline.set_parameter('off_diagonal_range', patric_blastn_conf['off_diagonal_range'])
            print(blastn_cline)
            stdout, stderr = blastn_cline()
            print(stdout)
            print(stderr)
            print('end blast patric_amr_db')

            print('begin handle_denovo_blast_result')
            scripts.handle_denovo_blast_result.handle(self.result_dir, self.result_dir+"/denovo/patric_amr_blast_out.xml", "patric_annotation")
            print('end handle_denovo_blast_result')

        if not self.conf['denovo']['card_blastn']['enable'] == 0:
            # blast card_prevalence_db
            print('begin card_blastn')
            card_blastn_conf = self.conf['denovo']['card_blastn']
            blastn_cline = NcbiblastnCommandline(query=assembly_result, db="card_prevalence_db", outfmt=5, out=self.result_dir+"/denovo/card_prevalence_blast_out.xml", num_threads=4, evalue=1e-5, num_alignments=10)
            if not card_blastn_conf['num_threads'] == 0:
                blastn_cline.set_parameter('num_threads', card_blastn_conf['num_threads'])
            if not card_blastn_conf['evalue'] == 0:
                blastn_cline.set_parameter('evalue', card_blastn_conf['evalue'])
            if not card_blastn_conf['num_alignments'] == 0:
                blastn_cline.set_parameter('num_alignments', card_blastn_conf['num_alignments'])
            if not card_blastn_conf['task'] == 0:
                blastn_cline.set_parameter('task', card_blastn_conf['task'])
            if not card_blastn_conf['penalty'] == 0:
                blastn_cline.set_parameter('penalty', card_blastn_conf['penalty'])
            if not card_blastn_conf['reward'] == 0:
                blastn_cline.set_parameter('reward', card_blastn_conf['reward'])
            if not card_blastn_conf['dust'] == 0:
                blastn_cline.set_parameter('dust', card_blastn_conf['dust'])
            if not card_blastn_conf['filtering_db'] == 0:
                blastn_cline.set_parameter('filtering_db', card_blastn_conf['filtering_db'])
            if not card_blastn_conf['window_masker_taxid'] == 0:
                blastn_cline.set_parameter('window_masker_taxid', card_blastn_conf['window_masker_taxid'])
            if not card_blastn_conf['no_greedy'] == 0:
                blastn_cline.set_parameter('no_greedy', card_blastn_conf['no_greedy'])
            if not card_blastn_conf['min_raw_gapped_score'] == 0:
                blastn_cline.set_parameter('min_raw_gapped_score', card_blastn_conf['min_raw_gapped_score'])
            if not card_blastn_conf['ungapped'] == 0:
                blastn_cline.set_parameter('ungapped', card_blastn_conf['ungapped'])
            if not card_blastn_conf['off_diagonal_range'] == 0:
                blastn_cline.set_parameter('off_diagonal_range', card_blastn_conf['off_diagonal_range'])
            print(blastn_cline)
            stdout, stderr = blastn_cline()
            print(stdout)
            print(stderr)

            scripts.handle_denovo_blast_result.handle(self.result_dir, self.result_dir+"/denovo/card_prevalence_blast_out.xml", "card_annotation")

            print('end card_blastn')

        if not self.conf['denovo']['card_blastn']['enable'] == 0:
            # blast drugbank_db
            print('begin drugbank_db')
            drugbank_blastn_conf = self.conf['denovo']['card_blastn']
            blastn_cline = NcbiblastnCommandline(query=assembly_result, db="drugbank_db", outfmt=5, out=self.result_dir+"/denovo/drugbank_blast_out.xml", num_threads=4, evalue=1e-5, num_alignments=10)
            if not drugbank_blastn_conf['num_threads'] == 0:
                blastn_cline.set_parameter('num_threads', drugbank_blastn_conf['num_threads'])
            if not drugbank_blastn_conf['evalue'] == 0:
                blastn_cline.set_parameter('evalue', drugbank_blastn_conf['evalue'])
            if not drugbank_blastn_conf['num_alignments'] == 0:
                blastn_cline.set_parameter('num_alignments', drugbank_blastn_conf['num_alignments'])
            if not drugbank_blastn_conf['task'] == 0:
                blastn_cline.set_parameter('task', drugbank_blastn_conf['task'])
            if not drugbank_blastn_conf['penalty'] == 0:
                blastn_cline.set_parameter('penalty', drugbank_blastn_conf['penalty'])
            if not drugbank_blastn_conf['reward'] == 0:
                blastn_cline.set_parameter('reward', drugbank_blastn_conf['reward'])
            if not drugbank_blastn_conf['dust'] == 0:
                blastn_cline.set_parameter('dust', drugbank_blastn_conf['dust'])
            if not drugbank_blastn_conf['filtering_db'] == 0:
                blastn_cline.set_parameter('filtering_db', drugbank_blastn_conf['filtering_db'])
            if not drugbank_blastn_conf['window_masker_taxid'] == 0:
                blastn_cline.set_parameter('window_masker_taxid', drugbank_blastn_conf['window_masker_taxid'])
            if not drugbank_blastn_conf['no_greedy'] == 0:
                blastn_cline.set_parameter('no_greedy', drugbank_blastn_conf['no_greedy'])
            if not drugbank_blastn_conf['min_raw_gapped_score'] == 0:
                blastn_cline.set_parameter('min_raw_gapped_score', drugbank_blastn_conf['min_raw_gapped_score'])
            if not drugbank_blastn_conf['ungapped'] == 0:
                blastn_cline.set_parameter('ungapped', drugbank_blastn_conf['ungapped'])
            if not drugbank_blastn_conf['off_diagonal_range'] == 0:
                blastn_cline.set_parameter('off_diagonal_range', patric_blastn_conf['off_diagonal_range'])
            print(blastn_cline)
            stdout, stderr = blastn_cline()
            print(stdout)
            print(stderr)

            scripts.handle_denovo_blast_result.handle(self.result_dir, self.result_dir+"/denovo/drugbank_blast_out.xml", "drugbank_annotation")
        
            print('end drugbank_db')

        if not self.conf['denovo']['prokka']['enable'] == 0:
            scripts.prokka_annotation.annotation(assembly_result, self.result_dir, self.conf)
        if not self.conf['denovo']['prodigal']['enable'] == 0:
            scripts.prodigal_prediction.prediction(assembly_result, self.result_dir)