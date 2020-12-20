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
        assembly_obj = assembly.Assembly(self.result_dir+'/denovo/assembly', assembly_conf, input_file=self.input_file)
        if assembly_conf['megahit']['enable']:
            assembly_obj.megahit_single()
            assembly_result = self.result_dir+'/denovo/assembly/megahit_out/final.contigs.fa'
        elif assembly_conf['spades']['enable']:
            assembly_obj.spades_single()
            assembly_result = self.result_dir+'/denovo/assembly/spades_out/contigs.fasta'
        elif assembly_conf['velvet']['enable']:
            assembly_obj.velvet_single()
            assembly_result = self.result_dir+'/denovo/assembly/velvet_output/contigs.fa'
        else:
            assembly_obj.megahit_single()
            assembly_result = self.result_dir+'/denovo/assembly/megahit_out/final.contigs.fa'
        print('end assembly')

        temp_file = open(os.path.abspath('.') + '/' + self.result_dir+'/Summary_of_results.html', 'a+')
        temp_file.write('<ul>\n')
        temp_file.write('<li>assembly result is in %s</li>\n' % './denovo/assembly')
        temp_file.write('<li><a href="%s">click to assembly result</a></li>\n' % ('./denovo/assembly'))
        temp_file.write('<li>assembly qc result is in %s</li>\n' % './denovo/quast_out')
        temp_file.write('<li><a href="%s">click to assembly qc result</a></li>\n' % ('./denovo/quast_out/report.html'))
        temp_file.write('</ul>\n')
        temp_file.close()

        print("Begin QUAST")
        subprocess.run('quast '+assembly_result+' --min-contig 50 -o '+self.result_dir+'/denovo/quast_out', shell=True, check=True)
        print("End QUAST")

        self.blastn_and_annotation(assembly_result)

        print('end single')

    def run_paired(self):
        print('denovo paired')
        print('begin assembly')
        assembly_conf = self.conf['denovo']['assembly']
        os.mkdir(os.path.abspath('.')+'/'+self.result_dir+'/denovo/assembly')
        assembly_obj = assembly.Assembly(self.result_dir+'/denovo/assembly', assembly_conf, input_file_1=self.input_file_1, input_file_2=self.input_file_2)
        if assembly_conf['megahit']['enable']:
            assembly_obj.megahit_paired()
            assembly_result = self.result_dir+'/denovo/assembly/megahit_out/final.contigs.fa'
        elif assembly_conf['spades']['enable']:
            assembly_obj.spades_paired()
            assembly_result = self.result_dir+'/denovo/assembly/spades_out/contigs.fasta'
        elif assembly_conf['velvet']['enable']:
            assembly_obj.velvet_paired()
            assembly_result = self.result_dir+'/denovo/assembly/velvet_output/contigs.fa'
        else:
            assembly_obj.megahit_paired()
            assembly_result = self.result_dir+'/denovo/assembly/megahit_out/final.contigs.fa'
        print('end assembly')

        print("Begin QUAST")
        subprocess.run('quast '+assembly_result+' --min-contig 50 -o '+self.result_dir+'/denovo/quast_out', shell=True, check=True)
        print("End QUAST")

        self.blastn_and_annotation(assembly_result)

        print('end paired')

    def run_interleaved(self):
        print('denovo interleaved')
        assembly_conf = self.conf['denovo']['assembly']
        os.mkdir(os.path.abspath('.')+'/'+self.result_dir+'/denovo/assembly')
        assembly_obj = assembly.Assembly(self.result_dir+'/denovo/assembly', assembly_conf, input_file_12=self.input_file_12)
        if assembly_conf['megahit']['enable']:
            assembly_obj.megahit_interleaved()
            assembly_result = self.result_dir+'/denovo/assembly/megahit_out/final.contigs.fa'
        elif assembly_conf['spades']['enable']:
            assembly_obj.spades_interlaced()
            assembly_result = self.result_dir+'/denovo/assembly/spades_out/contigs.fasta'
        elif assembly_conf['velvet']['enable']:
            assembly_obj.velvet_interleaved()
            assembly_result = self.result_dir+'/denovo/assembly/velvet_output/contigs.fa'
        else:
            assembly_obj.megahit_interleaved()
            assembly_result = self.result_dir+'/denovo/assembly/megahit_out/final.contigs.fa'

        print("Begin QUAST")
        subprocess.run('quast '+assembly_result+' --min-contig 50 -o '+self.result_dir+'/denovo/quast_out', shell=True, check=True)
        print("End QUAST")

        self.blastn_and_annotation(assembly_result)

        print('end interleaved')

    def run_3gs(self):
        print('begin denovo 3gs')

        assembly_conf = self.conf['denovo']['assembly']
        if assembly_conf['enable']:
            os.mkdir(os.path.abspath('.')+'/'+self.result_dir+'/denovo/assembly')
            assembly_obj = assembly.Assembly(self.result_dir+'/denovo/assembly', assembly_conf, input_file=self.input_file)
            if assembly_conf['canu']['enable']:
                assembly_obj.canu()
                assembly_result = self.result_dir+'/denovo/assembly/canu_output/canu_assembly_result.contigs.fasta'
            elif assembly_conf['spades']['enable']:
                assembly_obj.spades_3gs()
                assembly_result = self.result_dir+'/denovo/assembly/spades_out/contigs.fasta'
            else:
                assembly_obj.canu()
                assembly_result = self.result_dir+'/denovo/assembly/canu_output/canu_assembly_result.contigs.fasta'

            print("Begin QUAST")
            subprocess.run('quast '+assembly_result+' --min-contig 50 -o '+self.result_dir+'/denovo/quast_out', shell=True, check=True)
            print("End QUAST")
        else:
            scripts.fastq2fasta.fastq2fasta(self.input_file, self.result_dir+"/denovo")
            assembly_result = self.result_dir+"/denovo/input.fasta"

        self.blastn_and_annotation(assembly_result)

        print('end denovo 3gs')

    def blastn_and_annotation(self, assembly_result):
        if self.conf['denovo']['ncbi_bacteria_blastn']['enable']:
            print('begin blast ncbi_bacteria_db')
            blastn_conf = self.conf['denovo']['ncbi_bacteria_blastn']
            blastn_cline = NcbiblastnCommandline(cmd='blastn', query=assembly_result, db=self.conf["denovo"]["ncbi_bacteria_blastn"]["blast_db_path"], outfmt=5, out=self.result_dir+"/denovo/ncbi_bacteria_blast_out.xml")
            blastn_cline.set_parameter('num_threads', int(blastn_conf['num_threads']))
            blastn_cline.set_parameter('num_alignments', int(blastn_conf['num_alignments']))
            if blastn_conf['evalue'] != None:
                blastn_cline.set_parameter('evalue', float(blastn_conf['evalue']))
            if blastn_conf['penalty'] != None:
                blastn_cline.set_parameter('penalty', int(blastn_conf['penalty']))
            if blastn_conf['reward'] != None:
                blastn_cline.set_parameter('reward', int(blastn_conf['reward']))
            print(blastn_cline)
            stdout, stderr = blastn_cline()
            print(stdout)
            print(stderr)
            print('end blast ncbi_bacteria_db')

            print('begin handle_denovo_blast_result')
            scripts.handle_denovo_blast_result.handle(self.result_dir, self.result_dir+"/denovo/ncbi_bacteria_blast_out.xml", "ncbi_bacteria_annotation")
            print('end handle_denovo_blast_result')

            temp_file = open(os.path.abspath('.') + '/' + self.result_dir+'/Summary_of_results.html', 'a+')
            temp_file.write('<ul>\n')
            temp_file.write('<li>ncbi bacteria annotation is in %s</li>\n' % ('./denovo/ncbi_bacteria_annotation'))
            temp_file.write('<li><a href="%s">click to ncbi bacteria annotation</a></li>\n' % ('./denovo/ncbi_bacteria_annotation'))
            temp_file.write('</ul>\n')
            temp_file.close()

        # blast patric_amr_db
        if self.conf['denovo']['patric_blastn']['enable']:
            print('begin blast patric_amr_db')
            patric_blastn_conf = self.conf['denovo']['patric_blastn']
            blastn_cline = NcbiblastnCommandline(cmd='blastn', query=assembly_result, db=self.conf["denovo"]["patric_blastn"]["patric_db_path"], outfmt=5, out=self.result_dir+"/denovo/patric_amr_blast_out.xml")
            blastn_cline.set_parameter('num_threads', int(patric_blastn_conf['num_threads']))
            blastn_cline.set_parameter('num_alignments', int(patric_blastn_conf['num_alignments']))
            if patric_blastn_conf['evalue'] != None:
                blastn_cline.set_parameter('evalue', float(patric_blastn_conf['evalue']))
            if patric_blastn_conf['penalty'] != None:
                blastn_cline.set_parameter('penalty', int(patric_blastn_conf['penalty']))
            if patric_blastn_conf['reward'] != None:
                blastn_cline.set_parameter('reward', int(patric_blastn_conf['reward']))
            print(blastn_cline)
            stdout, stderr = blastn_cline()
            print(stdout)
            print(stderr)

            print('begin handle_denovo_blast_result')
            scripts.handle_denovo_blast_result.handle(self.result_dir, self.result_dir+"/denovo/patric_amr_blast_out.xml", "patric_annotation")
            print('end handle_denovo_blast_result')
            print('end blast patric_amr_db')

            temp_file = open(os.path.abspath('.') + '/' + self.result_dir+'/Summary_of_results.html', 'a+')
            temp_file.write('<ul>\n')
            temp_file.write('<li>patric annotation is in %s</li>\n' % ('./denovo/patric_annotation'))
            temp_file.write('<li><a href="%s">click to patric annotation</a></li>\n' % ('./denovo/patric_annotation'))
            temp_file.write('</ul>\n')
            temp_file.close()

        if self.conf['denovo']['card_blastn']['enable']:
            # blast card_prevalence_db
            print('begin card_blastn')
            card_blastn_conf = self.conf['denovo']['card_blastn']
            blastn_cline = NcbiblastnCommandline(cmd='blastn', query=assembly_result, db=self.conf["denovo"]["card_blastn"]["card_db_path"], outfmt=5, out=self.result_dir+"/denovo/card_prevalence_blast_out.xml")
            blastn_cline.set_parameter('num_threads', int(card_blastn_conf['num_threads']))
            blastn_cline.set_parameter('num_alignments', int(card_blastn_conf['num_alignments']))
            if card_blastn_conf['evalue'] != None:
                blastn_cline.set_parameter('evalue', float(card_blastn_conf['evalue']))
            if card_blastn_conf['penalty'] != None:
                blastn_cline.set_parameter('penalty', int(card_blastn_conf['penalty']))
            if card_blastn_conf['reward'] != None:
                blastn_cline.set_parameter('reward', int(card_blastn_conf['reward']))
            print(blastn_cline)
            stdout, stderr = blastn_cline()
            print(stdout)
            print(stderr)

            scripts.handle_denovo_blast_result.handle(self.result_dir, self.result_dir+"/denovo/card_prevalence_blast_out.xml", "card_annotation")

            print('end card_blastn')

            temp_file = open(os.path.abspath('.') + '/' + self.result_dir+'/Summary_of_results.html', 'a+')
            temp_file.write('<ul>\n')
            temp_file.write('<li>card annotation is in %s</li>\n' % ('./denovo/card_annotation'))
            temp_file.write('<li><a href="%s">click to card annotation</a></li>\n' % ('./denovo/card_annotation'))
            temp_file.write('</ul>\n')
            temp_file.close()

        if self.conf['denovo']['drugbank_blastn']['enable']:
            # blast drugbank_db
            print('begin drugbank_db')
            drugbank_blastn_conf = self.conf['denovo']['card_blastn']
            blastn_cline = NcbiblastnCommandline(cmd='blastn', query=assembly_result, db=self.conf["denovo"]["drugbank_blastn"]["drugbank_db_path"], outfmt=5, out=self.result_dir+"/denovo/drugbank_blast_out.xml")
            blastn_cline.set_parameter('num_threads', int(drugbank_blastn_conf['num_threads']))
            blastn_cline.set_parameter('num_alignments', int(drugbank_blastn_conf['num_alignments']))
            if drugbank_blastn_conf['evalue'] != None:
                blastn_cline.set_parameter('evalue', float(drugbank_blastn_conf['evalue']))
            if drugbank_blastn_conf['penalty'] != None:
                blastn_cline.set_parameter('penalty', int(drugbank_blastn_conf['penalty']))
            if drugbank_blastn_conf['reward'] != None:
                blastn_cline.set_parameter('reward', int(drugbank_blastn_conf['reward']))
            print(blastn_cline)
            stdout, stderr = blastn_cline()
            print(stdout)
            print(stderr)

            scripts.handle_denovo_blast_result.handle(self.result_dir, self.result_dir+"/denovo/drugbank_blast_out.xml", "drugbank_annotation")
        
            print('end drugbank_db')
            temp_file = open(os.path.abspath('.') + '/' + self.result_dir+'/Summary_of_results.html', 'a+')
            temp_file.write('<ul>\n')
            temp_file.write('<li>drugbank annotation is in %s</li>\n' % ('./denovo/drugbank_annotation'))
            temp_file.write('<li><a href="%s">click to drugbank annotation</a></li>\n' % ('./denovo/drugbank_annotation'))
            temp_file.write('</ul>\n')
            temp_file.close()

        if self.conf['denovo']['prokka']['enable']:
            scripts.prokka_annotation.annotation(assembly_result, self.result_dir, self.conf)
        if self.conf['denovo']['prodigal']['enable']:
            scripts.prodigal_prediction.prediction(assembly_result, self.result_dir)