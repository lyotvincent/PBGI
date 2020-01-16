from Bio.Blast.Applications import NcbiblastnCommandline
import subprocess
import os
import time

import assembly
import scripts.fasta2fastq
import scripts.fastq2fasta
import scripts.handle_blast_xml_result
import scripts.download_nucleotide
import scripts.download_gb
import scripts.snap
import scripts.genbank2gff3
import scripts.annotation_generator

class Resequencing:

    result_dir = None
    conf = None
    input_file = None
    input_file_paired = None

    def __init__(self, result_dir, conf, input_file, input_file_paired=None):
        self.result_dir = result_dir
        self.conf = conf
        self.input_file = input_file
        self.input_file_paired = input_file_paired
    
    def run(self):

        print('Begin Resequencing')

        print('input_file='+str(self.input_file))
        
        assembly_conf = self.conf['resequencing']['assembly']
        if not assembly_conf['enable'] == 0:
            if self.input_file_paired == None:
                print('begin assembly')
                os.mkdir(os.path.abspath('.')+'/'+self.result_dir+'/resequencing/assembly')
                assembly_obj = assembly.Assembly(self.result_dir+'/resequencing/assembly', self.conf, input_file=self.input_file)
                print(assembly_conf['megahit']['enable'], assembly_conf['spades']['enable'], assembly_conf['velvet']['enable'])
                if not assembly_conf['megahit']['enable'] == 0:
                    assembly_obj.megahit_single()
                    assembly_result = self.result_dir+'/resequencing/assembly/megahit_out/final.contigs.fa'
                elif not assembly_conf['spades']['enable'] == 0:
                    assembly_obj.spades_single()
                    assembly_result = self.result_dir+'/resequencing/assembly/spades_out/contigs.fasta'
                elif not assembly_conf['velvet']['enable'] == 0:
                    assembly_obj.velvet_single()
                    assembly_result = self.result_dir+'/resequencing/assembly/velvet_output/contigs.fa'
                else:
                    assembly_obj.megahit_single()
                    assembly_result = self.result_dir+'/resequencing/assembly/megahit_out/final.contigs.fa'
                print('end assembly')
                print("Begin QUAST")
                subprocess.run('quast.py '+assembly_result+' --min-contig 50 -o '+self.result_dir+'/resequencing/quast_out', shell=True, check=True)
                print("End QUAST")
            else:
                print('begin assembly')
                os.mkdir(os.path.abspath('.')+'/'+self.result_dir+'/resequencing/assembly')
                assembly_obj = assembly.Assembly(self.result_dir+'/resequencing/assembly', self.conf, input_file_1=self.input_file, input_file_2=self.input_file_paired)
                if not assembly_conf['megahit']['enable'] == 0:
                    assembly_obj.megahit_paired()
                    assembly_result = self.result_dir+'/resequencing/assembly/megahit_out/final.contigs.fa'
                elif not assembly_conf['spades']['enable'] == 0:
                    assembly_obj.spades_paired()
                    assembly_result = self.result_dir+'/resequencing/assembly/spades_out/contigs.fasta'
                elif not assembly_conf['velvet']['enable'] == 0:
                    assembly_obj.velvet_paired()
                    assembly_result = self.result_dir+'/resequencing/assembly/velvet_output/contigs.fa'
                else:
                    assembly_obj.megahit_paired()
                    assembly_result = self.result_dir+'/resequencing/assembly/megahit_out/final.contigs.fa'
                print('end assembly')
                print("Begin QUAST")
                subprocess.run('quast.py '+assembly_result+' --min-contig 50 -o '+self.result_dir+'/resequencing/quast_out', shell=True, check=True)
                print("End QUAST")
            blast_input = assembly_result
        else:
            scripts.fastq2fasta.fastq2fasta(self.input_file, self.result_dir+"/resequencing")
            blast_input = self.result_dir+"/resequencing/input.fasta"

        # time.sleep(1000)
        
        blastn_conf = self.conf['resequencing']['blastn']
        blastn_cline = NcbiblastnCommandline(query=blast_input, db="ncbi_bacteria_db", outfmt=5, out=self.result_dir+"/resequencing/ncbi_bacteria_blast_out.xml", num_threads=4, evalue=1e-5, num_alignments=10)
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

        accession_version_dict = scripts.handle_blast_xml_result.handle_blast_xml_result(self.result_dir+"/resequencing/ncbi_bacteria_blast_out.xml")
        # print("accession_version_dict="+str(accession_version_dict))

        accession_version_list = []
        i = 0
        for accession_version in accession_version_dict:
            if i < self.conf['resequencing']['number_of_candidate_similar_genome']:
                i += 1
            else:
                break
            accession_version_list.append(accession_version[0])
            scripts.download_nucleotide.download_by_accession_version(self.result_dir, accession_version[0])
        
        closest_accession_version = scripts.snap.snap(self.result_dir, accession_version_list, self.input_file, self.input_file_paired)
        print("closest_accession_version="+str(closest_accession_version))

        alignment_tool = 'bowtie2' if not self.conf['resequencing']['alignment_tool']['bowtie2'] == 0 else "snap"
        if alignment_tool == 'bowtie2':
            subprocess.run('bowtie2-build '+self.result_dir+'/resequencing/'+closest_accession_version+'.fasta '+self.result_dir+'/resequencing/'+closest_accession_version+'_bowtie2_index', shell=True, check=True)
        else:
            subprocess.run('snap-aligner index '+self.result_dir+'/resequencing/'+closest_accession_version+'.fasta '+self.result_dir+'/resequencing/'+closest_accession_version+'_index', shell=True, check=True)
        if self.input_file_paired == None:
            assembly_conf = self.conf['resequencing']['assembly']
            if not assembly_conf['enable'] == 0:
                if alignment_tool == 'bowtie2':
                    subprocess.run('bowtie2 -x '+self.result_dir+'/resequencing/'+closest_accession_version+'_bowtie2_index -f -U '+assembly_result+' -S '+self.result_dir+'/resequencing/alignment_result.sam', shell=True, check=True)
                else:
                    scripts.fasta2fastq.fasta_to_fastq(assembly_result, self.result_dir)
                    assembly_result = self.result_dir+"/resequencing/assembly_result.fastq"
                    subprocess.run('snap-aligner single '+self.result_dir+'/resequencing/'+closest_accession_version+'_index '+assembly_result+' -o '+self.result_dir+'/resequencing/alignment_result.sam', shell=True, check=True)
            else:
                if alignment_tool == 'bowtie2':
                    subprocess.run('bowtie2 -x '+self.result_dir+'/resequencing/'+closest_accession_version+'_bowtie2_index -U '+self.input_file+' -S '+self.result_dir+'/resequencing/alignment_result.sam', shell=True, check=True)
                else:
                    subprocess.run('snap-aligner single '+self.result_dir+'/resequencing/'+closest_accession_version+'_index '+self.input_file+' -o '+self.result_dir+'/resequencing/alignment_result.sam', shell=True, check=True)
            # subprocess.run('bowtie2 -x '+self.result_dir+'/resequencing/'+closest_accession_version+'_bowtie2_index -U '+self.input_file+' -S '+self.result_dir+'/resequencing/alignment_result.sam', shell=True, check=True)
        else:
            assembly_conf = self.conf['resequencing']['assembly']
            if not assembly_conf['enable'] == 0:
                if alignment_tool == 'bowtie2':
                    subprocess.run('bowtie2 -x '+self.result_dir+'/resequencing/'+closest_accession_version+'_bowtie2_index -f -U '+assembly_result+' -S '+self.result_dir+'/resequencing/alignment_result.sam', shell=True, check=True)
                else:
                    scripts.fasta2fastq.fasta_to_fastq(assembly_result, self.result_dir)
                    assembly_result = self.result_dir+"/resequencing/assembly_result.fastq"
                    subprocess.run('snap-aligner single '+self.result_dir+'/resequencing/'+closest_accession_version+'_index '+assembly_result+' -o '+self.result_dir+'/resequencing/alignment_result.sam', shell=True, check=True)
            else:
                if alignment_tool == 'bowtie2':
                    subprocess.run('bowtie2 -x '+self.result_dir+'/resequencing/'+closest_accession_version+'_bowtie2_index -1 '+self.input_file+' -2 '+self.input_file_paired+' -S '+self.result_dir+'/resequencing/alignment_result.sam', shell=True, check=True)
                else:
                    subprocess.run('snap-aligner paired '+self.result_dir+'/resequencing/'+closest_accession_version+'_index '+self.input_file+' '+self.input_file_paired+' -o '+self.result_dir+'/resequencing/alignment_result.sam', shell=True, check=True)
            # subprocess.run('bowtie2 -x '+self.result_dir+'/resequencing/'+closest_accession_version+'_bowtie2_index -1 '+self.input_file+' -2 '+self.input_file_paired+' -S '+self.result_dir+'/resequencing/alignment_result.sam', shell=True, check=True)

        print("begin sort_sam")
        subprocess.run('samtools view -bS '+self.result_dir+'/resequencing/alignment_result.sam > '+self.result_dir+'/resequencing/output.bam', shell=True, check=True)
        subprocess.run('samtools sort '+self.result_dir+'/resequencing/output.bam -o '+self.result_dir+'/resequencing/alignment_result.sorted.sam -O sam', shell=True, check=True)
        print("end sort_sam")

        scripts.download_gb.download_by_accession_version(self.result_dir, closest_accession_version)
        # subprocess.run('perl scripts/download_genbank.pl '+closest_accession_version+' '+self.result_dir+'/resequencing', shell=True, check=True)
        transform = scripts.genbank2gff3.Genbank2gff3(closest_accession_version+'.gb', self.result_dir)
        # transform.genbank2gff3_by_bioperl()
        transform.genbank2gff3_by_bcbio()

        annotaion_generator = scripts.annotation_generator.Annotation_generator(self.result_dir+'/resequencing/alignment_result.sorted.sam', self.result_dir+'/resequencing/'+closest_accession_version+'.gff', self.result_dir)
        annotaion_generator.look_for_annotation()
        
    def run_3gs(self):

        print('Begin Resequencing 3gs')

        scripts.fastq2fasta.fastq2fasta(self.input_file, self.result_dir+"/resequencing")
        
        time.sleep(900)

        blastn_conf = self.conf['resequencing']['blastn']
        blastn_cline = NcbiblastnCommandline(query=self.result_dir+"/resequencing/input.fasta", db="ncbi_bacteria_db", outfmt=5, out=self.result_dir+"/resequencing/ncbi_bacteria_blast_out.xml", num_threads=4, evalue=1e-5, num_alignments=10)
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
        # stdout, stderr = blastn_cline()
        # print(stdout)
        # print(stderr)

        accession_version_dict = scripts.handle_blast_xml_result.handle_blast_xml_result(self.result_dir+"/resequencing/ncbi_bacteria_blast_out.xml")
        # print(accession_version_dict)

        # accession_version_list = []
        # # i = 0
        # for accession_version in accession_version_dict:
        #     # if i < 5:
        #     #     i += 1
        #     # else:
        #     #     break
        #     accession_version_list.append(accession_version[0])
        #     break
        #     # scripts.download_nucleotide.download_by_accession_version(self.result_dir, accession_version[0])
        
        # # closest_accession_version = scripts.snap.snap(self.result_dir, accession_version_list, self.input_file, self.input_file_paired)
        # scripts.download_nucleotide.download_by_accession_version(self.result_dir, accession_version[0])
        # closest_accession_version = accession_version_list[0]

        accession_version_list = []
        i = 0
        for accession_version in accession_version_dict:
            if i < self.conf['resequencing']['number_of_candidate_similar_genome']:
                i += 1
            else:
                break
            accession_version_list.append(accession_version[0])
            scripts.download_nucleotide.download_by_accession_version(self.result_dir, accession_version[0])
        
        # minimap2
        minimap2_conf = self.conf['resequencing']['minimap2']
        command_line = "minimap2 -a "
        if not minimap2_conf['-H'] == 0:
            command_line += "-H "
        if not minimap2_conf['-k'] == 0:
            command_line += "-k "+str(minimap2_conf['-k'])+' '
        if not minimap2_conf['-w'] == 0:
            command_line += "-w "+str(minimap2_conf['-w'])+' '
        if not minimap2_conf['-I'] == 0:
            command_line += "-I "+str(minimap2_conf['-I'])+' '
        if not minimap2_conf['-f'] == 0:
            command_line += "-f "+str(minimap2_conf['-f'])+' '
        if not minimap2_conf['-g'] == 0:
            command_line += "-g "+str(minimap2_conf['-g'])+' '
        if not minimap2_conf['-G'] == 0:
            command_line += "-G "+str(minimap2_conf['-G'])+' '
        if not minimap2_conf['-F'] == 0:
            command_line += "-F "+str(minimap2_conf['-F'])+' '
        if not minimap2_conf['-r'] == 0:
            command_line += "-r "+str(minimap2_conf['-r'])+' '
        if not minimap2_conf['-n'] == 0:
            command_line += "-n "+str(minimap2_conf['-n'])+' '
        if not minimap2_conf['-m'] == 0:
            command_line += "-m "+str(minimap2_conf['-m'])+' '
        if not minimap2_conf['-X'] == 0:
            command_line += "-X "
        if not minimap2_conf['-p'] == 0:
            command_line += "-p "+str(minimap2_conf['-p'])+' '
        if not minimap2_conf['-N'] == 0:
            command_line += "-N "+str(minimap2_conf['-N'])+' '
        if not minimap2_conf['-A'] == 0:
            command_line += "-A "+str(minimap2_conf['-A'])+' '
        if not minimap2_conf['-B'] == 0:
            command_line += "-B "+str(minimap2_conf['-B'])+' '
        if not minimap2_conf['-O'] == 0:
            command_line += "-O "+str(minimap2_conf['-O'])+' '
        if not minimap2_conf['-E'] == 0:
            command_line += "-E "+str(minimap2_conf['-E'])+' '
        if not minimap2_conf['-z'] == 0:
            command_line += "-z "+str(minimap2_conf['-z'])+' '
        if not minimap2_conf['-s'] == 0:
            command_line += "-s "+str(minimap2_conf['-s'])+' '
        if not minimap2_conf['-u'] == 0:
            command_line += "-u "+str(minimap2_conf['-u'])+' '
        if not minimap2_conf['-x'] == 0:
            command_line += "-x "+str(minimap2_conf['-x'])+' '
        command_line += self.result_dir+"/resequencing/%s.fasta "+self.input_file+" > "+self.result_dir+"/resequencing/minimap2_%s_result.sam"


        closest_accession_version = accession_version_list[0]
        max_similarity = 0
        for accession_version in accession_version_list:
            print(command_line % (accession_version, accession_version))
            subprocess.run(command_line % (accession_version, accession_version), shell=True, check=True)
            subprocess.run('samtools sort '+self.result_dir+'/resequencing/minimap2_'+accession_version+'_result.sam -o '+self.result_dir+'/resequencing/minimap2_'+accession_version+'_result.sorted.bam -O bam', shell=True, check=True)
            res = subprocess.Popen('bedtools genomecov -ibam '+self.result_dir+'/resequencing/minimap2_'+accession_version+'_result.sorted.bam', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
            send_messages = res.stdout.readlines()
            print(send_messages[0].decode('utf-8'))
            for message in send_messages:
                message = message.decode('utf-8')
                if message.startswith('genome\t0'):
                    print(message)
                    line = message.split()
                    if line[-1] == '\n':
                        new_similarity = 1 - float(line[-2])
                    else:
                        new_similarity = 1 - float(line[-1])
                    break
            print('new_similarity='+str(new_similarity))
            if new_similarity > max_similarity:
                max_similarity = new_similarity
                closest_accession_version = accession_version
        print("closest_accession_version="+str(closest_accession_version))


        print("begin sort_sam")
        # subprocess.run('samtools view -bS '+self.result_dir+'/resequencing/minimap2_result.sam > '+self.result_dir+'/resequencing/output.bam', shell=True, check=True)
        # subprocess.run('samtools sort '+self.result_dir+'/resequencing/output.bam -o '+self.result_dir+'/resequencing/minimap2_result.sorted.sam -O sam', shell=True, check=True)
        subprocess.run('samtools sort '+self.result_dir+'/resequencing/minimap2_'+closest_accession_version+'_result.sam -o '+self.result_dir+'/resequencing/minimap2_'+closest_accession_version+'_result.sorted.sam -O sam', shell=True, check=True)
        print("end sort_sam")
        
        scripts.download_gb.download_by_accession_version(self.result_dir, closest_accession_version)
        # subprocess.run('perl scripts/download_genbank.pl '+closest_accession_version+' '+self.result_dir+'/resequencing', shell=True, check=True)
        transform = scripts.genbank2gff3.Genbank2gff3(closest_accession_version+'.gb', self.result_dir)
        # transform.genbank2gff3_by_bioperl()
        transform.genbank2gff3_by_bcbio()

        annotaion_generator = scripts.annotation_generator.Annotation_generator(self.result_dir+'/resequencing/minimap2_'+closest_accession_version+'_result.sorted.sam', self.result_dir+'/resequencing/'+closest_accession_version+'.gff', self.result_dir)
        annotaion_generator.look_for_annotation()