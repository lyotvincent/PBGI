from Bio.Blast.Applications import NcbiblastnCommandline
import subprocess
import os
import time

import assembly
import calculating_genome_similarity
# import scripts.fasta2fastq
import scripts.fastq2fasta
import scripts.handle_blast_xml_result
import scripts.handle_denovo_blast_result
import scripts.download_nucleotide
import scripts.download_gb
# import scripts.snap
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
        # self.data_path = data_path
        self.input_file = input_file
        self.input_file_paired = input_file_paired
    
    def run(self):

        print('Begin Resequencing')

        print('input_file='+str(self.input_file))
        
        assembly_conf = self.conf['resequencing']['assembly']
        if assembly_conf['enable']:
            if self.input_file_paired == None:
                print('begin assembly')
                os.mkdir(os.path.abspath('.')+'/'+self.result_dir+'/resequencing/assembly')
                assembly_obj = assembly.Assembly(self.result_dir+'/resequencing/assembly', assembly_conf, input_file=self.input_file)
                print(assembly_conf['megahit']['enable'], assembly_conf['spades']['enable'], assembly_conf['velvet']['enable'])
                if assembly_conf['megahit']['enable']:
                    assembly_obj.megahit_single()
                    assembly_result = self.result_dir+'/resequencing/assembly/megahit_out/final.contigs.fa'
                elif assembly_conf['spades']['enable']:
                    assembly_obj.spades_single()
                    assembly_result = self.result_dir+'/resequencing/assembly/spades_out/contigs.fasta'
                elif assembly_conf['velvet']['enable']:
                    assembly_obj.velvet_single()
                    assembly_result = self.result_dir+'/resequencing/assembly/velvet_output/contigs.fa'
                else:
                    assembly_obj.megahit_single()
                    assembly_result = self.result_dir+'/resequencing/assembly/megahit_out/final.contigs.fa'
                print('end assembly')
                print("Begin QUAST")
                subprocess.run('quast '+assembly_result+' --min-contig 50 -o '+self.result_dir+'/resequencing/quast_out', shell=True, check=True)
                print("End QUAST")
            else:
                print('begin assembly')
                os.mkdir(os.path.abspath('.')+'/'+self.result_dir+'/resequencing/assembly')
                assembly_obj = assembly.Assembly(self.result_dir+'/resequencing/assembly', assembly_conf, input_file_1=self.input_file, input_file_2=self.input_file_paired)
                if assembly_conf['megahit']['enable']:
                    assembly_obj.megahit_paired()
                    assembly_result = self.result_dir+'/resequencing/assembly/megahit_out/final.contigs.fa'
                elif assembly_conf['spades']['enable']:
                    assembly_obj.spades_paired()
                    assembly_result = self.result_dir+'/resequencing/assembly/spades_out/contigs.fasta'
                elif assembly_conf['velvet']['enable']:
                    assembly_obj.velvet_paired()
                    assembly_result = self.result_dir+'/resequencing/assembly/velvet_output/contigs.fa'
                else:
                    assembly_obj.megahit_paired()
                    assembly_result = self.result_dir+'/resequencing/assembly/megahit_out/final.contigs.fa'
                print('end assembly')
                print("Begin QUAST")
                subprocess.run('quast '+assembly_result+' --min-contig 50 -o '+self.result_dir+'/resequencing/quast_out', shell=True, check=True)
                print("End QUAST")
            blast_input = assembly_result
            
            temp_file = open(os.path.abspath('.') + '/' + self.result_dir+'/Summary_of_results.html', 'a+')
            temp_file.write('<ul>\n')
            temp_file.write('<li>assembly result is in %s</li>\n' % './resequencing/assembly')
            temp_file.write('<li><a href="%s">click to assembly result</a></li>\n' % ('./resequencing/assembly'))
            temp_file.write('<li>assembly qc result is in %s</li>\n' % './resequencing/quast_out')
            temp_file.write('<li><a href="%s">click to assembly qc result</a></li>\n' % ('./resequencing/quast_out/report.html'))
            temp_file.write('</ul>\n')
            temp_file.close()
        else:
            scripts.fastq2fasta.fastq2fasta(self.input_file, self.result_dir+"/resequencing")
            blast_input = self.result_dir+"/resequencing/input.fasta"

        if assembly_conf['enable']:
            blastn_conf = self.conf['resequencing']['blastn']
            
            # blastn_cline = NcbiblastnCommandline(cmd=os.path.dirname(os.path.realpath(__file__))+'/external_tools/blastn', query=blast_input, db=self.conf["resequencing"]['blastn']["blast_db_path"], outfmt=7, out=self.result_dir+"/resequencing/ncbi_bacteria_blast_out.xml")
            # blastn_cline.set_parameter('num_threads', int(blastn_conf['num_threads']))
            # blastn_cline.set_parameter('num_alignments', int(blastn_conf['num_alignments']))
            # if blastn_conf['evalue'] != None:
            #     blastn_cline.set_parameter('evalue', float(blastn_conf['evalue']))
            # # if blastn_conf['task'] != None:
            # #     blastn_cline.set_parameter('task', blastn_conf['task'])
            # if blastn_conf['penalty'] != None:
            #     blastn_cline.set_parameter('penalty', int(blastn_conf['penalty']))
            # if blastn_conf['reward'] != None:
            #     blastn_cline.set_parameter('reward', int(blastn_conf['reward']))
            # print(blastn_cline)
            # stdout, stderr = blastn_cline()
            # print(stdout)
            # print(stderr)

            blastn_line = """%s -query %s -db %s -outfmt 7 -num_threads %s -num_alignments %s -evalue %s""" % ('blastn', assembly_result, blastn_conf["blast_db_path"], blastn_conf['num_threads'], blastn_conf['num_alignments'], blastn_conf['evalue'])
            if blastn_conf['penalty'] != None:
                blastn_line += 'penalty %s' % blastn_conf['penalty']
            if blastn_conf['reward'] != None:
                blastn_line += 'reward %s' % blastn_conf['reward']
            try:
                print('blastn_line = %s' % blastn_line)
                res = os.popen(blastn_line).read()
            except Exception as e:
                print('blastn exception')
                print(e)
        else:
            # 不组装就用blat
            blastn_conf = self.conf['resequencing']['blastn']
            blat_line = """blat %s %s %s -out=blast8""" % (blastn_conf["blast_db_path"].replace("bacteria_db", "bacteria_sequences.fasta"), self.input_file, self.result_dir+"/resequencing/blat_out")
            
            try:
                print('blat_line = %s' % blat_line)
                os.popen(blat_line)
                res = open(self.result_dir+"/resequencing/blat_out", 'r').read()
            except Exception as e:
                print('blat exception')
                print(e)

        # accession_version_list = scripts.handle_blast_xml_result.handle_blast_xml_result_outfmt7(self.result_dir+"/resequencing/ncbi_bacteria_blast_out.xml")
        accession_version_list = scripts.handle_blast_xml_result.handle_blast_xml_result_outfmt7_v2([list(x.split()) for x in res.splitlines()])
        # print("accession_version_list="+str(accession_version_list))
        if int(self.conf['resequencing']['number_of_candidate_similar_genome']) <= len(accession_version_list):
            accession_version_list = [accession_version_list[i][0] for i in range(int(self.conf['resequencing']['number_of_candidate_similar_genome']))]
        else:
            accession_version_list = [accession_version_list[i][0] for i in range(len(accession_version_list))]

        # accession_version_list = []
        # i = 0
        # for accession_version in accession_version_dict:
        #     if i < int(self.conf['resequencing']['number_of_candidate_similar_genome']):
        #         i += 1
        #     else:
        #         break
        #     accession_version_list.append(accession_version[0])
        #     scripts.download_nucleotide.download_by_accession_version(self.result_dir, accession_version[0])

        # closest_accession_version = scripts.snap.snap(self.result_dir, accession_version_list, self.input_file, self.input_file_paired)
        # print("closest_accession_version="+str(closest_accession_version))

        temp_file = open(os.path.abspath('.') + '/' + self.result_dir+'/Summary_of_results.html', 'a+')
        temp_file.write('<ul>\n')
        temp_file.write('<li>candiate closest genome list %s \n' % (accession_version_list))
        temp_file.write('</ul>\n')
        temp_file.close()

        for accession_version in accession_version_list:
            time.sleep(1)
            scripts.download_nucleotide.download_by_accession_version(self.result_dir, accession_version)

        genome_similarity_calculator = calculating_genome_similarity.GenomeSimilarityCalculator(self.result_dir, blast_input, accession_version_list)
        closest_accession_version, max_jaccard = genome_similarity_calculator.search_similar_genome_from_bacteria_dataset()
        print("closest_accession_version = %s, max_jaccard = %s" % (closest_accession_version, max_jaccard))

        temp_file = open(os.path.abspath('.') + '/' + self.result_dir+'/Summary_of_results.html', 'a+')
        temp_file.write('<ul>\n')
        temp_file.write('<li>closest accession version is %s in %s Jaccard similarity coefficient</li>\n' % (closest_accession_version, max_jaccard))
        temp_file.write('</ul>\n')
        temp_file.close()

        alignment_tool = 'bowtie2' if self.conf['resequencing']['alignment_tool']['bowtie2'] != None else "snap"
        if assembly_conf['enable']:
            subprocess.run('bowtie2-build '+self.result_dir+'/resequencing/'+closest_accession_version+'.fasta '+self.result_dir+'/resequencing/'+closest_accession_version+'_bowtie2_index', shell=True, check=True)
        else:
            if alignment_tool == 'bowtie2':
                subprocess.run('bowtie2-build '+self.result_dir+'/resequencing/'+closest_accession_version+'.fasta '+self.result_dir+'/resequencing/'+closest_accession_version+'_bowtie2_index', shell=True, check=True)
            else:
                subprocess.run('snap-aligner index '+self.result_dir+'/resequencing/'+closest_accession_version+'.fasta '+self.result_dir+'/resequencing/'+closest_accession_version+'_index', shell=True, check=True)
        
        if self.input_file_paired == None:
            assembly_conf = self.conf['resequencing']['assembly']
            if assembly_conf['enable']:
                # if alignment_tool == 'bowtie2':
                subprocess.run('bowtie2 -x '+self.result_dir+'/resequencing/'+closest_accession_version+'_bowtie2_index -f -U '+assembly_result+' -S '+self.result_dir+'/resequencing/alignment_result.sam', shell=True, check=True)
                # else:
                #     scripts.fasta2fastq.fasta_to_fastq(assembly_result, self.result_dir)
                #     assembly_result = self.result_dir+"/resequencing/assembly_result.fastq"
                #     subprocess.run('snap-aligner single '+self.result_dir+'/resequencing/'+closest_accession_version+'_index '+assembly_result+' -o '+self.result_dir+'/resequencing/alignment_result.sam', shell=True, check=True)
            else:
                if alignment_tool == 'bowtie2':
                    subprocess.run('bowtie2 -x '+self.result_dir+'/resequencing/'+closest_accession_version+'_bowtie2_index -U '+self.input_file+' -S '+self.result_dir+'/resequencing/alignment_result.sam', shell=True, check=True)
                else:
                    subprocess.run('snap-aligner single '+self.result_dir+'/resequencing/'+closest_accession_version+'_index '+self.input_file+' -o '+self.result_dir+'/resequencing/alignment_result.sam', shell=True, check=True)
            # subprocess.run('bowtie2 -x '+self.result_dir+'/resequencing/'+closest_accession_version+'_bowtie2_index -U '+self.input_file+' -S '+self.result_dir+'/resequencing/alignment_result.sam', shell=True, check=True)
        else:
            assembly_conf = self.conf['resequencing']['assembly']
            if assembly_conf['enable']:
                # if alignment_tool == 'bowtie2':
                subprocess.run('bowtie2 -x '+self.result_dir+'/resequencing/'+closest_accession_version+'_bowtie2_index -f -U '+assembly_result+' -S '+self.result_dir+'/resequencing/alignment_result.sam', shell=True, check=True)
                # else:
                #     scripts.fasta2fastq.fasta_to_fastq(assembly_result, self.result_dir)
                #     assembly_result = self.result_dir+"/resequencing/assembly_result.fastq"
                #     subprocess.run('snap-aligner single '+self.result_dir+'/resequencing/'+closest_accession_version+'_index '+assembly_result+' -o '+self.result_dir+'/resequencing/alignment_result.sam', shell=True, check=True)
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

        temp_file = open(os.path.abspath('.') + '/' + self.result_dir+'/Summary_of_results.html', 'a+')
        temp_file.write('<ul>\n')
        temp_file.write('<li>CDS annotation is in %s</li>\n' % ('resequencing/comparative_annotation'))
        temp_file.write('<li><a href="%s">click to CDS annotation result</a></li>\n' % ('./resequencing/comparative_annotation'))
        temp_file.write('</ul>\n')
        temp_file.close()

        if assembly_conf['enable']:
            self.blastn_and_annotation(assembly_result)
        else:
            self.blastn_and_annotation(self.input_file)

    def run_3gs(self):

        print('Begin Resequencing 3gs')

        assembly_conf = self.conf['resequencing']['assembly']
        if assembly_conf['enable']:
            print('begin assembly')
            os.mkdir(os.path.abspath('.')+'/'+self.result_dir+'/resequencing/assembly')
            assembly_obj = assembly.Assembly(self.result_dir+'/resequencing/assembly', assembly_conf, input_file=self.input_file)
            if assembly_conf['canu']['enable']:
                assembly_obj.canu()
                assembly_result = self.result_dir+'/resequencing/assembly/canu_output/canu_assembly_result.contigs.fasta'
            elif assembly_conf['spades']['enable']:
                assembly_obj.spades_3gs()
                assembly_result = self.result_dir+'/resequencing/assembly/spades_out/contigs.fasta'
            else:
                assembly_obj.canu()
                assembly_result = self.result_dir+'/resequencing/assembly/canu_output/canu_assembly_result.contigs.fasta'
            print('end assembly')

        if assembly_conf['enable']:
            blastn_conf = self.conf['resequencing']['blastn']

            # blastn_cline = NcbiblastnCommandline(cmd=os.path.dirname(os.path.realpath(__file__))+'/external_tools/blastn', query=assembly_result, db=self.conf["resequencing"]['blastn']["blast_db_path"], outfmt=7, out=self.result_dir+"/resequencing/ncbi_bacteria_blast_out.xml")
            # blastn_cline.set_parameter('num_threads', int(blastn_conf['num_threads']))
            # blastn_cline.set_parameter('num_alignments', int(blastn_conf['num_alignments']))
            # if blastn_conf['evalue'] != None:
            #     blastn_cline.set_parameter('evalue', float(blastn_conf['evalue']))
            # # if blastn_conf['task'] != None:
            # #     blastn_cline.set_parameter('task', blastn_conf['task'])
            # if blastn_conf['penalty'] != None:
            #     blastn_cline.set_parameter('penalty', int(blastn_conf['penalty']))
            # if blastn_conf['reward'] != None:
            #     blastn_cline.set_parameter('reward', int(blastn_conf['reward']))
            # print(blastn_cline)
            # stdout, stderr = blastn_cline()
            # print(stdout)
            # print(stderr)

            blastn_line = """%s -query %s -db %s -outfmt 7 -num_threads %s -num_alignments %s -evalue %s""" % ('blastn', assembly_result, self.conf["resequencing"]['blastn']["blast_db_path"], blastn_conf['num_threads'], blastn_conf['num_alignments'], blastn_conf['evalue'])
            if blastn_conf['penalty'] != None:
                blastn_line += 'penalty %s' % blastn_conf['penalty']
            if blastn_conf['reward'] != None:
                blastn_line += 'reward %s' % blastn_conf['reward']
            try:
                print('blastn_line = %s' % blastn_line)
                res = os.popen(blastn_line).read()
            except Exception as e:
                print('blastn exception')
                print(e)
        else:
            # 不组装就用blat
            blastn_conf = self.conf['resequencing']['blastn']
            blat_line = """blat %s %s %s -out=blast8""" % (blastn_conf["blast_db_path"].replace("bacteria_db", "bacteria_sequences.fasta"), self.input_file, self.result_dir+"/resequencing/blat_out")
            
            try:
                print('blat_line = %s' % blat_line)
                os.popen(blat_line)
                res = open(self.result_dir+"/resequencing/blat_out", 'r').read()
            except Exception as e:
                print('blat exception')
                print(e)

        # accession_version_list = scripts.handle_blast_xml_result.handle_blast_xml_result_outfmt7(self.result_dir+"/resequencing/ncbi_bacteria_blast_out.xml")
        accession_version_list = scripts.handle_blast_xml_result.handle_blast_xml_result_outfmt7_v2([list(x.split()) for x in res.splitlines()])
        # print("accession_version_list="+str(accession_version_list))
        if int(self.conf['resequencing']['number_of_candidate_similar_genome']) <= len(accession_version_list):
            accession_version_list = [accession_version_list[i][0] for i in range(int(self.conf['resequencing']['number_of_candidate_similar_genome']))]
        else:
            accession_version_list = [accession_version_list[i][0] for i in range(len(accession_version_list))]

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

        # accession_version_list = []
        # i = 0
        # for accession_version in accession_version_dict:
        #     if i < int(self.conf['resequencing']['number_of_candidate_similar_genome']):
        #         i += 1
        #     else:
        #         break
        #     accession_version_list.append(accession_version[0])
        #     scripts.download_nucleotide.download_by_accession_version(self.result_dir, accession_version[0])

        temp_file = open(os.path.abspath('.') + '/' + self.result_dir+'/Summary_of_results.html', 'a+')
        temp_file.write('<ul>\n')
        temp_file.write('<li>candiate closest genome list %s \n' % (accession_version_list))
        temp_file.write('</ul>\n')
        temp_file.close()

        for accession_version in accession_version_list:
            time.sleep(5)
            scripts.download_nucleotide.download_by_accession_version(self.result_dir, accession_version)

        temp_file = open(os.path.abspath('.') + '/' + self.result_dir+'/Summary_of_results.html', 'a+')
        temp_file.write('<ul>\n')
        temp_file.write('<li>reference genome list: %s </li>\n' % (str(accession_version_list)))
        temp_file.write('</ul>\n')
        temp_file.close()

        genome_similarity_calculator = calculating_genome_similarity.GenomeSimilarityCalculator(self.result_dir, assembly_result, accession_version_list)
        closest_accession_version, max_jaccard = genome_similarity_calculator.search_similar_genome_from_bacteria_dataset()
        print("closest_accession_version = %s, max_jaccard = %s" % (closest_accession_version, max_jaccard))

        temp_file = open(os.path.abspath('.') + '/' + self.result_dir+'/Summary_of_results.html', 'a+')
        temp_file.write('<ul>\n')
        temp_file.write('<li>closest genome (accession version) is %s in %s Jaccard similarity coefficient</li>\n' % (closest_accession_version, max_jaccard))
        temp_file.write('</ul>\n')
        temp_file.close()

        # minimap2
        minimap2_conf = self.conf['resequencing']['minimap2']
        command_line = 'minimap2 -a '
        if minimap2_conf['-H'] != None:
            command_line += "-H "
        if minimap2_conf['-k'] != None:
            command_line += "-k "+str(minimap2_conf['-k'])+' '
        if minimap2_conf['-w'] != None:
            command_line += "-w "+str(minimap2_conf['-w'])+' '
        if minimap2_conf['-I'] != None:
            command_line += "-I "+str(minimap2_conf['-I'])+' '
        if minimap2_conf['-f'] != None:
            command_line += "-f "+str(minimap2_conf['-f'])+' '
        if minimap2_conf['-g'] != None:
            command_line += "-g "+str(minimap2_conf['-g'])+' '
        if minimap2_conf['-G'] != None:
            command_line += "-G "+str(minimap2_conf['-G'])+' '
        if minimap2_conf['-F'] != None:
            command_line += "-F "+str(minimap2_conf['-F'])+' '
        if minimap2_conf['-r'] != None:
            command_line += "-r "+str(minimap2_conf['-r'])+' '
        if minimap2_conf['-n'] != None:
            command_line += "-n "+str(minimap2_conf['-n'])+' '
        if minimap2_conf['-m'] != None:
            command_line += "-m "+str(minimap2_conf['-m'])+' '
        if minimap2_conf['-X'] != None:
            command_line += "-X "
        if minimap2_conf['-p'] != None:
            command_line += "-p "+str(minimap2_conf['-p'])+' '
        if minimap2_conf['-N'] != None:
            command_line += "-N "+str(minimap2_conf['-N'])+' '
        if minimap2_conf['-A'] != None:
            command_line += "-A "+str(minimap2_conf['-A'])+' '
        if minimap2_conf['-B'] != None:
            command_line += "-B "+str(minimap2_conf['-B'])+' '
        if minimap2_conf['-O'] != None:
            command_line += "-O "+str(minimap2_conf['-O'])+' '
        if minimap2_conf['-E'] != None:
            command_line += "-E "+str(minimap2_conf['-E'])+' '
        if minimap2_conf['-z'] != None:
            command_line += "-z "+str(minimap2_conf['-z'])+' '
        if minimap2_conf['-s'] != None:
            command_line += "-s "+str(minimap2_conf['-s'])+' '
        if minimap2_conf['-u'] != None:
            command_line += "-u "+str(minimap2_conf['-u'])+' '
        if minimap2_conf['-x'] != None:
            command_line += "-x "+str(minimap2_conf['-x'])+' '
        if assembly_conf['enable']:
            command_line += self.result_dir+"/resequencing/%s.fasta "+assembly_result+" > "+self.result_dir+"/resequencing/minimap2_%s_result.sam"
        else:
            command_line += self.result_dir+"/resequencing/%s.fasta "+self.input_file+" > "+self.result_dir+"/resequencing/minimap2_%s_result.sam"


        # closest_accession_version = accession_version_list[0]
        # max_similarity = 0
        # for accession_version in accession_version_list:
        #     print(command_line % (accession_version, accession_version))
        #     subprocess.run(command_line % (accession_version, accession_version), shell=True, check=True)
        #     subprocess.run(os.path.dirname(os.path.realpath(__file__))+'/external_tools/samtools-1.10/samtools sort '+self.result_dir+'/resequencing/minimap2_'+accession_version+'_result.sam -o '+self.result_dir+'/resequencing/minimap2_'+accession_version+'_result.sorted.bam -O bam', shell=True, check=True)
        #     res = subprocess.Popen('bedtools genomecov -ibam '+self.result_dir+'/resequencing/minimap2_'+accession_version+'_result.sorted.bam', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
        #     send_messages = res.stdout.readlines()
        #     print(send_messages[0].decode('utf-8'))
        #     for message in send_messages:
        #         message = message.decode('utf-8')
        #         if message.startswith('genome\t0'):
        #             print(message)
        #             line = message.split()
        #             if line[-1] == '\n':
        #                 new_similarity = 1 - float(line[-2])
        #             else:
        #                 new_similarity = 1 - float(line[-1])
        #             break
        #     print('new_similarity='+str(new_similarity))
        #     if new_similarity > max_similarity:
        #         max_similarity = new_similarity
        #         closest_accession_version = accession_version
        # print("closest_accession_version="+str(closest_accession_version))

        subprocess.run(command_line % (closest_accession_version, closest_accession_version), shell=True, check=True)

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

        temp_file = open(os.path.abspath('.') + '/' + self.result_dir+'/Summary_of_results.html', 'a+')
        temp_file.write('<ul>\n')
        temp_file.write('<li>CDS annotation is in %s</li>\n' % ('resequencing/comparative_annotation'))
        temp_file.write('<li><a href="%s">click to CDS annotation result</a></li>\n' % ('./resequencing/comparative_annotation'))
        temp_file.write('</ul>\n')
        temp_file.close()
        
        if assembly_conf['enable']:
            self.blastn_and_annotation(assembly_result)
        else:
            self.blastn_and_annotation(self.input_file)

    def blastn_and_annotation(self, assembly_result):
        # blast patric_amr_db
        if self.conf['resequencing']['patric_blastn']['enable']:
            print('begin blast patric_amr_db')
            patric_blastn_conf = self.conf['resequencing']['patric_blastn']

            blastn_line = """%s -query %s -db %s -outfmt 5 -num_threads %s -num_alignments %s -evalue %s""" % ('blastn', assembly_result, patric_blastn_conf["blast_db_path"], patric_blastn_conf['num_threads'], patric_blastn_conf['num_alignments'], patric_blastn_conf['evalue'])
            if patric_blastn_conf['penalty'] != None:
                blastn_line += 'penalty %s' % patric_blastn_conf['penalty']
            if patric_blastn_conf['reward'] != None:
                blastn_line += 'reward %s' % patric_blastn_conf['reward']
            try:
                print('blastn_line = %s' % blastn_line)
                res = os.popen(blastn_line).read()
            except Exception as e:
                print('blastn exception')
                print(e)

            scripts.handle_blast_xml_result.handle(self.result_dir, res, "patric_annotation")

            print('end blast patric_amr_db')

            temp_file = open(os.path.abspath('.') + '/' + self.result_dir+'/Summary_of_results.html', 'a+')
            temp_file.write('<ul>\n')
            temp_file.write('<li>patric annotation is in %s</li>\n' % ('./resequencing/patric_annotation'))
            temp_file.write('<li><a href="%s">click to patric annotation</a></li>\n' % ('./resequencing/patric_annotation'))
            temp_file.write('</ul>\n')
            temp_file.close()

        if self.conf['resequencing']['card_blastn']['enable']:
            # blast card_prevalence_db
            print('begin card_blastn')
            card_blastn_conf = self.conf['resequencing']['card_blastn']

            blastn_line = """%s -query %s -db %s -outfmt 5 -num_threads %s -num_alignments %s -evalue %s""" % ('blastn', assembly_result, card_blastn_conf["blast_db_path"], card_blastn_conf['num_threads'], card_blastn_conf['num_alignments'], card_blastn_conf['evalue'])
            if card_blastn_conf['penalty'] != None:
                blastn_line += 'penalty %s' % card_blastn_conf['penalty']
            if card_blastn_conf['reward'] != None:
                blastn_line += 'reward %s' % card_blastn_conf['reward']
            try:
                print('blastn_line = %s' % blastn_line)
                res = os.popen(blastn_line).read()
            except Exception as e:
                print('blastn exception')
                print(e)

            scripts.handle_blast_xml_result.handle(self.result_dir, res, "card_annotation")

            print('end card_blastn')

            temp_file = open(os.path.abspath('.') + '/' + self.result_dir+'/Summary_of_results.html', 'a+')
            temp_file.write('<ul>\n')
            temp_file.write('<li>card annotation is in %s</li>\n' % ('./resequencing/card_annotation'))
            temp_file.write('<li><a href="%s">click to card annotation</a></li>\n' % ('./resequencing/card_annotation'))
            temp_file.write('</ul>\n')
            temp_file.close()

        if self.conf['resequencing']['drugbank_blastn']['enable']:
            # blast drugbank_db
            print('begin drugbank_db')
            drugbank_blastn_conf = self.conf['resequencing']['card_blastn']

            blastn_line = """%s -query %s -db %s -outfmt 5 -num_threads %s -num_alignments %s -evalue %s""" % ('blastn', assembly_result, drugbank_blastn_conf["blast_db_path"], drugbank_blastn_conf['num_threads'], drugbank_blastn_conf['num_alignments'], drugbank_blastn_conf['evalue'])
            if drugbank_blastn_conf['penalty'] != None:
                blastn_line += 'penalty %s' % drugbank_blastn_conf['penalty']
            if drugbank_blastn_conf['reward'] != None:
                blastn_line += 'reward %s' % drugbank_blastn_conf['reward']
            try:
                print('blastn_line = %s' % blastn_line)
                res = os.popen(blastn_line).read()
            except Exception as e:
                print('blastn exception')
                print(e)

            scripts.handle_denovo_blast_result.handle(self.result_dir, res, "drugbank_annotation")
        
            print('end drugbank_db')
            
            temp_file = open(os.path.abspath('.') + '/' + self.result_dir+'/Summary_of_results.html', 'a+')
            temp_file.write('<ul>\n')
            temp_file.write('<li>drugbank annotation is in %s</li>\n' % ('./resequencing/drugbank_annotation'))
            temp_file.write('<li><a href="%s">click to drugbank annotation</a></li>\n' % ('./resequencing/drugbank_annotation'))
            temp_file.write('</ul>\n')
            temp_file.close()
