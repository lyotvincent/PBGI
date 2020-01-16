import subprocess
import os

class Assembly:

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

    def megahit_single(self):
        print("begin megahit_single_end")
        megahit_conf = self.conf['denovo']['assembly']['megahit']
        command_line = 'megahit '
        if not megahit_conf['--min-count'] == 0:
            command_line += '--min-count '+str(megahit_conf['--min-count'])+' '
        if not megahit_conf['--k-list'] == 0:
            command_line += '--k-list '+str(megahit_conf['--k-list'])+' '
        if not megahit_conf['--no-mercy'] == 0:
            command_line += '--no-mercy '
        if not megahit_conf['--bubble-level'] == 0:
            command_line += '--bubble-level '+str(megahit_conf['--bubble-level'])+' '
        if not megahit_conf['--merge-level'] == 0:
            command_line += '--merge-level '+str(megahit_conf['--merge-level'])+' '
        if not megahit_conf['--prune-level'] == 0:
            command_line += '--prune-level '+str(megahit_conf['--prune-level'])+' '
        if not megahit_conf['--prune-depth'] == 0:
            command_line += '--prune-depth '+str(megahit_conf['--prune-depth'])+' '
        if not megahit_conf['--low-local-ratio'] == 0:
            command_line += '--low-local-ratio '+str(megahit_conf['--low-local-ratio'])+' '
        if not megahit_conf['--max-tip-len'] == 0:
            command_line += '--max-tip-len '+str(megahit_conf['--max-tip-len'])+' '
        if not megahit_conf['--no-local'] == 0:
            command_line += '--no-local '
        if not megahit_conf['--kmin-1pass'] == 0:
            command_line += '--kmin-1pass '
        if not megahit_conf['-m'] == 0:
            command_line += '-m '+str(megahit_conf['-m'])+' '
        if not megahit_conf['--mem-flag'] == 0:
            command_line += '--mem-flag '+str(megahit_conf['--mem-flag'])+' '
        if not megahit_conf['-t'] == 0:
            command_line += '-t '+str(megahit_conf['-t'])+' '
        if not megahit_conf['--no-hw-accel'] == 0:
            command_line += '--no-hw-accel '
        if not megahit_conf['--min-contig-len'] == 0:
            command_line += '--min-contig-len '+str(megahit_conf['--min-contig-len'])+' '
        command_line += '-r '+self.input_file+' -o '+self.result_dir+'/megahit_out'
        subprocess.run(command_line, shell=True, check=True)
        print("end megahit_single_end")
    
    def megahit_paired(self):
        print("begin megahit_paired_end")
        megahit_conf = self.conf['denovo']['assembly']['megahit']
        command_line = 'megahit '
        if not megahit_conf['--min-count'] == 0:
            command_line += '--min-count '+str(megahit_conf['--min-count'])+' '
        if not megahit_conf['--k-list'] == 0:
            command_line += '--k-list '+str(megahit_conf['--k-list'])+' '
        if not megahit_conf['--no-mercy'] == 0:
            command_line += '--no-mercy '
        if not megahit_conf['--bubble-level'] == 0:
            command_line += '--bubble-level '+str(megahit_conf['--bubble-level'])+' '
        if not megahit_conf['--merge-level'] == 0:
            command_line += '--merge-level '+str(megahit_conf['--merge-level'])+' '
        if not megahit_conf['--prune-level'] == 0:
            command_line += '--prune-level '+str(megahit_conf['--prune-level'])+' '
        if not megahit_conf['--prune-depth'] == 0:
            command_line += '--prune-depth '+str(megahit_conf['--prune-depth'])+' '
        if not megahit_conf['--low-local-ratio'] == 0:
            command_line += '--low-local-ratio '+str(megahit_conf['--low-local-ratio'])+' '
        if not megahit_conf['--max-tip-len'] == 0:
            command_line += '--max-tip-len '+str(megahit_conf['--max-tip-len'])+' '
        if not megahit_conf['--no-local'] == 0:
            command_line += '--no-local '
        if not megahit_conf['--kmin-1pass'] == 0:
            command_line += '--kmin-1pass '
        if not megahit_conf['-m'] == 0:
            command_line += '-m '+str(megahit_conf['-m'])+' '
        if not megahit_conf['--mem-flag'] == 0:
            command_line += '--mem-flag '+str(megahit_conf['--mem-flag'])+' '
        if not megahit_conf['-t'] == 0:
            command_line += '-t '+str(megahit_conf['-t'])+' '
        if not megahit_conf['--no-hw-accel'] == 0:
            command_line += '--no-hw-accel '
        if not megahit_conf['--min-contig-len'] == 0:
            command_line += '--min-contig-len '+str(megahit_conf['--min-contig-len'])+' '
        command_line += '-1 '+self.input_file_1+' -2 '+self.input_file_2+' -o '+self.result_dir+'/megahit_out'
        subprocess.run(command_line, shell=True, check=True)
        print("end megahit_paired_end")

    def megahit_interleaved(self):
        print("begin megahit_interleaved")
        megahit_conf = self.conf['denovo']['assembly']['megahit']
        command_line = 'megahit '
        if not megahit_conf['--min-count'] == 0:
            command_line += '--min-count '+str(megahit_conf['--min-count'])+' '
        if not megahit_conf['--k-list'] == 0:
            command_line += '--k-list '+str(megahit_conf['--k-list'])+' '
        if not megahit_conf['--no-mercy'] == 0:
            command_line += '--no-mercy '
        if not megahit_conf['--bubble-level'] == 0:
            command_line += '--bubble-level '+str(megahit_conf['--bubble-level'])+' '
        if not megahit_conf['--merge-level'] == 0:
            command_line += '--merge-level '+str(megahit_conf['--merge-level'])+' '
        if not megahit_conf['--prune-level'] == 0:
            command_line += '--prune-level '+str(megahit_conf['--prune-level'])+' '
        if not megahit_conf['--prune-depth'] == 0:
            command_line += '--prune-depth '+str(megahit_conf['--prune-depth'])+' '
        if not megahit_conf['--low-local-ratio'] == 0:
            command_line += '--low-local-ratio '+str(megahit_conf['--low-local-ratio'])+' '
        if not megahit_conf['--max-tip-len'] == 0:
            command_line += '--max-tip-len '+str(megahit_conf['--max-tip-len'])+' '
        if not megahit_conf['--no-local'] == 0:
            command_line += '--no-local '
        if not megahit_conf['--kmin-1pass'] == 0:
            command_line += '--kmin-1pass '
        if not megahit_conf['-m'] == 0:
            command_line += '-m '+str(megahit_conf['-m'])+' '
        if not megahit_conf['--mem-flag'] == 0:
            command_line += '--mem-flag '+str(megahit_conf['--mem-flag'])+' '
        if not megahit_conf['-t'] == 0:
            command_line += '-t '+str(megahit_conf['-t'])+' '
        if not megahit_conf['--no-hw-accel'] == 0:
            command_line += '--no-hw-accel '
        if not megahit_conf['--min-contig-len'] == 0:
            command_line += '--min-contig-len '+str(megahit_conf['--min-contig-len'])+' '
        command_line += '-12 '+self.input_file_12+' -o '+self.result_dir+'/megahit_out'
        subprocess.run(command_line, shell=True, check=True)
        print("end megahit_interleaved")

    def spades_single(self):
        print("begin spades_single_end")
        spades_conf = self.conf['denovo']['assembly']['spades']
        command_line = 'spades.py -s '+self.input_file+' '
        if not spades_conf['--iontorrent'] == 0:
            command_line += '--iontorrent '
        if not spades_conf['-t'] == 0:
            command_line += '-t '+str(spades_conf['-t'])+' '
        if not spades_conf['-m'] == 0:
            command_line += '-m '+str(spades_conf['-m'])+' '
        if not spades_conf['-k'] == 0:
            command_line += '-k '+str(spades_conf['-k'])+' '
        if not spades_conf['--cov-cutoff'] == 0:
            command_line += '--cov-cutoff '+str(spades_conf['--cov-cutoff'])+' '
        if not spades_conf['--phred-offset'] == 0:
            command_line += '--phred-offset '+str(spades_conf['--phred-offset'])+' '
        command_line += '-o '+self.result_dir+'/spades_out/'
        subprocess.run(command_line, shell=True, check=True)
        print("end spades_single_end")
    
    def spades_paired(self):
        print("begin spades_paired_end")
        spades_conf = self.conf['denovo']['assembly']['spades']
        command_line = 'spades.py -1 '+self.input_file_1+' -2 '+self.input_file_2+' '
        if not spades_conf['--iontorrent'] == 0:
            command_line += '--iontorrent '
        if not spades_conf['-t'] == 0:
            command_line += '-t '+str(spades_conf['-t'])+' '
        if not spades_conf['-m'] == 0:
            command_line += '-m '+str(spades_conf['-m'])+' '
        if not spades_conf['-k'] == 0:
            command_line += '-k '+str(spades_conf['-k'])+' '
        if not spades_conf['--cov-cutoff'] == 0:
            command_line += '--cov-cutoff '+str(spades_conf['--cov-cutoff'])+' '
        if not spades_conf['--phred-offset'] == 0:
            command_line += '--phred-offset '+str(spades_conf['--phred-offset'])+' '
        command_line = '-o '+self.result_dir+'/spades_out/'
        subprocess.run(command_line, shell=True, check=True)
        print("end spades_paired_end")
    
    def spades_interlaced(self):
        print("begin spades_interlaced")
        spades_conf = self.conf['denovo']['assembly']['spades']
        command_line = 'spades.py -12 '+self.input_file_12+' '
        if not spades_conf['--iontorrent'] == 0:
            command_line += '--iontorrent '
        if not spades_conf['-t'] == 0:
            command_line += '-t '+str(spades_conf['-t'])+' '
        if not spades_conf['-m'] == 0:
            command_line += '-m '+str(spades_conf['-m'])+' '
        if not spades_conf['-k'] == 0:
            command_line += '-k '+str(spades_conf['-k'])+' '
        if not spades_conf['--cov-cutoff'] == 0:
            command_line += '--cov-cutoff '+str(spades_conf['--cov-cutoff'])+' '
        if not spades_conf['--phred-offset'] == 0:
            command_line += '--phred-offset '+str(spades_conf['--phred-offset'])+' '
        command_line += '-o '+self.result_dir+'/spades_out/'
        subprocess.run(command_line, shell=True, check=True)
        print("end spades_interlaced")
    
    def spades_3gs(self):
        print("begin spades_single_end")
        spades_conf = self.conf['denovo']['assembly']['spades']
        if not spades_conf['--pacbio'] == 0:
            command_line = 'spades.py --pacbio '+self.input_file+' '
        elif not spades_conf['--nanopore'] == 0:
            command_line = 'spades.py --nanopore '+self.input_file+' '
        else:
            print('please choose pacbio or nanopore in conf.json. use pacbio by default.')
            command_line = 'spades.py --pacbio '+self.input_file+' '
        if not spades_conf['-t'] == 0:
            command_line += '-t '+str(spades_conf['-t'])+' '
        if not spades_conf['-m'] == 0:
            command_line += '-m '+str(spades_conf['-m'])+' '
        if not spades_conf['-k'] == 0:
            command_line += '-k '+str(spades_conf['-k'])+' '
        if not spades_conf['--cov-cutoff'] == 0:
            command_line += '--cov-cutoff '+str(spades_conf['--cov-cutoff'])+' '
        if not spades_conf['--phred-offset'] == 0:
            command_line += '--phred-offset '+str(spades_conf['--phred-offset'])+' '
        command_line += '-o '+self.result_dir+'/spades_out/'
        subprocess.run(command_line, shell=True, check=True)
        print("end spades_single_end")

    def velvet_single(self):
        print('begin velvet_single_end')
        velveth_conf = self.conf['denovo']['assembly']['velvet']['velveth']
        command_line = 'velveth '+self.result_dir+'/velvet_output/ '
        if not velveth_conf['hash_length'] == 0:
            command_line += str(velveth_conf['hash_length'])+' -fastq '
        else:
            command_line += "31 -fastq "
        if not velveth_conf['-short'] == 0:
            command_line += "-short "
        elif not velveth_conf['-short2'] == 0:
            command_line += '-short2 '
        elif not velveth_conf['-long'] == 0:
            command_line += '-long'
        else:
            command_line += "-short "
        command_line += self.input_file
        subprocess.run(command_line, shell=True, check=True)
        
        velvetg_conf = self.conf['denovo']['assembly']['velvet']['velvetg']
        command_line = 'velvetg '+self.result_dir+'/velvet_output/ '
        if not velvetg_conf['-cov_cutoff'] == 0:
            command_line += '-cov_cutoff '+str(velvetg_conf['-cov_cutoff'])+' '
        if not velvetg_conf['-ins_length'] == 0:
            command_line += '-ins_length '+str(velvetg_conf['-ins_length'])+' '
        if not velvetg_conf['-read_trkg'] == 0:
            command_line += '-read_trkg '+str(velvetg_conf['-read_trkg'])+' '
        if not velvetg_conf['-min_contig_lgth'] == 0:
            command_line += '-min_contig_lgth '+str(velvetg_conf['-min_contig_lgth'])+' '
        if not velvetg_conf['-amos_file'] == 0:
            command_line += '-amos_file '+str(velvetg_conf['-amos_file'])+' '
        if not velvetg_conf['-exp_cov'] == 0:
            command_line += '-exp_cov '+str(velvetg_conf['-exp_cov'])+' '
        if not velvetg_conf['-long_cov_cutoff'] == 0:
            command_line += '-long_cov_cutoff '+str(velvetg_conf['-long_cov_cutoff'])+' '
        if not velvetg_conf['-ins_length'] == 0:
            command_line += '-ins_length '+str(velvetg_conf['-ins_length'])+' '
        if not velvetg_conf['-ins_length_long'] == 0:
            command_line += '-ins_length_long '+str(velvetg_conf['-ins_length_long'])+' '
        if not velvetg_conf['-ins_length*_sd'] == 0:
            command_line += '-ins_length*_sd '+str(velvetg_conf['-ins_length*_sd'])+' '
        if not velvetg_conf['-scaffolding'] == 0:
            command_line += '-scaffolding '+str(velvetg_conf['-scaffolding'])+' '
        if not velvetg_conf['-max_branch_length'] == 0:
            command_line += '-max_branch_length '+str(velvetg_conf['-max_branch_length'])+' '
        if not velvetg_conf['-max_divergence'] == 0:
            command_line += '-max_divergence '+str(velvetg_conf['-max_divergence'])+' '
        if not velvetg_conf['-max_gap_count'] == 0:
            command_line += '-max_gap_count '+str(velvetg_conf['-max_gap_count'])+' '
        if not velvetg_conf['-min_pair_count'] == 0:
            command_line += '-min_pair_count '+str(velvetg_conf['-min_pair_count'])+' '
        if not velvetg_conf['-max_coverage'] == 0:
            command_line += '-max_coverage '+str(velvetg_conf['-max_coverage'])+' '
        if not velvetg_conf['-coverage_mask'] == 0:
            command_line += '-coverage_mask '+str(velvetg_conf['-coverage_mask'])+' '
        if not velvetg_conf['-long_mult_cutoff'] == 0:
            command_line += '-long_mult_cutoff '+str(velvetg_conf['-long_mult_cutoff'])+' '
        if not velvetg_conf['-alignments'] == 0:
            command_line += '-alignments '+str(velvetg_conf['-alignments'])+' '
        if not velvetg_conf['-exportFiltered'] == 0:
            command_line += '-exportFiltered '+str(velvetg_conf['-exportFiltered'])+' '
        if not velvetg_conf['-clean'] == 0:
            command_line += '-clean '+str(velvetg_conf['-clean'])+' '
        if not velvetg_conf['-very_clean'] == 0:
            command_line += '-very_clean '+str(velvetg_conf['-very_clean'])+' '
        if not velvetg_conf['-paired_exp_fraction'] == 0:
            command_line += '-paired_exp_fraction '+str(velvetg_conf['-paired_exp_fraction'])+' '
        if not velvetg_conf['-conserveLong'] == 0:
            command_line += '-conserveLong '+str(velvetg_conf['-conserveLong'])
        subprocess.run(command_line, shell=True, check=True)
        print('end velvet_single_end')
    
    def velvet_paired(self):
        print('begin velvet_paired_end')
        velveth_conf = self.conf['denovo']['assembly']['velvet']['velveth']
        command_line = 'velveth '+self.result_dir+'/velvet_output/ '
        if not velveth_conf['hash_length'] == 0:
            command_line += str(velveth_conf['hash_length'])+' -fastq '
        else:
            command_line += "31 -fastq "
        if not velveth_conf['-short'] == 0:
            command_line += "-short "
        elif not velveth_conf['-short2'] == 0:
            command_line += '-short2 '
        elif not velveth_conf['-long'] == 0:
            command_line += '-long'
        else:
            command_line += "-short "
        command_line += '-separate '
        command_line += self.input_file_1+' '+self.input_file_2
        subprocess.run(command_line, shell=True, check=True)
        
        velvetg_conf = self.conf['denovo']['assembly']['velvet']['velvetg']
        command_line = 'velvetg '+self.result_dir+'/velvet_output/ '
        if not velvetg_conf['-cov_cutoff'] == 0:
            command_line += '-cov_cutoff '+str(velvetg_conf['-cov_cutoff'])+' '
        if not velvetg_conf['-ins_length'] == 0:
            command_line += '-ins_length '+str(velvetg_conf['-ins_length'])+' '
        if not velvetg_conf['-read_trkg'] == 0:
            command_line += '-read_trkg '+str(velvetg_conf['-read_trkg'])+' '
        if not velvetg_conf['-min_contig_lgth'] == 0:
            command_line += '-min_contig_lgth '+str(velvetg_conf['-min_contig_lgth'])+' '
        if not velvetg_conf['-amos_file'] == 0:
            command_line += '-amos_file '+str(velvetg_conf['-amos_file'])+' '
        if not velvetg_conf['-exp_cov'] == 0:
            command_line += '-exp_cov '+str(velvetg_conf['-exp_cov'])+' '
        if not velvetg_conf['-long_cov_cutoff'] == 0:
            command_line += '-long_cov_cutoff '+str(velvetg_conf['-long_cov_cutoff'])+' '
        if not velvetg_conf['-ins_length'] == 0:
            command_line += '-ins_length '+str(velvetg_conf['-ins_length'])+' '
        if not velvetg_conf['-ins_length_long'] == 0:
            command_line += '-ins_length_long '+str(velvetg_conf['-ins_length_long'])+' '
        if not velvetg_conf['-ins_length*_sd'] == 0:
            command_line += '-ins_length*_sd '+str(velvetg_conf['-ins_length*_sd'])+' '
        if not velvetg_conf['-scaffolding'] == 0:
            command_line += '-scaffolding '+str(velvetg_conf['-scaffolding'])+' '
        if not velvetg_conf['-max_branch_length'] == 0:
            command_line += '-max_branch_length '+str(velvetg_conf['-max_branch_length'])+' '
        if not velvetg_conf['-max_divergence'] == 0:
            command_line += '-max_divergence '+str(velvetg_conf['-max_divergence'])+' '
        if not velvetg_conf['-max_gap_count'] == 0:
            command_line += '-max_gap_count '+str(velvetg_conf['-max_gap_count'])+' '
        if not velvetg_conf['-min_pair_count'] == 0:
            command_line += '-min_pair_count '+str(velvetg_conf['-min_pair_count'])+' '
        if not velvetg_conf['-max_coverage'] == 0:
            command_line += '-max_coverage '+str(velvetg_conf['-max_coverage'])+' '
        if not velvetg_conf['-coverage_mask'] == 0:
            command_line += '-coverage_mask '+str(velvetg_conf['-coverage_mask'])+' '
        if not velvetg_conf['-long_mult_cutoff'] == 0:
            command_line += '-long_mult_cutoff '+str(velvetg_conf['-long_mult_cutoff'])+' '
        if not velvetg_conf['-alignments'] == 0:
            command_line += '-alignments '+str(velvetg_conf['-alignments'])+' '
        if not velvetg_conf['-exportFiltered'] == 0:
            command_line += '-exportFiltered '+str(velvetg_conf['-exportFiltered'])+' '
        if not velvetg_conf['-clean'] == 0:
            command_line += '-clean '+str(velvetg_conf['-clean'])+' '
        if not velvetg_conf['-very_clean'] == 0:
            command_line += '-very_clean '+str(velvetg_conf['-very_clean'])+' '
        if not velvetg_conf['-paired_exp_fraction'] == 0:
            command_line += '-paired_exp_fraction '+str(velvetg_conf['-paired_exp_fraction'])+' '
        if not velvetg_conf['-conserveLong'] == 0:
            command_line += '-conserveLong '+str(velvetg_conf['-conserveLong'])
        subprocess.run(command_line, shell=True, check=True)
        print('end velvet_paired_end')
    
    def velvet_interleaved(self):
        print('begin velvet_interleaved')
        velveth_conf = self.conf['denovo']['assembly']['velvet']['velveth']
        command_line = 'velveth '+self.result_dir+'/velvet_output/ '
        if not velveth_conf['hash_length'] == 0:
            command_line += str(velveth_conf['hash_length'])+' -fastq '
        else:
            command_line += "31 -fastq "
        if not velveth_conf['-short'] == 0:
            command_line += "-short "
        elif not velveth_conf['-short2'] == 0:
            command_line += '-short2 '
        elif not velveth_conf['-long'] == 0:
            command_line += '-long'
        else:
            command_line += "-short "
        command_line += '-interleaved '
        command_line += self.input_file_12
        subprocess.run(command_line, shell=True, check=True)

        velvetg_conf = self.conf['denovo']['assembly']['velvet']['velvetg']
        command_line = 'velvetg '+self.result_dir+'/velvet_output/ '
        if not velvetg_conf['-cov_cutoff'] == 0:
            command_line += '-cov_cutoff '+str(velvetg_conf['-cov_cutoff'])+' '
        if not velvetg_conf['-ins_length'] == 0:
            command_line += '-ins_length '+str(velvetg_conf['-ins_length'])+' '
        if not velvetg_conf['-read_trkg'] == 0:
            command_line += '-read_trkg '+str(velvetg_conf['-read_trkg'])+' '
        if not velvetg_conf['-min_contig_lgth'] == 0:
            command_line += '-min_contig_lgth '+str(velvetg_conf['-min_contig_lgth'])+' '
        if not velvetg_conf['-amos_file'] == 0:
            command_line += '-amos_file '+str(velvetg_conf['-amos_file'])+' '
        if not velvetg_conf['-exp_cov'] == 0:
            command_line += '-exp_cov '+str(velvetg_conf['-exp_cov'])+' '
        if not velvetg_conf['-long_cov_cutoff'] == 0:
            command_line += '-long_cov_cutoff '+str(velvetg_conf['-long_cov_cutoff'])+' '
        if not velvetg_conf['-ins_length'] == 0:
            command_line += '-ins_length '+str(velvetg_conf['-ins_length'])+' '
        if not velvetg_conf['-ins_length_long'] == 0:
            command_line += '-ins_length_long '+str(velvetg_conf['-ins_length_long'])+' '
        if not velvetg_conf['-ins_length*_sd'] == 0:
            command_line += '-ins_length*_sd '+str(velvetg_conf['-ins_length*_sd'])+' '
        if not velvetg_conf['-scaffolding'] == 0:
            command_line += '-scaffolding '+str(velvetg_conf['-scaffolding'])+' '
        if not velvetg_conf['-max_branch_length'] == 0:
            command_line += '-max_branch_length '+str(velvetg_conf['-max_branch_length'])+' '
        if not velvetg_conf['-max_divergence'] == 0:
            command_line += '-max_divergence '+str(velvetg_conf['-max_divergence'])+' '
        if not velvetg_conf['-max_gap_count'] == 0:
            command_line += '-max_gap_count '+str(velvetg_conf['-max_gap_count'])+' '
        if not velvetg_conf['-min_pair_count'] == 0:
            command_line += '-min_pair_count '+str(velvetg_conf['-min_pair_count'])+' '
        if not velvetg_conf['-max_coverage'] == 0:
            command_line += '-max_coverage '+str(velvetg_conf['-max_coverage'])+' '
        if not velvetg_conf['-coverage_mask'] == 0:
            command_line += '-coverage_mask '+str(velvetg_conf['-coverage_mask'])+' '
        if not velvetg_conf['-long_mult_cutoff'] == 0:
            command_line += '-long_mult_cutoff '+str(velvetg_conf['-long_mult_cutoff'])+' '
        if not velvetg_conf['-alignments'] == 0:
            command_line += '-alignments '+str(velvetg_conf['-alignments'])+' '
        if not velvetg_conf['-exportFiltered'] == 0:
            command_line += '-exportFiltered '+str(velvetg_conf['-exportFiltered'])+' '
        if not velvetg_conf['-clean'] == 0:
            command_line += '-clean '+str(velvetg_conf['-clean'])+' '
        if not velvetg_conf['-very_clean'] == 0:
            command_line += '-very_clean '+str(velvetg_conf['-very_clean'])+' '
        if not velvetg_conf['-paired_exp_fraction'] == 0:
            command_line += '-paired_exp_fraction '+str(velvetg_conf['-paired_exp_fraction'])+' '
        if not velvetg_conf['-conserveLong'] == 0:
            command_line += '-conserveLong '+str(velvetg_conf['-conserveLong'])
        subprocess.run(command_line, shell=True, check=True)
        print('end velvet_interleaved')
    
    def canu(self):
        print('begin canu')
        canu_conf = self.conf['denovo']['assembly']['canu']
        command_line = 'canu -p canu_assembly_result -d '+self.result_dir+'/canu_output/ '
        if not canu_conf['genomeSize='] == 0:
            command_line += 'genomeSize='+canu_conf['-p']+' '
        else:
            command_line += 'genomeSize=4.8m '
        if not canu_conf['-pacbio-raw'] == 0:
            command_line += '-pacbio-raw '
        elif not canu_conf['-pacbio-corrected'] == 0:
            command_line += '-pacbio-corrected '
        elif not canu_conf['-nanopore-raw'] == 0:
            command_line += '-nanopore-raw '
        elif not canu_conf['-nanopore-corrected'] == 0:
            command_line += '-nanopore-corrected '
        else:
            command_line += '-pacbio-raw '
        command_line += self.input_file
        subprocess.run(command_line, shell=True, check=True)
        print('end canu')
    
def assembly_single_end(input_file, result_dir):
    print("Begin assembly")
    print("input_file="+input_file)
    # result = subprocess.run('megahit -r '+input_file+' -o '+result_dir+'/assembly', shell=True, check=True)
    print("Begin megahit")
    subprocess.run('megahit -r '+input_file+' -o '+result_dir+'/assembly/megahit_out', shell=True, check=True)
    # print("result of megahit:")
    # print(result.stdout)
    print("End megahit")
    print("Begin QUAST")
    subprocess.run('quast.py '+result_dir+'/assembly/megahit_out/final.contigs.fa -o '+result_dir+'/assembly/quast_out', shell=True, check=True)
    print("End QUAST")
    print("End assembly")

def assembly_paired_end(input_file_1, input_file_2, result_dir):
    print("Begin assembly")
    print("input_file="+input_file_1+input_file_2)
    # result = subprocess.run('megahit -r '+input_file+' -o '+result_dir+'/assembly', shell=True, check=True)
    print("Begin megahit")
    subprocess.run('megahit -1 '+input_file_1+' -2 '+input_file_2+' -o '+result_dir+'/assembly/megahit_out', shell=True, check=True)
    # print("result of megahit:")
    # print(result.stdout)
    print("End megahit")
    print("Begin QUAST")
    subprocess.run('quast.py '+result_dir+'/assembly/megahit_out/final.contigs.fa -o '+result_dir+'/assembly/quast_out', shell=True, check=True)
    print("End QUAST")
    print("End assembly")
