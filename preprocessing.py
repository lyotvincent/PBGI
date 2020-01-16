import subprocess

class Preprocessing:

    result_dir = None
    conf = None
    input_file = None
    input_file_1 = None
    input_file_2 = None

    def __init__(self, result_dir, conf, input_file=None, input_file_1=None, input_file_2=None):
        self.result_dir = result_dir
        self.conf = conf
        self.input_file = input_file
        self.input_file_1 = input_file_1
        self.input_file_2 = input_file_2

    def fastqc_single_end(self, folder_name):
        print("begin fastqc single end")
        fastqc_conf = self.conf['preprocessing']['fastqc']
        com = 'fastqc -o '+self.result_dir+'/preprocessing/'+folder_name+' '
        if not fastqc_conf['--casava'] == 0:
            com += "--casava "
        if not fastqc_conf['--nofilter'] == 0:
            com += "--nofilter "
        if not fastqc_conf['--nogroup'] == 0:
            com += "--nogroup "
        if not fastqc_conf['-t'] == 0:
            com += "-t "+fastqc_conf['-t']+' '
        if not fastqc_conf['-c'] == 0:
            com += "-c "+fastqc_conf['-c']+' '
        if not fastqc_conf['-a'] == 0:
            com += "-a "+fastqc_conf['-a']+' '
        if not fastqc_conf['-l'] == 0:
            com += "-l "+fastqc_conf['-l']+' '
        if not fastqc_conf['-k'] == 0:
            com += "-k "+str(fastqc_conf['-t'])+' '
        com += self.input_file
        try:
            subprocess.run(com, shell=True, check=True)
        except Exception as e:
            print(e)
        print("end fastqc single end")
    
    def fastqc_paired_end(self, folder_name):
        print("begin fastqc paired end")
        fastqc_conf = self.conf['preprocessing']['fastqc']
        com = 'fastqc -o '+self.result_dir+'/preprocessing/'+folder_name+' '
        if not fastqc_conf['--casava'] == 0:
            com += "--casava "
        if not fastqc_conf['--nofilter'] == 0:
            com += "--nofilter "
        if not fastqc_conf['--nogroup'] == 0:
            com += "--nogroup "
        if not fastqc_conf['-t'] == 0:
            com += "-t "+fastqc_conf['-t']+' '
        if not fastqc_conf['-c'] == 0:
            com += "-c "+fastqc_conf['-c']+' '
        if not fastqc_conf['-a'] == 0:
            com += "-a "+fastqc_conf['-a']+' '
        if not fastqc_conf['-l'] == 0:
            com += "-l "+fastqc_conf['-l']+' '
        if not fastqc_conf['-k'] == 0:
            com += "-k "+str(fastqc_conf['-t'])+' '
        com += self.input_file_1+" "+self.input_file_2
        try:
            subprocess.run(com, shell=True, check=True)
        except Exception as e:
            print(e)
        print("end fastqc paired end")
    
    def fastp_single_end(self):
        print("begin fastp_single_end")
        fastp_conf = self.conf['preprocessing']['fastp']
        com = 'fastp -i '+self.input_file+' -o '+self.result_dir+'/preprocessing/fastp_output.fastq -j '+self.result_dir+'/preprocessing/fastp.json -h '+self.result_dir+'/preprocessing/fastp.html '
        if not fastp_conf['-6'] == 0:
            com += '-6 '
        if not fastp_conf['--interleaved_in'] == 0:
            com += '--interleaved_in '
        if not fastp_conf['-V'] == 0:
            com += '-V '
        if not fastp_conf['-A'] == 0:
            com += '-A '
        if not fastp_conf['--adapter_sequence'] == 0:
            com += '--adapter_sequence '+fastp_conf['--adapter_sequence']+' '
        # if not fastp_conf['--adapter_sequence_r2'] == 0:
        #     com += '--adapter_sequence_r2 '+fastp_conf['--adapter_sequence_r2']+' '
        if not fastp_conf['--adapter_fasta'] == 0:
            com += '--adapter_fasta '+fastp_conf['--adapter_fasta']+' '
        # if not fastp_conf['--detect_adapter_for_pe'] == 0:
        #     com += '--detect_adapter_for_pe '
        if not fastp_conf['-f'] == 0:
            com += '-f '+str(fastp_conf['-f'])+' '
        if not fastp_conf['-t'] == 0:
            com += '-t '+str(fastp_conf['-t'])+' '
        if not fastp_conf['-b'] == 0:
            com += '-b '+str(fastp_conf['-b'])+' '
        # if not fastp_conf['-F'] == 0:
        #     com += '-F '+str(fastp_conf['-F'])+' '
        # if not fastp_conf['-T'] == 0:
        #     com += '-T '+str(fastp_conf['-T'])+' '
        # if not fastp_conf['-B'] == 0:
        #     com += '-B '+str(fastp_conf['-B'])+' '
        if not fastp_conf['--trim_poly_g'] == 0:
            com += '--trim_poly_g '
        if not fastp_conf['--poly_g_min_len'] == 0:
            com += '--poly_g_min_len '+str(fastp_conf['--poly_g_min_len'])+' '
        if not fastp_conf['-G'] == 0:
            com += '-G '
        if not fastp_conf['--trim_poly_x'] == 0:
            com += '--trim_poly_x '
        if not fastp_conf['--poly_x_min_len'] == 0:
            com += '--poly_x_min_len '+str(fastp_conf['--poly_x_min_len'])+' '
        if not fastp_conf['-5'] == 0:
            com += '-5 '
        if not fastp_conf['-3'] == 0:
            com += '-3 '
        if not fastp_conf['-r'] == 0:
            com += '-r '
        if not fastp_conf['-W'] == 0:
            com += "-W "+str(fastp_conf['-W'])+' '
        if not fastp_conf['--cut_mean_quality'] == 0:
            com += "--cut_mean_quality "+str(fastp_conf['--cut_mean_quality'])+' '
        if not fastp_conf['--cut_front_window_size'] == 0:
            com += "--cut_front_window_size "+str(fastp_conf['--cut_front_window_size'])+' '
        if not fastp_conf['--cut_front_mean_quality'] == 0:
            com += "--cut_front_mean_quality "+str(fastp_conf['--cut_front_mean_quality'])+' '
        if not fastp_conf['--cut_tail_window_size'] == 0:
            com += "--cut_tail_window_size "+str(fastp_conf['--cut_tail_window_size'])+' '
        if not fastp_conf['--cut_tail_mean_quality'] == 0:
            com += "--cut_tail_mean_quality "+str(fastp_conf['--cut_tail_mean_quality'])+' '
        if not fastp_conf['--cut_right_window_size'] == 0:
            com += "--cut_right_window_size "+str(fastp_conf['--cut_right_window_size'])+' '
        if not fastp_conf['--cut_right_mean_quality'] == 0:
            com += "--cut_right_mean_quality "+str(fastp_conf['--cut_right_mean_quality'])+' '
        if not fastp_conf['-Q'] == 0:
            com += '-Q '
        if not fastp_conf['-q'] == 0:
            com += "-q "+str(fastp_conf['-q'])+' '
        if not fastp_conf['-u'] == 0:
            com += "-u "+str(fastp_conf['-u'])+' '
        if not fastp_conf['-n'] == 0:
            com += "-n "+str(fastp_conf['-n'])+' '
        if not fastp_conf['-e'] == 0:
            com += "-e "+str(fastp_conf['-e'])+' '
        if not fastp_conf['-L'] == 0:
            com += "-L "
        if not fastp_conf['--length_required'] == 0:
            com += "--length_required "+str(fastp_conf['--length_required'])+' '
        if not fastp_conf['--length_limit'] == 0:
            com += "--length_limit "+str(fastp_conf['--length_limit'])+' '
        if not fastp_conf['--low_complexity_filter'] == 0:
            com += "--low_complexity_filter "
        if not fastp_conf['--complexity_threshold'] == 0:
            com += "--complexity_threshold "+str(fastp_conf['--complexity_threshold'])+' '
        if not fastp_conf['--filter_by_index1'] == 0:
            com += "--filter_by_index1 "+str(fastp_conf['--filter_by_index1'])+' '
        if not fastp_conf['--filter_by_index2'] == 0:
            com += "--filter_by_index2 "+str(fastp_conf['--filter_by_index2'])+' '
        if not fastp_conf['--filter_by_index_threshold'] == 0:
            com += "--filter_by_index_threshold "+str(fastp_conf['--filter_by_index_threshold'])+' '
        # if not fastp_conf['--correction'] == 0:
        #     com += "--correction "
        # if not fastp_conf['--overlap_len_require'] == 0:
        #     com += "--overlap_len_require "+str(fastp_conf['--overlap_len_require'])+' '
        # if not fastp_conf['--overlap_diff_limit'] == 0:
        #     com += "--overlap_diff_limit "+str(fastp_conf['--overlap_diff_limit'])+' '
        # if not fastp_conf['--overlap_diff_percent_limit'] == 0:
        #     com += "--overlap_diff_percent_limit "+str(fastp_conf['--overlap_diff_percent_limit'])+' '
        if not fastp_conf['--umi'] == 0:
            com += '--umi '
        if not fastp_conf['--umi_loc'] == 0:
            com += '--umi_loc '+fastp_conf['--umi_loc']+' '
        if not fastp_conf['--umi_len'] == 0:
            com += '--umi_len '+str(fastp_conf['--umi_len'])+' '
        if not fastp_conf['--umi_prefix'] == 0:
            com += '--umi_prefix '+fastp_conf['--umi_prefix']+' '
        if not fastp_conf['--umi_skip'] == 0:
            com += '--umi_skip '+str(fastp_conf['--umi_skip'])+' '
        if not fastp_conf['-p'] == 0:
            com += '-p '
        if not fastp_conf['-P'] == 0:
            com += '-P '+str(fastp_conf['-P'])+' '
        if not fastp_conf['-w'] == 0:
            com += '-w '+str(fastp_conf['-w'])
        subprocess.run(com, shell=True, check=True)
        print("end fastp_single_end")
    
    def fastp_paired_end(self):
        print("begin fastp_paired_end")
        fastp_conf = self.conf['preprocessing']['fastp']
        com = 'fastp --in1 '+self.input_file_1+" --in2 "+self.input_file_2+' --out1 '+self.result_dir+'/preprocessing/fastp_output_1.fastq --out2 '+self.result_dir+'/preprocessing/fastp_output_2.fastq -j '+self.result_dir+'/preprocessing/fastp.json -h '+self.result_dir+'/preprocessing/fastp.html '
        if not fastp_conf['-6'] == 0:
            com += '-6 '
        if not fastp_conf['--interleaved_in'] == 0:
            com += '--interleaved_in '
        if not fastp_conf['-V'] == 0:
            com += '-V '
        if not fastp_conf['-A'] == 0:
            com += '-A '
        if not fastp_conf['--adapter_sequence'] == 0:
            com += '--adapter_sequence '+fastp_conf['--adapter_sequence']+' '
        if not fastp_conf['--adapter_sequence_r2'] == 0:
            com += '--adapter_sequence_r2 '+fastp_conf['--adapter_sequence_r2']+' '
        if not fastp_conf['--adapter_fasta'] == 0:
            com += '--adapter_fasta '+fastp_conf['--adapter_fasta']+' '
        if not fastp_conf['--detect_adapter_for_pe'] == 0:
            com += '--detect_adapter_for_pe '
        if not fastp_conf['-f'] == 0:
            com += '-f '+str(fastp_conf['-f'])+' '
        if not fastp_conf['-t'] == 0:
            com += '-t '+str(fastp_conf['-t'])+' '
        if not fastp_conf['-b'] == 0:
            com += '-b '+str(fastp_conf['-b'])+' '
        if not fastp_conf['-F'] == 0:
            com += '-F '+str(fastp_conf['-F'])+' '
        if not fastp_conf['-T'] == 0:
            com += '-T '+str(fastp_conf['-T'])+' '
        if not fastp_conf['-B'] == 0:
            com += '-B '+str(fastp_conf['-B'])+' '
        if not fastp_conf['--trim_poly_g'] == 0:
            com += '--trim_poly_g '
        if not fastp_conf['--poly_g_min_len'] == 0:
            com += '--poly_g_min_len '+str(fastp_conf['--poly_g_min_len'])+' '
        if not fastp_conf['-G'] == 0:
            com += '-G '
        if not fastp_conf['--trim_poly_x'] == 0:
            com += '--trim_poly_x '
        if not fastp_conf['--poly_x_min_len'] == 0:
            com += '--poly_x_min_len '+str(fastp_conf['--poly_x_min_len'])+' '
        if not fastp_conf['-5'] == 0:
            com += '-5 '
        if not fastp_conf['-3'] == 0:
            com += '-3 '
        if not fastp_conf['-r'] == 0:
            com += '-r '
        if not fastp_conf['-W'] == 0:
            com += "-W "+str(fastp_conf['-W'])+' '
        if not fastp_conf['--cut_mean_quality'] == 0:
            com += "--cut_mean_quality "+str(fastp_conf['--cut_mean_quality'])+' '
        if not fastp_conf['--cut_front_window_size'] == 0:
            com += "--cut_front_window_size "+str(fastp_conf['--cut_front_window_size'])+' '
        if not fastp_conf['--cut_front_mean_quality'] == 0:
            com += "--cut_front_mean_quality "+str(fastp_conf['--cut_front_mean_quality'])+' '
        if not fastp_conf['--cut_tail_window_size'] == 0:
            com += "--cut_tail_window_size "+str(fastp_conf['--cut_tail_window_size'])+' '
        if not fastp_conf['--cut_tail_mean_quality'] == 0:
            com += "--cut_tail_mean_quality "+str(fastp_conf['--cut_tail_mean_quality'])+' '
        if not fastp_conf['--cut_right_window_size'] == 0:
            com += "--cut_right_window_size "+str(fastp_conf['--cut_right_window_size'])+' '
        if not fastp_conf['--cut_right_mean_quality'] == 0:
            com += "--cut_right_mean_quality "+str(fastp_conf['--cut_right_mean_quality'])+' '
        if not fastp_conf['-Q'] == 0:
            com += '-Q '
        if not fastp_conf['-q'] == 0:
            com += "-q "+str(fastp_conf['-q'])+' '
        if not fastp_conf['-u'] == 0:
            com += "-u "+str(fastp_conf['-u'])+' '
        if not fastp_conf['-n'] == 0:
            com += "-n "+str(fastp_conf['-n'])+' '
        if not fastp_conf['-e'] == 0:
            com += "-e "+str(fastp_conf['-e'])+' '
        if not fastp_conf['-L'] == 0:
            com += "-L "
        if not fastp_conf['--length_required'] == 0:
            com += "--length_required "+str(fastp_conf['--length_required'])+' '
        if not fastp_conf['--length_limit'] == 0:
            com += "--length_limit "+str(fastp_conf['--length_limit'])+' '
        if not fastp_conf['--low_complexity_filter'] == 0:
            com += "--low_complexity_filter "
        if not fastp_conf['--complexity_threshold'] == 0:
            com += "--complexity_threshold "+str(fastp_conf['--complexity_threshold'])+' '
        if not fastp_conf['--filter_by_index1'] == 0:
            com += "--filter_by_index1 "+str(fastp_conf['--filter_by_index1'])+' '
        if not fastp_conf['--filter_by_index2'] == 0:
            com += "--filter_by_index2 "+str(fastp_conf['--filter_by_index2'])+' '
        if not fastp_conf['--filter_by_index_threshold'] == 0:
            com += "--filter_by_index_threshold "+str(fastp_conf['--filter_by_index_threshold'])+' '
        if not fastp_conf['--correction'] == 0:
            com += "--correction "
        if not fastp_conf['--overlap_len_require'] == 0:
            com += "--overlap_len_require "+str(fastp_conf['--overlap_len_require'])+' '
        if not fastp_conf['--overlap_diff_limit'] == 0:
            com += "--overlap_diff_limit "+str(fastp_conf['--overlap_diff_limit'])+' '
        if not fastp_conf['--overlap_diff_percent_limit'] == 0:
            com += "--overlap_diff_percent_limit "+str(fastp_conf['--overlap_diff_percent_limit'])+' '
        if not fastp_conf['--umi'] == 0:
            com += '--umi '
        if not fastp_conf['--umi_loc'] == 0:
            com += '--umi_loc '+fastp_conf['--umi_loc']+' '
        if not fastp_conf['--umi_len'] == 0:
            com += '--umi_len '+str(fastp_conf['--umi_len'])+' '
        if not fastp_conf['--umi_prefix'] == 0:
            com += '--umi_prefix '+fastp_conf['--umi_prefix']+' '
        if not fastp_conf['--umi_skip'] == 0:
            com += '--umi_skip '+str(fastp_conf['--umi_skip'])+' '
        if not fastp_conf['-p'] == 0:
            com += '-p '
        if not fastp_conf['-P'] == 0:
            com += '-P '+str(fastp_conf['-P'])+' '
        if not fastp_conf['-w'] == 0:
            com += '-w '+str(fastp_conf['-w'])
        subprocess.run(com, shell=True, check=True)
        print("end fastp_paired_end")
    
    def trimmomatic_single_end(self):
        print("begin trimmomatic_single_end")
        trimmomatic_conf = self.conf['preprocessing']['trimmomatic']
        com = 'java -jar trimmomatic-0.38.jar SE '
        if not trimmomatic_conf['-threads'] == 0:
            com += "-threads "+str(trimmomatic_conf['-threads'])+' '
        if not trimmomatic_conf['-phred33'] == 0:
            com += "-phred33 "
        elif not trimmomatic_conf['-phred64'] == 0:
            com += "-phred64 "
        if not trimmomatic_conf['-validatePairs'] == 0:
            com += "-validatePairs "
        com += self.input_file+" "+self.result_dir+'/preprocessing/trimmomatic_output.fastq '
        if not trimmomatic_conf['ILLUMINACLIP'] == 0:
            com += "ILLUMINACLIP:"+trimmomatic_conf['ILLUMINACLIP']+' '
        if not trimmomatic_conf['LEADING'] == 0:
            com += "LEADING:"+str(trimmomatic_conf['LEADING'])+' '
        else:
            com += "LEADING:3 "
        if not trimmomatic_conf['TRAILING'] == 0:
            com += "TRAILING:"+str(trimmomatic_conf['TRAILING'])+' '
        else:
            com += "TRAILING:3 "
        if not trimmomatic_conf['SLIDINGWINDOW'] == 0:
            com += "SLIDINGWINDOW:"+trimmomatic_conf['SLIDINGWINDOW']+' '
        else:
            com += "SLIDINGWINDOW:4:15 "
        if not trimmomatic_conf['MINLEN'] == 0:
            com += "MINLEN:"+str(trimmomatic_conf['MINLEN'])
        else:
            com += "MINLEN:36"
        subprocess.run(com, shell=True, check=True)
        print("end trimmomatic_single_end")
        
    def trimmomatic_paired_end(self):
        print("begin trimmomatic_paired_end")
        trimmomatic_conf = self.conf['preprocessing']['trimmomatic']
        com = 'java -jar trimmomatic-0.38.jar PE '
        if not trimmomatic_conf['-threads'] == 0:
            com += "-threads "+str(trimmomatic_conf['-threads'])+' '
        if not trimmomatic_conf['-phred33'] == 0:
            com += "-phred33 "
        elif not trimmomatic_conf['-phred64'] == 0:
            com += "-phred64 "
        if not trimmomatic_conf['-validatePairs'] == 0:
            com += "-validatePairs "
        com += self.input_file_1+" "+self.input_file_2+" "+self.result_dir+'/preprocessing/paired_output_1.fastq '+self.result_dir+'/preprocessing/unpaired_output_1.fastq '+self.result_dir+'/preprocessing/paired_output_2.fastq '+self.result_dir+'/preprocessing/unpaired_output_2.fastq '
        if not trimmomatic_conf['ILLUMINACLIP'] == 0:
            com += "ILLUMINACLIP:"+trimmomatic_conf['ILLUMINACLIP']+' '
        if not trimmomatic_conf['LEADING'] == 0:
            com += "LEADING:"+str(trimmomatic_conf['LEADING'])+' '
        else:
            com += "LEADING:3 "
        if not trimmomatic_conf['TRAILING'] == 0:
            com += "TRAILING:"+str(trimmomatic_conf['TRAILING'])+' '
        else:
            com += "TRAILING:3 "
        if not trimmomatic_conf['SLIDINGWINDOW'] == 0:
            com += "SLIDINGWINDOW:"+trimmomatic_conf['SLIDINGWINDOW']+' '
        else:
            com += "SLIDINGWINDOW:4:15 "
        if not trimmomatic_conf['MINLEN'] == 0:
            com += "MINLEN:"+str(trimmomatic_conf['MINLEN'])
        else:
            com += "MINLEN:36"
        subprocess.run(com, shell=True, check=True)
        print("end trimmomatic_paired_end")

    def cutadapt_single_end(self):
        print("begin cutadapt_single_end")
        cutadapt_conf = self.conf['preprocessing']['cutadapt']
        com = 'cutadapt '
        if not cutadapt_conf['-a'] == 0:
            com += "-a "+cutadapt_conf['-a']+' '
        if not cutadapt_conf['-g'] == 0:
            com += "-g "+cutadapt_conf['-g']+' '
        if not cutadapt_conf['-b'] == 0:
            com += "-b "+cutadapt_conf['-b']+' '
        if not cutadapt_conf['-e'] == 0:
            com += "-e "+str(cutadapt_conf['-e'])+' '
        if not cutadapt_conf['--no-indels'] == 0:
            com += "--no-indels "
        if not cutadapt_conf['-n'] == 0:
            com += "-n "+str(cutadapt_conf['-n'])+' '
        if not cutadapt_conf['-O'] == 0:
            com += "-O "+str(cutadapt_conf['-O'])+' '
        if not cutadapt_conf['--match-read-wildcards'] == 0:
            com += "--match-read-wildcards "
        if not cutadapt_conf['-N'] == 0:
            com += "-N "
        if not cutadapt_conf['-u'] == 0:
            com += "-u "+str(cutadapt_conf['-u'])+' '
        if not cutadapt_conf['--nextseq-trim'] == 0:
            com += "--nextseq-trim "+str(cutadapt_conf['--nextseq-trim'])+' '
        if not cutadapt_conf['-q'] == 0:
            com += "-q "+str(cutadapt_conf['-q'])+' '
        if not cutadapt_conf['--quality-base'] == 0:
            com += "--quality-base "+str(cutadapt_conf['--quality-base'])+' '
        if not cutadapt_conf['--length'] == 0:
            com += "--length "+str(cutadapt_conf['--length'])+' '
        if not cutadapt_conf['--trim-n'] == 0:
            com += "--trim-n "
        if not cutadapt_conf['--length-tag'] == 0:
            com += "--length-tag "+cutadapt_conf['--length-tag']+' '
        if not cutadapt_conf['--strip-suffix'] == 0:
            com += "--strip-suffix "+cutadapt_conf['--strip-suffix']+' '
        if not cutadapt_conf['-x'] == 0:
            com += "-x "+cutadapt_conf['-x']+' '
        if not cutadapt_conf['-y'] == 0:
            com += "-y "+cutadapt_conf['-y']+' '
        if not cutadapt_conf['--zero-cap'] == 0:
            com += "--zero-cap "+cutadapt_conf['--zero-cap']+' '
        if not cutadapt_conf['-m'] == 0:
            com += "-m "+str(cutadapt_conf['-m'])+' '
        if not cutadapt_conf['-M'] == 0:
            com += "-M "+str(cutadapt_conf['-M'])+' '
        if not cutadapt_conf['--max-n'] == 0:
            com += "--max-n "+str(cutadapt_conf['--max-n'])+' '
        if not cutadapt_conf['--discard-trimmed'] == 0:
            com += "--discard-trimmed "
        if not cutadapt_conf['--discard-untrimmed'] == 0:
            com += "--discard-untrimmed "
        if not cutadapt_conf['--discard-casava'] == 0:
            com += "--discard-casava "
        com += '-o '+self.result_dir+'/preprocessing/cutadapt_output.fastq '+self.input_file
        subprocess.run(com, shell=True, check=True)
        print("end cutadapt_single_end")

    def cutadapt_paired_end(self):
        print("begin cutadapt_paired_end")
        cutadapt_conf = self.conf['preprocessing']['cutadapt']
        com = 'cutadapt '
        if not cutadapt_conf['-A'] == 0:
            com += "-A "+cutadapt_conf['-A']+' '
        if not cutadapt_conf['-G'] == 0:
            com += "-G "+cutadapt_conf['-G']+' '
        if not cutadapt_conf['-B'] == 0:
            com += "-B "+cutadapt_conf['-B']+' '
        if not cutadapt_conf['-U'] == 0:
            com += "-U "+str(cutadapt_conf['-U'])+' '
        if not cutadapt_conf['-e'] == 0:
            com += "-e "+str(cutadapt_conf['-e'])+' '
        if not cutadapt_conf['--no-indels'] == 0:
            com += "--no-indels "
        if not cutadapt_conf['-n'] == 0:
            com += "-n "+str(cutadapt_conf['-n'])+' '
        if not cutadapt_conf['-O'] == 0:
            com += "-O "+str(cutadapt_conf['-O'])+' '
        if not cutadapt_conf['--match-read-wildcards'] == 0:
            com += "--match-read-wildcards "
        if not cutadapt_conf['-N'] == 0:
            com += "-N "
        if not cutadapt_conf['-u'] == 0:
            com += "-u "+str(cutadapt_conf['-u'])+' '
        if not cutadapt_conf['--nextseq-trim'] == 0:
            com += "--nextseq-trim "+str(cutadapt_conf['--nextseq-trim'])+' '
        if not cutadapt_conf['-q'] == 0:
            com += "-q "+str(cutadapt_conf['-q'])+' '
        if not cutadapt_conf['--quality-base'] == 0:
            com += "--quality-base "+str(cutadapt_conf['--quality-base'])+' '
        if not cutadapt_conf['--length'] == 0:
            com += "--length "+str(cutadapt_conf['--length'])+' '
        if not cutadapt_conf['--trim-n'] == 0:
            com += "--trim-n "
        if not cutadapt_conf['--length-tag'] == 0:
            com += "--length-tag "+cutadapt_conf['--length-tag']+' '
        if not cutadapt_conf['--strip-suffix'] == 0:
            com += "--strip-suffix "+cutadapt_conf['--strip-suffix']+' '
        if not cutadapt_conf['-x'] == 0:
            com += "-x "+cutadapt_conf['-x']+' '
        if not cutadapt_conf['-y'] == 0:
            com += "-y "+cutadapt_conf['-y']+' '
        if not cutadapt_conf['--zero-cap'] == 0:
            com += "--zero-cap "+cutadapt_conf['--zero-cap']+' '
        if not cutadapt_conf['-m'] == 0:
            com += "-m "+str(cutadapt_conf['-m'])+' '
        if not cutadapt_conf['-M'] == 0:
            com += "-M "+str(cutadapt_conf['-M'])+' '
        if not cutadapt_conf['--max-n'] == 0:
            com += "--max-n "+str(cutadapt_conf['--max-n'])+' '
        if not cutadapt_conf['--discard-trimmed'] == 0:
            com += "--discard-trimmed "
        if not cutadapt_conf['--discard-untrimmed'] == 0:
            com += "--discard-untrimmed "
        if not cutadapt_conf['--discard-casava'] == 0:
            com += "--discard-casava "
        if not cutadapt_conf['--pair-adapters'] == 0:
            com += "--pair-adapters "
        if not cutadapt_conf['--pair-filter'] == 0:
            com += "--pair-filter "+cutadapt_conf['--pair-filter']+' '
        if not cutadapt_conf['--interleaved'] == 0:
            com += "--interleaved "
        com += '-o '+self.result_dir+'/preprocessing/cutadapt_output_1.fastq -p '+self.result_dir+'/preprocessing/cutadapt_output_2.fastq '+self.input_file_1+' '+self.input_file_2
        subprocess.run(com, shell=True, check=True)
        print("end cutadapt_paired_end")

def preprocessing_single_end(input_file, result_dir):
    print("Begin preprocessing")
    print("input_file="+input_file)
    subprocess.run('fastqc -o '+result_dir+'/preprocessing/ '+input_file, shell=True, check=True)
    subprocess.run('java -jar trimmomatic-0.38.jar SE '+input_file+" "+result_dir+'/preprocessing/output.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:2 MINLEN:25', shell=True, check=True)
    print("End preprocessing")

def preprocessing_paired_end(input_file_1, input_file_2, result_dir):
    print("Begin preprocessing")
    print("input_file="+input_file_1+input_file_2)
    subprocess.run('fastqc -o '+result_dir+'/preprocessing/ '+input_file_1+" "+input_file_2, shell=True, check=True)
    subprocess.run('java -jar trimmomatic-0.38.jar PE '+input_file_1+" "+input_file_2+" "+result_dir+'/preprocessing/paired_output_1.fq.gz '+result_dir+'/preprocessing/unpaired_output_1.fq.gz '+result_dir+'/preprocessing/paired_output_2.fq.gz '+result_dir+'/preprocessing/unpaired_output_2.fq.gz '+' ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36', shell=True, check=True)
    print("End preprocessing")
