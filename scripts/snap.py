import subprocess

def snap(result_dir, accession_version_list, input_file_1, input_file_2=None):
    closest_accession_version = accession_version_list[0]
    max_similarity = 0
    if input_file_2 == None:
        for accession_version in accession_version_list:
            subprocess.run('snap-aligner index '+result_dir+'/resequencing/'+accession_version+'.fasta '+result_dir+'/resequencing/'+accession_version+'_index', shell=True, check=True)
            subprocess.run('snap-aligner single '+result_dir+'/resequencing/'+accession_version+'_index '+input_file_1+' -o '+result_dir+'/resequencing/identification_'+accession_version+'.bam', shell=True, check=True)
            subprocess.run('samtools sort '+result_dir+'/resequencing/identification_'+accession_version+'.bam -o '+result_dir+'/resequencing/identification_'+accession_version+'.sorted.bam -O bam', shell=True, check=True)
            res = subprocess.Popen('bedtools genomecov -ibam '+result_dir+'/resequencing/identification_'+accession_version+'.sorted.bam', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
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
            # print('send_message='+str(send_message))
            # print('send_message[-1]='+str(send_message[-1]))
            # if send_message[-1].decode('utf-8') == 'Aligning.\n':
            #     break
            # new_aligned = float(str(send_message[-1]).split()[2][1:-2])
            # if new_aligned > aligned:
            #     aligned = new_aligned
            #     closest_accession_version = accession_version
            # subprocess.run('rm -rf '+accession_version+'_index', shell=True, check=True)
    else:
        for accession_version in accession_version_list:
            subprocess.run('snap-aligner index '+result_dir+'/resequencing/'+accession_version+'.fasta '+result_dir+'/resequencing/'+accession_version+'_index', shell=True, check=True)
            subprocess.run('snap-aligner paired '+result_dir+'/resequencing/'+accession_version+'_index '+input_file_1+' '+input_file_2+' -o '+result_dir+'/resequencing/identification_'+accession_version+'.bam', shell=True, check=True)
            subprocess.run('samtools sort '+result_dir+'/resequencing/identification_'+accession_version+'.bam -o '+result_dir+'/resequencing/identification_'+accession_version+'.sorted.bam -O bam', shell=True, check=True)
            res = subprocess.Popen('bedtools genomecov -ibam '+result_dir+'/resequencing/identification_'+accession_version+'.sorted.bam', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
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
            # res = subprocess.Popen('snap-aligner paired '+result_dir+'/resequencing/'+accession_version+'_index '+input_file_1+' '+input_file_2+' -o '+result_dir+'/resequencing/'+accession_version+'.sam', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
            # send_message = res.stdout.readlines()
            # print('send_message='+str(send_message))
            # if send_message[-1].decode('utf-8') == 'Aligning.\n':
            #     break
            # new_aligned = float(str(send_message[-1]).split()[2][1:-2])
            # if new_aligned > aligned:
            #     aligned = new_aligned
            #     closest_accession_version = accession_version
            # subprocess.run('rm -rf '+accession_version+'_index', shell=True, check=True)
    return closest_accession_version

            
            