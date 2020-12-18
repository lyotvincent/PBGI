import scripts.hit
import scripts.hsp
import traceback

def handle(result_dir, file_name, subdir):
    try:
        print('begin handle blast result')
        reader = open(file_name, 'r')
        out = open(result_dir+'/denovo/'+subdir, 'w')
        
        hit = scripts.hit.Hit()
        hsp = scripts.hsp.Hsp()

        l = reader.readline()
        while l:

            if not l.find('<BlastOutput_version>') == -1:
                tmp = l[l.find('<BlastOutput_version>')+21: l.find('</BlastOutput_version>')]
                out.write('#'+tmp+'\n')
            if not l.find('<BlastOutput_db>') == -1:
                tmp = l[l.find('<BlastOutput_db>')+16: l.find('</BlastOutput_db>')]
                out.write('#'+tmp+'\n')
            if not l.find('<Parameters_expect>') == -1:
                tmp = l[l.find('<Parameters_expect>')+19: l.find('</Parameters_expect>')]
                out.write('#'+tmp+'\n')
            
            if not l.find('<Iteration_query-def>') == -1:
                tmp = l[l.find('<Iteration_query-def>')+21: l.find('</Iteration_query-def>')]
                out.write('@'+tmp+'\n')
            
            
            if not l.find('<Hit_id>') == -1:
                tmp = l[l.find('<Hit_id>')+8: l.find('</Hit_id>')]
                hit.id = tmp
            
            if not l.find('<Hit_def>') == -1:
                tmp = l[l.find('<Hit_def>')+9: l.find('</Hit_def>')]
                hit.desc = tmp
            
            if not l.find('<Hit_len>') == -1:
                tmp = l[l.find('<Hit_len>')+9: l.find('</Hit_len>')]
                hit.length = int(tmp)
            # bit score
            if not l.find('<Hsp_bit-score>') == -1:
                tmp = l[l.find('<Hsp_bit-score>')+15: l.find('</Hsp_bit-score>')]
                hsp.score = tmp
            
            if not l.find('<Hsp_evalue>') == -1:
                tmp = l[l.find('<Hsp_evalue>')+12: l.find('</Hsp_evalue>')]
                hsp.evalue = tmp
            
            if not l.find('<Hsp_query-from>') == -1:
                tmp = l[l.find('<Hsp_query-from>')+16: l.find('</Hsp_query-from>')]
                hsp.q_start = int(tmp)
            
            if not l.find('<Hsp_query-to>') == -1:
                tmp = l[l.find('<Hsp_query-to>')+14: l.find('</Hsp_query-to>')]
                hsp.q_end = int(tmp)
            
            if not l.find('<Hsp_hit-from>') == -1:
                tmp = l[l.find('<Hsp_hit-from>')+14: l.find('</Hsp_hit-from>')]
                hsp.s_start = int(tmp)
            
            if not l.find('<Hsp_hit-to>') == -1:
                tmp = l[l.find('<Hsp_hit-to>')+12: l.find('</Hsp_hit-to>')]
                hsp.s_end = int(tmp)
            
            # end of hit
            if not l.find('</Hit>') == -1:
                out.write(hit.getHit()+'\n')
                hit = scripts.hit.Hit()

            # end of hsp
            if not l.find('</Hsp>') == -1:
                hit.hsps.append(hsp)
                hsp = scripts.hsp.Hsp()

            l = reader.readline()
    except Exception as e:
        # print(e)
        traceback.print_exc()
    finally:
        reader.close()
        out.close()
        print('end handle blast result')
