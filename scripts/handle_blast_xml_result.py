import os
import sys
from Bio import SearchIO
import scripts.hit
import scripts.hsp
import traceback


def handle_blast_xml_result(filename):
    print('begin handle_blast_xml_result')
    accession_version_dict = {}
    # if os.path.isfile(filename):
    #     print("Handling...")
    #     blast_qresult = SearchIO.parse(filename, "blast-xml")
    #     for hit in blast_qresult:
    #         print(str(hit.id))
    #         accession_version = str(hit.id).split('|')[1]
    #         if accession_version in accession_version_dict:
    #             accession_version_dict[accession_version] += 1
    #         else:
    #             accession_version_dict[accession_version] = 1
    #     print("Handled")
    blast_result = open(filename, 'r')
    line = blast_result.readline()
    while line:
        if '<Hit_id>' in line and '</Hit_id>' in line:
            accession_version = str(line).split('|')[1]
            if accession_version in accession_version_dict:
                accession_version_dict[accession_version] += 1
            else:
                accession_version_dict[accession_version] = 1
        line = blast_result.readline()
    blast_result.close()
    print('end handle_blast_xml_result')
    return sorted(accession_version_dict.items(), key=lambda d: d[1], reverse=True)


def handle_blast_xml_result_outfmt7(filename):
    blast_result_file = open(filename, 'r')
    blast_lines = blast_result_file.readlines()
    blast_result_file.close()

    accession_version_dict = {}

    for blast_line in blast_lines:
        if not blast_line.startswith('#'):
            # items[1] = subject acc.ver; items[2] = % identity; items[3] = alignment length
            items = blast_line.split()
            items[2] = float(items[2])
            items[3] = float(items[3])
            if not accession_version_dict.__contains__(items[1]):
                accession_version_dict[items[1]] = items[2]*0.01*items[3]
            else:
                accession_version_dict[items[1]] += items[2]*0.01*items[3]

    accession_version_list = []
    for key in accession_version_dict:
        # map_identification_list.append([key, family_weight_dict[key][0] / family_weight_dict[key][1]])
        accession_version_list.append([key, accession_version_dict[key]])

    accession_version_list = sorted(accession_version_list, key=lambda x:x[1], reverse=True)
    print(accession_version_list[:50])
    return accession_version_list


def handle_blast_xml_result_outfmt7_v2(blast_lines):

    accession_version_dict = {}

    for blast_line in blast_lines:
        if not blast_line[0].startswith('#'):
            blast_line[2] = float(blast_line[2])
            blast_line[3] = float(blast_line[3])
            if not accession_version_dict.__contains__(blast_line[1]):
                accession_version_dict[blast_line[1]] = blast_line[2]*0.01*blast_line[3]
            else:
                accession_version_dict[blast_line[1]] += blast_line[2]*0.01*blast_line[3]

    accession_version_list = []
    for key in accession_version_dict:
        accession_version_list.append([key, accession_version_dict[key]])

    accession_version_list = sorted(accession_version_list, key=lambda x:x[1], reverse=True)
    print(accession_version_list[:50])
    return accession_version_list


def handle(result_dir, res, subdir):
    try:
        print('begin handle blast result')
        out = open(result_dir+'/resequencing/'+subdir, 'w')
        
        hit = scripts.hit.Hit()
        hsp = scripts.hsp.Hsp()

        for l in res:

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

    except Exception as e:
        # print(e)
        traceback.print_exc()
    finally:
        out.close()
        print('end handle blast result')

# def handle_blast_xml_result(filename):
#     if not os.path.isfile(filename):
#         gi_dict = {}
#         print("Handling...")
#         blast_qresult = SearchIO.read(filename, "blast-xml")
#         for hit in blast_qresult:
#             gi = str(hit.id).split('|')[1]
#             if gi in gi_dict:
#                 gi_dict[gi] += 1
#             else:
#                 gi_dict[gi] = 1
#         print("Handled")
#     return sorted(gi_dict.items(), key=lambda d: d[1], reverse=True)

# def handle_blast_xml_result():
#     argv = sys.argv
#     filename = argv[1]
#     if not os.path.isfile(filename):
#         gi_dict = {}
#         print("Handling...")
#         blast_qresult = SearchIO.read(filename, "blast-xml")
#         for hit in blast_qresult:
#             gi = str(hit.id).split('|')[1]
#             if gi in gi_dict:
#                 gi_dict[gi] += 1
#             else:
#                 gi_dict[gi] = 1
#         print("Handled")
#     return sorted(gi_dict.items(), key=lambda d: d[1], reversed=True) # [("a", 4), ("b", 2), ("c", 1)]

# if __name__ == "__main__":
#     handle_blast_xml_result()
