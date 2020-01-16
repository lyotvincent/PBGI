"""
author: sun jialiang
environment: python3.6
how to run: 
"""

import os
import sys
import time
import json

import next_generation_sequencing
import third_generation_sequencing
import preprocessing
import assembly

#
def main():
    time_note2 = time.time()
    time_note = time.process_time()
    argv = sys.argv
    print(str(argv))
    # for i in argv:
    #     print(str(argv[argv.index(i)]))

    with open('./conf.json', 'r') as load_f:
        conf = json.load(load_f)

    if len(argv) == 1 or "-h" in argv or "-help" in argv or "--h" in argv or "--help" in argv:
        print_help()
        exit()

    # input file
    if "-f" in argv and not "-1" in argv and not "-2" in argv and not "-12" in argv:
        input_file = argv[argv.index("-f") + 1]
        mode = 1
    elif "-12" in argv and not "-f" in argv and not "-1" in argv and not "-2" in argv:
        input_file_12 = argv[argv.index("-12") + 1]
        mode = 12
    elif "-1" in argv and "-2" in argv and not "-f" in argv and not "-12" in argv:
        input_file_1 = argv[argv.index("-1") + 1]
        input_file_2 = argv[argv.index("-2") + 1]
        mode = 2
    else:
        print('Please input -f or -1 and -2. And dont input -f, -1 and -2 at the same time.')
        exit()

    ngs = None
    thirdgs = None
    if '-ngs' in argv and not '-3gs' in argv:
        ngs = 1
    elif '-3gs' in argv and not '-ngs' in argv:
        thirdgs = 1
    else:
        print('please choose ngs or 3gs')

    if '-o' in argv:
        result_dir = argv[argv.index("-o") + 1] + '_'
    else:
        result_dir = ''
    
    start_time = time.strftime("%Y_%m_%d_%H:%M:%S", time.localtime())
    result_dir += 'result_' + start_time
    print("result_dir=" + os.path.abspath('.') + '/' + result_dir)

    os.mkdir(os.path.abspath('.') + '/' + result_dir)

    if ngs == 1:
        if mode == 1:
            ngs = next_generation_sequencing.NextGenerationSequencing(result_dir, conf, input_file=input_file)
            ngs.run_single()
        elif mode == 12:
            ngs = next_generation_sequencing.NextGenerationSequencing(result_dir, conf, input_file_12=input_file_12)
            ngs.run_interleaved()
        elif mode == 2:
            ngs = next_generation_sequencing.NextGenerationSequencing(result_dir, conf, input_file_1=input_file_1, input_file_2=input_file_2)
            ngs.run_paired()
    if thirdgs == 1:
        print('3gs')
        thirdgs = third_generation_sequencing.ThirdGenerationSequencing(result_dir, conf, input_file)
        thirdgs.run()

    print('starttime='+start_time)
    print('endtime='+time.strftime("%Y_%m_%d_%H:%M:%S", time.localtime()))
    print('process_time='+str(time.process_time()-time_note))
    print('time='+str(time.time()-time_note2))
    print("End")


def print_help():
    print("help:")
    print("-f input single end file")
    print("-1 and -2 input paired end files")
    print("-o out directory name")
    print("-h or -help or --h or --help look help")


if __name__ == "__main__":
    main()
