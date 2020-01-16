import preprocessing
import resequencing
import denovo
import os

class ThirdGenerationSequencing:

    input_file = None
    result_dir = None
    conf = None

    def __init__(self, result_dir, conf, input_file):
        self.result_dir = result_dir
        self.conf = conf
        self.input_file = input_file

    def run(self):
        if not self.conf['preprocessing']['enable'] == 0:
            os.mkdir(os.path.abspath('.') + '/' + self.result_dir + '/' + 'preprocessing')
            preprocessing_obj = preprocessing.Preprocessing(self.result_dir, self.conf, input_file=self.input_file)
            if not self.conf['preprocessing']['fastqc']['enable'] == 0:
                preprocessing_obj.fastqc_single_end()

        if not self.conf['resequencing']['enable'] == 0:
            os.mkdir(os.path.abspath('.') + '/' + self.result_dir + '/' + 'resequencing')
            resequencing_obj = resequencing.Resequencing(self.result_dir, self.conf, self.input_file)
            resequencing_obj.run_3gs()
        if not self.conf['denovo']['enable'] == 0:
            os.mkdir(os.path.abspath('.') + '/' + self.result_dir + '/' + 'denovo')
            denovo_obj = denovo.Denovo(self.result_dir, self.conf, input_file=self.input_file)
            denovo_obj.run_3gs()
