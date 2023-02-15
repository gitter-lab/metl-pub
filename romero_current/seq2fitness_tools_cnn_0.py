# Import necessary libraries here
from os import listdir
from os.path import isfile, join, exists
import numpy as np
from simple_inference import encoding, inference

class seq2fitness_handler:
    def __init__(self):
        pid = '8852736'
        model_num = 0

        self.model_path = 'models_{}/model_{}'.format(pid, model_num)
        if (exists(self.model_path)):
            for file_name in listdir(self.model_path):
                if '.pb' in file_name:
                    model_name = file_name
                    print(file_name)
                    self.model_path = self.model_path+'/'+model_name

        self.model = inference.restore_sess_from_pb(self.model_path)


    def seq2fitness(self, seq):
        encoded_seq = encoding.encode(encoding="one_hot,aa_index", char_seqs=[seq])
        return inference.run_inference(encoded_data=encoded_seq, sess=self.model)