import os
import torch
import shutil

import sys
import numpy as np

'''Save best model to checkpoint'''
def save_checkpoint(state, is_best, checkpoint='checkpoint', filename='checkpoint.pth.tar'):
    filepath = os.path.join(checkpoint, filename)
    torch.save(state, filepath)
    if is_best:
        shutil.copyfile(filepath, os.path.join(checkpoint, 'model_best.pth.tar'))
        
'''Save Dataset for repeat experiment'''
def save_dataset(dataset,  checkpoint='checkpoint', filename='dataset.pth.tar'):
    filepath = os.path.join(checkpoint, filepath)
    torch.save(dataset)
    
    

'''Save learning curves during training and testing'''
def save_learncurve(state, curve='curve', filename='curve.pt'):
    filepath = os.path.join(curve, filename)
    torch.save(state, filepath)  
         
class Logger(object):
    '''Save training process to log file'''
    def __init__(self, fpath, title=None): 
        self.file = None
        self.title = '' if title == None else title
        if fpath is not None:
            self.file = open(fpath, 'w')

    def set_names(self, names):
        # initialize numbers as empty list
        self.numbers = {}
        self.names = names
        for _, name in enumerate(self.names):
            self.file.write(name)
            self.file.write('\t')
            self.numbers[name] = []
        self.file.write('\n')
        self.file.flush()

    def append(self, numbers):
        assert len(self.names) == len(numbers), 'Numbers do not match names'
        for index, num in enumerate(numbers):
            self.file.write("{0:.8f}".format(num))
            self.file.write('\t')
            self.numbers[self.names[index]].append(num)
        self.file.write('\n')
        self.file.flush()

    def close(self):
        if self.file is not None:
            self.file.close()  

class LossTracker:
    def __init__(self):
        self.loss_intrr = 0
        self.loss_bndry = 0
        self.loss_intfc = 0
        self.loss_symtry = 0

    def update(self, loss_intrr, loss_bndry, loss_intfc, loss_symtry, intrr_size, bndry_size, intfc_size):
        self.loss_intrr = loss_intrr.item()
        self.loss_bndry = loss_bndry.item()
        self.loss_intfc = loss_intfc.item()
        self.loss_symtry = loss_symtry.item()

        self.intrr_size = intrr_size
        self.bndry_size = bndry_size
        self.intfc_size = intfc_size

    def get_losses(self):
        return self.loss_intrr, self.loss_bndry, self.loss_intfc, self.loss_symtry, self.intrr_size, self.bndry_size, self.intfc_size


def mkdir_p(path):
    '''make dir if not exist'''
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise              
