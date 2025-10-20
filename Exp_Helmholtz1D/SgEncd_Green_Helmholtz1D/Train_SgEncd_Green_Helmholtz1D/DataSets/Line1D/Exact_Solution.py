import torch
import numpy as np

def c_Exact_Line1D(x):

    return (x-2) **2

def cx_Exact_Line1D(x):

    return 2 * (x-2)

def k_Exact_Line1D(x):

    return 15 * torch.sin(10 * x)

def augmented_variable(x, y):

    return torch.abs( x - y )
