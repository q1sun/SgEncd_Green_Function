import torch
import numpy as np

# exact Green's function
def Green_Exact_Line1D(x, y):

	return x * (1-y) * (x<=y) + y * (1-x) * (x>y)


def augmented_variable(x, y):
     
    return torch.abs(x - y)






