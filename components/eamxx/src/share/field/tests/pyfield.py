import numpy as np

# Simple function using numpy to set all entries of an array to 0,1,2,...
def set_iota(arr):
    arr[:] = np.arange(np.prod(arr.shape)).reshape(arr.shape)
