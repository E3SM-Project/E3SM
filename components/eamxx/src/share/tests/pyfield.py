import numpy as np

def set_iota(arr):
    arr[:] = np.arange(np.prod(arr.shape)).reshape(arr.shape)
