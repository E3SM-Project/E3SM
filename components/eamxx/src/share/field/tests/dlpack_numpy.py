import numpy as np
import ctypes


# Simple function using numpy to set all entries of an array to 0,1,2,...
def set_iota(dltens):
    print (f"np version: {np.__version__}")
    print (f"dltens type: {type(dltens)}")
    print (f"dltens dir: {dltens.__dir__()}")
    print (f"dltens repr: {dltens.__repr__()}")
    print (f"dltens subclass hook: {dltens.__subclasshook__()}")
    print (f"dltens class: {dltens.__class__()}")
    arr = np.from_dlpack(dltens)
    arr[:] = np.arange(np.prod(arr.shape)).reshape(arr.shape)
