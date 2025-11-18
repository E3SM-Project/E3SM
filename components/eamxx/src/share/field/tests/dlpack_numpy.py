import numpy as np


# Simple function using numpy to set all entries of an array to 0,1,2,...
def set_iota(dltens):
    print (f"np version: {np.__version__}")
    print (f"dltens type: {type(dltens)}")
    print (f"dltens dir: {dltens.__dir__()}")
    print (f"dltens repr: {dltens.__repr__()}")
    dlt = dltens.__dlpack__()
    print (f"dltens: {dlt.__dir__()}")
    arr = np.from_dlpack(dltens)
    print(f"flags: {arr.flags}")
    print(f"writable: {arr.flags['WRITEABLE']}")
    print(f"methods: {arr.__dir__()}")
    arr.setflags(write=True)
    #  arr.flags['WRITEABLE'] = True
    print (arr)
    arr[:] = np.arange(np.prod(arr.shape)).reshape(arr.shape)
