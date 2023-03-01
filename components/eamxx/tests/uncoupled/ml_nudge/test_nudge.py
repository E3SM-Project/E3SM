import numpy as np


def modify_view(data, Ncol, Nlev):
    data = np.reshape(data, (-1, Nlev))
    data[1, 10] += 0.1
