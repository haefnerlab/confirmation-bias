import os.path as path
import pickle as pk

def open_file(name):
    if path.exists(name):
        with open(name, 'rb') as handle:
            file = pk.load(handle)
    else:
        return -1
    return file
def save_file(name, file):
    with open(name, 'wb') as handle:
        pk.dump(file, handle, protocol=pk.HIGHEST_PROTOCOL)