import subprocess
import numpy as np

def read_pt2(filename, channel):
    d = subprocess.check_output(['extract_pt2_timestamps', '-c', str(channel)], stdin=open(filename))
    return np.frombuffer(d, dtype='u8')
