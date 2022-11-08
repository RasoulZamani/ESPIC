"""
This file (utils.py) contains some tools for using in other
pice of code(e.g for plotting and assistance functions)
"""

# plotter ____________________________________________________
import time
import sys
import numpy as np
from parameters import *
from matplotlib import pyplot as plt

# plotter _____________________________________________________
#....................
#   ...... TODO .....
# ...................






# timming ____________________________________________________
def timing (f):
    """
    USE @timing before function def to show exec time 
    when function is called
    """
    def wrap (*args):
        time1 = time.time()
        ret = f(*args)
        time2 = time.time()
        print('%s function took %0.3f ms' % (f.func_name,
                               (time2 - time1) * 1000.0) )
        return ret
    return wrap


# Print iterations progress___________________________________

# replase by tqdm 

# TODO

 
def printProgress (iteration, total, prefix = '', suffix = '', decimals = 2, barLength = 100):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : number of decimals in percent complete (Int)
        barLength   - Optional  : character length of bar (Int)
    """
    filledLength    = int(round(barLength * iteration / float(total)))
    percents        = round(100.00 * (iteration / float(total)), decimals)
    bar             = u'\u2588' * filledLength + '-' * (barLength - filledLength)
    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),
    sys.stdout.flush()
    if iteration == total:
        sys.stdout.write('\n')
        sys.stdout.flush()
        
