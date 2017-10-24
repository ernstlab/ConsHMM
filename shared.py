from __future__ import print_function
import os

def formatDir(directory):
    if not os.path.isdir(directory):
        print(directory, " not found!")
        exit(1)
    if directory[-1] != '/':
        return directory + '/'
    return directory
