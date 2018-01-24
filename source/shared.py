import os


def format_dir(directory):
    if not os.path.isdir(directory):
        print(directory, " not found!")
        exit(1)
    if directory[-1] != '/':
        return directory + '/'
    return directory


def print_dict_keys(dictionary):
    for key in dictionary.keys():
        print(key + " ")
