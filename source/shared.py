import os


def format_dir(directory):
    if not os.path.isdir(directory):
        print(directory, " directory not found. Attempting to create . . .")
        try:
            os.makedirs(directory, exist_ok=False)
        except OSError:
            print("Could not create ", directory, ".")
            exit(1)
        print("Done.")
    if directory[-1] != '/':
        return directory + '/'
    return directory


def print_dict_keys(dictionary):
    for key in dictionary.keys():
        print(key + " ")

if __name__ == '__main__':
    print ("Loaded shared.py module.")
