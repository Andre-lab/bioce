"""
1. Reads in VBI output and file_list
2. Selects non-zero models and outputs short_lists intenstites and
"""
import os
import sys
import numpy as np

def read_file_safe(filename, dtype="float64"):
    """
    Simple check if file exists
    :param filename:
    :return:
    """
    try:
        results = np.genfromtxt(filename, dtype=dtype)
    except IOError as err:
        print(os.strerror(err.errno))
    return results

def savetext(filename, string_array):
    """

    :param string_array:
    :return:
    """
    output_file = open(filename,'w')
    for name in string_array:
        output_file.write(name+" ")
    output_file.close()

if __name__=="__main__":

    intensities = read_file_safe(sys.argv[2])
    vbi_output = read_file_safe(sys.argv[1])
    last_output = vbi_output[-1][:-4]
    file_list = read_file_safe(sys.argv[3], 'unicode')

    nonzero_indexes = np.nonzero(last_output>0.0)
    print(nonzero_indexes)

    output_intensities = intensities[:,nonzero_indexes[0]]
    output_filelist = file_list[last_output>0.0]

    savetext('file_sub.txt', output_filelist)
    np.savetxt('TrmIntensites.txt', output_intensities)