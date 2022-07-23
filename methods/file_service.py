import re
import nums_from_string
import json
import numpy as np

def write_to_file(str, url):
    file = open(url, 'w')
    file.write(str)
    file.close()

def read_from_file(url):
    file = open(url, 'r')
    text = file.read()
    file.close()
    return text

def write_table_to_file(table:str, url):
    file_text = nums_from_string.get_nums(table)
    write_to_file(str(file_text), url)

def read_table_to_array(url):
    text = read_from_file(url)
    float_array = nums_from_string.get_nums(text)
    float_array = [float(i) for i in float_array]
    return float_array


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

