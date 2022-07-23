import numpy as np

from .models import Matrix
from .file_service import read_table_to_array, write_to_file

def get_matrix_parametrs():
    d = dict();
    d['stages'] = Matrix._meta.get_field('stages').value_from_object(Matrix.objects.first())
    d['components'] = Matrix._meta.get_field('components').value_from_object(Matrix.objects.first())
    d['experients'] = Matrix._meta.get_field('experients').value_from_object(Matrix.objects.first())
    d['step'] = Matrix._meta.get_field('step').value_from_object(Matrix.objects.first())
    return d

def array_to_matrix_dictionary():
    array = read_table_to_array('methods/static/methods/json/tables.json')
    matrixParam = get_matrix_parametrs()

    d = dict();
    d['stehCoef'] = array[:matrixParam['stages'] * matrixParam['components']]
    d['pokazStep'] = array[matrixParam['stages'] * matrixParam['components']:matrixParam['stages'] * matrixParam['components'] * 2]
    d['expDat'] = array[matrixParam['stages'] * matrixParam['components'] * 2: matrixParam['stages'] * matrixParam['components'] * 2 + matrixParam['experients'] * (matrixParam['components']+1)]
    d['constSpeed'] = array[matrixParam['stages'] * matrixParam['components'] * 2 + matrixParam['experients'] * (matrixParam['components']+1):]
    return d

def system_diff_equation(y):
    dictionary = array_to_matrix_dictionary()
    matrixParam = get_matrix_parametrs()

    sum = np.zeros(matrixParam['components']);
    A = np.zeros(matrixParam['components'])

    for i in range(0,matrixParam['components']):
        for j in range(0,matrixParam['stages']):

            A[i] = dictionary['stehCoef'][i+(j*matrixParam['components'])] * dictionary['constSpeed'][j] * y[j] ** dictionary['pokazStep'][i+(j*matrixParam['components'])]
            sum[i] += A[i]

    return sum

def method_Euler():
    dictionary = array_to_matrix_dictionary()
    matrixParam = get_matrix_parametrs()

    t = dictionary['expDat'][1]
    h = matrixParam['step']
    n = int(dictionary['expDat'][(matrixParam['experients'] * (matrixParam['components'] + 1) - (matrixParam['components'] + 1) )])
    steps = int(n/h)

    y = np.zeros((steps+1,matrixParam['components']+1))
    x = 0

    for i in range (steps+1):
        y[i, 0] = x
        x += h

    for j in range(matrixParam['components']):
        y[0,j+1] = dictionary['expDat'][j+1]

    for i in range(0,steps,1):
        for j in range(matrixParam['components']):
            y[i+1,j+1] =  y[i,j+1] + h * system_diff_equation(y[i,1:])[j]

    return y

def method_RungeKuttaSecond():
    dictionary = array_to_matrix_dictionary()
    matrixParam = get_matrix_parametrs()

    t = dictionary['expDat'][1]
    h = matrixParam['step']
    n = int(dictionary['expDat'][
                (matrixParam['experients'] * (matrixParam['components'] + 1) - (matrixParam['components'] + 1))])
    steps = int(n / h)

    y = np.zeros((steps + 1, matrixParam['components'] + 1))
    x = 0
    for i in range(steps + 1):
        y[i, 0] = x
        x += h

    for j in range(matrixParam['components']):
        y[0, j + 1] = dictionary['expDat'][j + 1]

    for i in range(0,steps,1):
        for j in range(matrixParam['components']):
            y[i+1,j+1] =  y[i,j+1] + (h / 2) * (system_diff_equation(y[i,1:])[j] + system_diff_equation(y[i,1:] + h * system_diff_equation(y[i,1:])[j])[j])

    return y

def method_RungeKuttaFourth():
    dictionary = array_to_matrix_dictionary()
    matrixParam = get_matrix_parametrs()

    t = dictionary['expDat'][1]
    h = matrixParam['step']
    n = int(dictionary['expDat'][
                (matrixParam['experients'] * (matrixParam['components'] + 1) - (matrixParam['components'] + 1))])
    steps = int(n / h)

    y = np.zeros((steps + 1, matrixParam['components'] + 1))
    x = 0
    for i in range(steps + 1):
        y[i, 0] = x
        x += h

    for j in range(matrixParam['components']):
        y[0, j + 1] = dictionary['expDat'][j + 1]

    for i in range(0,steps,1):
        for j in range(matrixParam['components']):
            K1 = system_diff_equation(y[i,1:])[j]
            K2 = system_diff_equation(y[i,1:] + h * K1 / 2)[j]
            K3 = system_diff_equation(y[i,1:] + h * K2 / 2)[j]
            K4 = system_diff_equation(y[i,1:] + h * K3)[j]
            y[i+1,j+1] =  y[i,j+1] + (h / 6)  * (K1 + 2*K2 + 2*K3 + K4)
    return y


def experients_table():
    dictionary = array_to_matrix_dictionary()
    return dictionary['expDat']


