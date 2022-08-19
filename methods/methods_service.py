import numpy as np
from .models import Matrix
from .file_service import read_table_to_array


class FunctionCalculate:
    array = read_table_to_array('methods/static/methods/json/tables.json')
    matrixParam = dict()
    table = dict()

    matrixParam['stages'] = Matrix._meta.get_field('stages').value_from_object(Matrix.objects.first())
    matrixParam['components'] = Matrix._meta.get_field('components').value_from_object(Matrix.objects.first())
    matrixParam['experients'] = Matrix._meta.get_field('experients').value_from_object(Matrix.objects.first())
    matrixParam['step'] = Matrix._meta.get_field('step').value_from_object(Matrix.objects.first())

    table['stehCoef'] = array[:matrixParam['stages'] * matrixParam['components']]
    table['pokazStep'] = array[
                              matrixParam['stages'] * matrixParam['components']:matrixParam['stages'] * matrixParam[
                                  'components'] * 2]
    table['expDat'] = array[
                           matrixParam['stages'] * matrixParam['components'] * 2: matrixParam['stages'] * matrixParam[
                               'components'] * 2 + matrixParam['experients'] * (matrixParam['components'] + 1)]
    table['constSpeed'] = array[
                               matrixParam['stages'] * matrixParam['components'] * 2 + matrixParam['experients'] * (
                                       matrixParam['components'] + 1):]

    def system_diff_equation(self, y):

        A = np.zeros(self.matrixParam['components'])
        sumA = np.zeros(self.matrixParam['components'])
        r = np.zeros(self.matrixParam['stages'])
        sumR = np.zeros(self.matrixParam['stages'])

        for i in range(0, self.matrixParam['stages']):
            sumR[i] = self.table['constSpeed'][i]
            for j in range(0, self.matrixParam['components']):
                if self.table['pokazStep'][j + (i * self.matrixParam['components'])] != 0:
                    r[i] = y[j] ** self.table['pokazStep'][
                        j + (i * self.matrixParam['components'])]
                    sumR[i] *= r[i]

        for i in range(0, self.matrixParam['components']):
            for j in range(0, self.matrixParam['stages']):
                A[i] = self.table['stehCoef'][i + (j * self.matrixParam['components'])] * sumR[j]
                sumA[i] += A[i]
        return sumA


class BaseMethods(FunctionCalculate):
    def __init__(self):
        self.t = self.table['expDat'][1]
        self.h = self.matrixParam['step']
        self.n = int(self.table['expDat'][
                         (self.matrixParam['experients'] * (self.matrixParam['components'] + 1) - (
                                 self.matrixParam['components'] + 1))])
        self.steps = int(self.n / self.h)

    def method(self):
        x = 0
        y = np.zeros((self.steps + 1, self.matrixParam['components'] + 1))
        for i in range(self.steps + 1):
            y[i, 0] = x
            x += self.h

        for j in range(self.matrixParam['components']):
            y[0, j + 1] = self.table['expDat'][j + 1]
        return y

    def method_euler(self):
        y = self.method()
        for i in range(0, self.steps, 1):
            for j in range(self.matrixParam['components']):
                y[i + 1, j + 1] = y[i, j + 1] + self.h * self.system_diff_equation(y[i, 1:])[j]
        return y

    def method_runge_kutta_second(self):
        y = self.method()
        for i in range(0, self.steps, 1):
            for j in range(self.matrixParam['components']):
                y[i + 1, j + 1] = y[i, j + 1] + (self.h / 2) * \
                                  (self.system_diff_equation(y[i, 1:])[j] +
                                   self.system_diff_equation(
                                       y[i, 1:] + self.h * self.system_diff_equation(y[i, 1:])[j])[j])
        return y

    def method_runge_kutta_fourth(self):
        y = self.method()
        for i in range(0, self.steps, 1):
            for j in range(self.matrixParam['components']):
                K1 = self.system_diff_equation(y[i, 1:])[j]
                K2 = self.system_diff_equation(y[i, 1:] + self.h * K1 / 2)[j]
                K3 = self.system_diff_equation(y[i, 1:] + self.h * K2 / 2)[j]
                K4 = self.system_diff_equation(y[i, 1:] + self.h * K3)[j]
                y[i + 1, j + 1] = y[i, j + 1] + (self.h / 6) * (K1 + 2 * K2 + 2 * K3 + K4)
        return y

