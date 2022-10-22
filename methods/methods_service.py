import numpy as np
from methods.file_service import set_json, get_json

from .models import Matrix


class BaseMethods:
    matrixParam = dict()

    def __init__(self, arr):
        self.tables = arr
        self.matrixParam['stages'] = Matrix._meta.get_field('stages').value_from_object(Matrix.objects.first())
        self.matrixParam['components'] = Matrix._meta.get_field('components').value_from_object(Matrix.objects.first())
        self.matrixParam['experients'] = Matrix._meta.get_field('experients').value_from_object(Matrix.objects.first())
        self.matrixParam['step'] = Matrix._meta.get_field('step').value_from_object(Matrix.objects.first())
        self.stehCoef = self.tables["Матрица стехиометрических коэффициентов"]
        self.pokazStep = self.tables["Матрица показателей степени"]
        self.expDat = self.tables["Экспериментальные данные"]
        self.constSpeed = self.tables["Константы скорости"]
        self.selectedMethod = self.tables["Название метода"]
        self.t = self.expDat[0][0]
        self.h = self.matrixParam['step']
        self.n = self.expDat[self.matrixParam['experients'] - 1][0]
        self.steps = int(self.n / self.h)

    def system_diff_equation(self, y):
        A = np.zeros(self.matrixParam['components'])
        sumA = np.zeros(self.matrixParam['components'])
        r = np.zeros(self.matrixParam['stages'])
        sumR = np.zeros(self.matrixParam['stages'])

        for i in range(self.matrixParam['stages']):
            sumR[i] = self.constSpeed[i][0]
            for j in range(self.matrixParam['components']):
                if self.pokazStep[i][j] != 0:
                    r[i] = y[j] ** self.pokazStep[i][j]
                    sumR[i] *= r[i]
        for i in range(self.matrixParam['components']):
            for j in range(self.matrixParam['stages']):
                A[i] = self.stehCoef[j][i] * sumR[j]
                sumA[i] += A[i]
        return sumA

    def method(self):
        x = self.t
        y = np.zeros((self.steps + 1, self.matrixParam['components'] + 1))

        for i in range(self.steps + 1):
            y[i, 0] = x
            x += self.h

        for j in range(self.matrixParam['components'] + 1):
            y[0, j] = self.expDat[0][j]
        return y

    def method_euler(self):
        y = self.method()

        for i in range(self.steps):
            for j in range(self.matrixParam['components']):
                y[i + 1, j + 1] = y[i, j + 1] + self.h * self.system_diff_equation(y[i, 1:])[j]
        return y

    def implicit_method_euler(self):
        y = self.method()
        y_n = self.method()

        for i in range(self.steps):
            for j in range(self.matrixParam['components']):
                y_n[i + 1, j + 1] = y[i, j + 1] + self.h * self.system_diff_equation(y[i, 1:])[j]
                y[i + 1, j + 1] = y[i, j + 1] + self.h * self.system_diff_equation(y_n[i, 1:])[j]
        return y

    def method_trapezoid(self):
        y = self.method()
        y_n = self.method()

        for i in range(self.steps):
            for j in range(self.matrixParam['components']):
                y_n[i + 1, j + 1] = y[i, j + 1] + self.h * self.system_diff_equation(y[i, 1:])[j]
                y[i + 1, j + 1] = y[i, j + 1] + (self.h / 2) * (self.system_diff_equation(y_n[i + 1, 1:])[j] + self.system_diff_equation(y[i, 1:])[j])
        return y

    def method_middle_point(self):
        y = self.method()
        y_n = self.method()

        for i in range(self.steps):
            for j in range(self.matrixParam['components']):
                y_n[i + 1, j + 1] = y[i, j + 1] + self.h * self.system_diff_equation(y[i, 1:])[j]
                y[i + 1, j + 1] = y[i, j + 1] + self.h * self.system_diff_equation((y_n[i + 1, 1:] + y[i, 1:]) / 2)[j]
        return y

    def method_runge_kutta_second(self):
        y = self.method()
        for i in range(self.steps):
            for j in range(self.matrixParam['components']):
                y[i + 1, j + 1] = y[i, j + 1] + (self.h / 2) * \
                                  (self.system_diff_equation(y[i, 1:])[j] +
                                   self.system_diff_equation(
                                       y[i, 1:] + self.h * self.system_diff_equation(y[i, 1:])[j])[j])
        return y

    def method_runge_kutta_fourth(self):
        y = self.method()
        for i in range(self.steps):
            for j in range(self.matrixParam['components']):
                K1 = self.system_diff_equation(y[i, 1:])[j]
                K2 = self.system_diff_equation(y[i, 1:] + self.h * K1 / 2)[j]
                K3 = self.system_diff_equation(y[i, 1:] + self.h * K2 / 2)[j]
                K4 = self.system_diff_equation(y[i, 1:] + self.h * K3)[j]
                y[i + 1, j + 1] = y[i, j + 1] + (self.h / 6) * (K1 + 2 * K2 + 2 * K3 + K4)
        return y

    def calculation_all_methods(self):
        method_euler = self.method_euler()
        implicit_method_euler = self.implicit_method_euler()
        method_trapezoid = self.method_trapezoid()
        method_middle_point = self.method_middle_point()
        method_runge_kutta_second = self.method_runge_kutta_second()
        method_runge_kutta_fourth = self.method_runge_kutta_fourth()

        tdict = {
            "метод Эйлера": method_euler.tolist(),
            "Неявный метод Эйлера": implicit_method_euler.tolist(),
            "Метод трапеций": method_trapezoid.tolist(),
            "Метод средней точки": method_middle_point.tolist(),
            "Метод Рунге-Кутты 2-го порядка": method_runge_kutta_second.tolist(),
            "Метод Рунге-Кутты 4-го порядка": method_runge_kutta_fourth.tolist(),
        }
        set_json(tdict, "methods/static/methods/json/all_methods_result.json")

    def experient(self, name: str):
        methods = get_json("methods/static/methods/json/all_methods_result.json")
        mas = self.expDat.copy()
        j = 0
        for i in range(len(methods[name])):
            if round(self.expDat[j][0], 1) == round(methods[name][i][0], 1):
                mas[j] = methods[name][i][:]
                j += 1
        return mas

    def calculation_all_experients(self):

        tdict = {
            "метод Эйлера экспериментальные": self.experient("метод Эйлера"),
            "Неявный метод Эйлера экспериментальные": self.experient("Неявный метод Эйлера"),
            "Метод трапеций экспериментальные": self.experient("Метод трапеций"),
            "Метод средней точки экспериментальные": self.experient("Метод средней точки"),
            "Метод Рунге-Кутты 2-го порядка экспериментальные": self.experient("Метод Рунге-Кутты 2-го порядка"),
            "Метод Рунге-Кутты 4-го порядка экспериментальные": self.experient("Метод Рунге-Кутты 4-го порядка"),
        }
        set_json(tdict, "methods/static/methods/json/all_methods_exp.json")
