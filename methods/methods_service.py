import math

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
        self.stehCoef = self.tables["Матрица стехиометрических коэффициентов"].copy()
        self.pokazStep = self.tables["Матрица показателей степени"].copy()
        self.expDat = self.tables["Экспериментальные данные"].copy()
        self.constSpeed = self.tables["Константы скорости"].copy()
        self.selectedMethod = self.tables["Название метода"]
        self.t = self.expDat[0][0]
        self.h = self.matrixParam['step']
        self.n = self.expDat[self.matrixParam['experients'] - 1][0]
        self.steps = int(self.n / self.h)

    def system_diff_equation(self, y, const_speed):
        C = np.zeros(self.matrixParam['components'])
        sumC = np.zeros(self.matrixParam['components'])
        r = np.zeros(self.matrixParam['stages'])
        sumR = np.zeros(self.matrixParam['stages'])

        for i in range(self.matrixParam['stages']):
            sumR[i] = const_speed[i][0]
            for j in range(self.matrixParam['components']):
                if self.pokazStep[i][j] != 0:
                    r[i] = y[j] ** self.pokazStep[i][j]
                    sumR[i] *= r[i]
        for i in range(self.matrixParam['components']):
            for j in range(self.matrixParam['stages']):
                C[i] = self.stehCoef[j][i] * sumR[j]
                sumC[i] += C[i]
        return sumC

    def method(self):
        x = self.t
        y = np.zeros((self.steps + 1, self.matrixParam['components'] + 1))

        for i in range(self.steps + 1):
            y[i, 0] = x
            x += self.h

        for j in range(self.matrixParam['components'] + 1):
            y[0, j] = self.expDat[0][j]
        return y

    def method_euler(self, const_speed):
        y = self.method()

        for i in range(self.steps):
            for j in range(self.matrixParam['components']):
                y[i + 1, j + 1] = y[i, j + 1] + self.h * self.system_diff_equation(y[i, 1:], const_speed)[j]
        return y

    def implicit_method_euler(self, const_speed):
        y = self.method()
        y_euler = self.method_euler(const_speed)

        for i in range(self.steps):
            for j in range(self.matrixParam['components']):
                y[i + 1, j + 1] = y[i, j + 1] + self.h * self.system_diff_equation(y_euler[i + 1, 1:], const_speed)[j]
        return y

    def method_trapezoid(self, const_speed):
        y = self.method()
        y_euler = self.method_euler(const_speed)

        for i in range(self.steps):
            for j in range(self.matrixParam['components']):
                y[i + 1, j + 1] = y[i, j + 1] + (self.h / 2) * (
                            self.system_diff_equation(y_euler[i + 1, 1:], const_speed)[j] + self.system_diff_equation(y[i, 1:], const_speed)[j])
        return y

    def method_middle_point(self, const_speed):
        y = self.method()
        y_middle = self.method()

        for i in range(self.steps):
            for j in range(self.matrixParam['components']):
                y_middle[i + 1, j + 1] = y[i, j + 1] + self.h * self.system_diff_equation(y[i, 1:], const_speed)[j]
                y[i + 1, j + 1] = y[i, j + 1] + self.h * self.system_diff_equation((y_middle[i + 1, 1:] + y[i, 1:]) / 2, const_speed)[j]
        return y

    def method_runge_kutta_second(self, const_speed):
        y = self.method()
        for i in range(self.steps):
            for j in range(self.matrixParam['components']):
                y[i + 1, j + 1] = y[i, j + 1] + (self.h / 2) * \
                                  (self.system_diff_equation(y[i, 1:], const_speed)[j] +
                                   self.system_diff_equation(
                                       y[i, 1:] + self.h * self.system_diff_equation(y[i, 1:], const_speed)[j], const_speed)[j])
        return y

    def method_runge_kutta_fourth(self, const_speed):
        y = self.method()
        for i in range(self.steps):
            for j in range(self.matrixParam['components']):
                K1 = self.system_diff_equation(y[i, 1:], const_speed)[j]
                K2 = self.system_diff_equation(y[i, 1:] + self.h * K1 / 2, const_speed)[j]
                K3 = self.system_diff_equation(y[i, 1:] + self.h * K2 / 2, const_speed)[j]
                K4 = self.system_diff_equation(y[i, 1:] + self.h * K3, const_speed)[j]
                y[i + 1, j + 1] = y[i, j + 1] + (self.h / 6) * (K1 + 2 * K2 + 2 * K3 + K4)
        return y

    def method_kutta_merson(self, const_speed):
        y = self.method()
        h = self.h
        R = 0
        e = 1
        for i in range(self.steps):
            for j in range(self.matrixParam['components']):
                if R < e:
                    K1 = self.system_diff_equation(y[i, 1:], const_speed)[j]
                    K2 = self.system_diff_equation(y[i, 1:] + (h / 3 * K1), const_speed)[j]
                    K3 = self.system_diff_equation(y[i, 1:] + (h / 6 * K1) + (h / 6 * K2), const_speed)[j]
                    K4 = self.system_diff_equation(y[i, 1:] + (h / 8 * K1) + (3 * h / 8 * K2), const_speed)[j]
                    K5 = self.system_diff_equation(y[i, 1:] + (h / 2 * K1) - (3 * h / 2 * K3) + (2 * h * K4), const_speed)[j]
                    y[i + 1, j + 1] = y[i, j + 1] + h / 6 * (K1 + 4 * K4 + K5)
                elif R >= e:
                    h = h / 2
                elif R <= e / 64:
                    h = h * 2
        return y

    def rkf45(self, const_speed):
        y = self.method()
        h = self.h
        R = 0
        e = 1
        for i in range(self.steps):
            for j in range(self.matrixParam['components']):
                if R < e:
                    K1 = self.system_diff_equation(y[i, 1:], const_speed)[j]
                    K2 = self.system_diff_equation(y[i, 1:] + (h / 4 * K1), const_speed)[j]
                    K3 = self.system_diff_equation(y[i, 1:] + (3 * h / 32 * K1) + (9 * h / 32 * K2), const_speed)[j]
                    K4 = self.system_diff_equation(
                        y[i, 1:] + (1932 * h / 2197 * K1) - (7200 * h / 2197 * K2) + (7296 * h / 2197 * K3), const_speed)[j]
                    K5 = self.system_diff_equation(
                        y[i, 1:] + (439 * h / 216 * K1) - (8 * h * K2) + (3600 * h / 513 * K3) - (845 * h / 4104 * K4), const_speed)[
                        j]
                    K6 = self.system_diff_equation(
                        y[i, 1:] - (8 * h / 27 * K1) + (2 * h * K2) - (3544 * h / 2565 * K3) + (
                                    1859 * h / 4104 * K4) + (11 * h / 40 * K5), const_speed)[j]
                    y[i + 1, j + 1] = y[i, j + 1] + h * (
                                (16 * K1 / 135) + (6656 * K3 / 12825) + (28561 * K4 / 56430) - (9 * K6 / 50) + (
                                    2 * K5 / 55))

                elif R >= e:
                    h = h / 2
                elif R <= e / 64:
                    h = h * 2
        return y

    def explicit_second_adams_method(self, const_speed):
        y = self.method()
        for i in range(self.steps):
            for j in range(self.matrixParam['components']):
                y[i + 1, j + 1] = y[i, j + 1] + (self.h / 2) * (3 * self.system_diff_equation(y[i, 1:], const_speed)[j]
                                                                - self.system_diff_equation(y[i - 1, 1:], const_speed)[j])

        return y

    #def implicit_second_adams_method(self):

    def calculation_all_methods(self):
        method_euler = self.method_euler(self.constSpeed)
        implicit_method_euler = self.implicit_method_euler(self.constSpeed)
        method_trapezoid = self.method_trapezoid(self.constSpeed)
        method_middle_point = self.method_middle_point(self.constSpeed)
        method_runge_kutta_second = self.method_runge_kutta_second(self.constSpeed)
        method_runge_kutta_fourth = self.method_runge_kutta_fourth(self.constSpeed)
        method_kutta_merson = self.method_kutta_merson(self.constSpeed)
        rkf45 = self.rkf45(self.constSpeed)
        explicit_second_adams_method = self.explicit_second_adams_method(self.constSpeed)

        tdict = {
            "метод Эйлера": method_euler.tolist(),
            "Неявный метод Эйлера": implicit_method_euler.tolist(),
            "Метод трапеций": method_trapezoid.tolist(),
            "Метод средней точки": method_middle_point.tolist(),
            "Метод Рунге-Кутты 2-го порядка": method_runge_kutta_second.tolist(),
            "Метод Рунге-Кутты 4-го порядка": method_runge_kutta_fourth.tolist(),
            "Метод Кутты-Мерсона": method_kutta_merson.tolist(),
            "Метод Рунге-Кутты-Фельберга": rkf45.tolist(),
            "Явный двухшаговый метод Адамса": explicit_second_adams_method.tolist(),
        }
        set_json(tdict, "methods/static/methods/json/all_methods_result.json")

    def experient(self, name: str):
        methods = get_json("methods/static/methods/json/all_methods_result.json")
        mas = self.expDat.copy()

        j = 0
        for i in range(len(methods[name])):
            if round(mas[j][0], 1) == round(methods[name][i][0], 1):
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
            "Метод Кутты-Мерсона экспериментальные": self.experient("Метод Кутты-Мерсона"),
            "Метод Рунге-Кутты-Фельберга экспериментальные": self.experient("Метод Рунге-Кутты-Фельберга"),
            "Явный двухшаговый метод Адамса экспериментальные": self.experient("Явный двухшаговый метод Адамса")
        }
        set_json(tdict, "methods/static/methods/json/all_methods_exp.json")

    def relative_error(self, name: str):
        expDots = get_json("methods/static/methods/json/all_methods_exp.json")
        expErorr = expDots[name].copy()

        for i in range(0, len(expErorr), 1):
            for j in range(1, len(expErorr[i]), 1):
                if i == 0:
                    expErorr[i][j] = 0
                elif self.expDat[i][j] == 0:
                    expErorr[i][j] = 0
                else:
                    expErorr[i][j] = (math.fabs(expDots[name][i][j] - self.expDat[i][j]) / self.expDat[i][j]) * 100
        return expErorr

    def calculation_relative_error(self):
        tdict = {
            "метод Эйлера погрешность": self.relative_error("метод Эйлера экспериментальные"),
            "Неявный метод Эйлера погрешность": self.relative_error("Неявный метод Эйлера экспериментальные"),
            "Метод трапеций погрешность": self.relative_error("Метод трапеций экспериментальные"),
            "Метод средней точки погрешность": self.relative_error("Метод средней точки экспериментальные"),
            "Метод Рунге-Кутты 2-го порядка погрешность": self.relative_error(
                "Метод Рунге-Кутты 2-го порядка экспериментальные"),
            "Метод Рунге-Кутты 4-го порядка погрешность": self.relative_error(
                "Метод Рунге-Кутты 4-го порядка экспериментальные"),
            "Метод Кутты-Мерсона погрешность": self.relative_error("Метод Кутты-Мерсона экспериментальные"),
            "Метод Рунге-Кутты-Фельберга погрешность": self.relative_error(
                "Метод Рунге-Кутты-Фельберга экспериментальные"),
            "Явный двухшаговый метод Адамса погрешность": self.relative_error("Явный двухшаговый метод Адамса экспериментальные")
        }
        set_json(tdict, "methods/static/methods/json/all_methods_relative_error.json")

    def current_method(self, const_speed):
        if self.selectedMethod == "метод Эйлера":
            return self.method_euler(const_speed)
        if self.selectedMethod == "Неявный метод Эйлера":
            return self.implicit_method_euler(const_speed)
        if self.selectedMethod == "Метод трапеций":
            return self.method_trapezoid(const_speed)
        if self.selectedMethod == "Метод средней точки":
            return self.method_middle_point(const_speed)
        if self.selectedMethod == "Метод Рунге-Кутты 2-го порядка":
            return self.method_runge_kutta_second(const_speed)
        if self.selectedMethod == "Метод Рунге-Кутты 4-го порядка":
            return self.method_runge_kutta_fourth(const_speed)
        if self.selectedMethod == "Метод Кутты-Мерсона":
            return self.method_kutta_merson(const_speed)
        if self.selectedMethod == "Метод Рунге-Кутты-Фельберга":
            return self.rkf45(const_speed)
        if self.selectedMethod == "Явный двухшаговый метод Адамса":
            return self.explicit_second_adams_method(const_speed)
