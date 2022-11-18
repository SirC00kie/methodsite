from methods.methods_service import BaseMethods


class OptimizationCalculator:

    def __init__(self, arr):
        self.BaseMethods = BaseMethods(arr)

    def calculation_for_optimization(self, speed_const):
        return self.BaseMethods.current_method(speed_const)

    def print_result(self):
        speed_const = [[1.3], [2.5]]
        print(self.calculation_for_optimization(speed_const))
