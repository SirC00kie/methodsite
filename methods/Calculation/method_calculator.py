from methods.methods_service import BaseMethods


class MethodCalculator(BaseMethods):
    def __init__(self, arr):
        super().__init__(arr)

    def calculation_for_optimization(self):
        speed_const = [[13], [3.5], [3.5]]
        return self.current_method(speed_const)


