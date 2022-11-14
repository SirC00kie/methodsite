import time

from django.shortcuts import render, redirect
from django.views.generic.edit import CreateView
from django.urls import reverse_lazy, reverse
from methods.Calculation.method_calculator import MethodCalculator

from .file_service import write_to_file, read_from_file, get_json
from .methods_service import BaseMethods

from .forms import MatrixForm
from .models import Matrix


def index(request):
    return render(request, 'methods/index.html')


def tables(request):
    if request.method == 'POST':
        json_tables = request.POST.get('json_tables')
        write_to_file(json_tables, 'methods/static/methods/json/tables_base.json')
        return redirect(reverse('result'))
    else:
        matrix = Matrix.objects.all()
        context = {'matrix': matrix}
        return render(request, 'methods/tables.html', context)


def result(request):
    tables = get_json("methods/static/methods/json/tables_base.json")
    expData = tables["Экспериментальные данные"]
    methodName = tables["Название метода"]
    baseMethod = BaseMethods(tables)

    baseMethod.calculation_all_methods()
    baseMethod.calculation_all_experients()
    baseMethod.calculation_relative_error()

    opt_calculate = MethodCalculator(tables)
    print(opt_calculate.calculation_for_optimization())


    allMethodsData = read_from_file('methods/static/methods/json/all_methods_result.json')
    methods = get_json("methods/static/methods/json/all_methods_result.json")
    selectedMethodData = methods.get(methodName)
    allMethodsExpData = read_from_file("methods/static/methods/json/all_methods_exp.json")
    allMethodsRelativeError = read_from_file("methods/static/methods/json/all_methods_relative_error.json")

    time.sleep(5)
    context = {
        'expData': expData,
        'methodName': methodName,
        'allMethodsData': allMethodsData,
        'selectedMethodData': selectedMethodData,
        'allMethodsExpData': allMethodsExpData,
        'allMethodsRelativeError': allMethodsRelativeError,
    }
    return render(request, 'methods/result.html', context)


class MatrixCreateView(CreateView):
    template_name = 'methods/create.html'
    form_class = MatrixForm
    success_url = reverse_lazy('tables')

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        table_param = Matrix.objects.all()
        table_param.delete()
        return context
