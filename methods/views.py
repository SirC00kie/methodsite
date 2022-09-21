from django.shortcuts import render, redirect
from django.views.generic.edit import CreateView
from django.urls import reverse_lazy, reverse
import json

from .file_service import write_table_to_file, NumpyEncoder, write_to_file, read_from_file, read_table_to_array
from .methods_service import BaseMethods

from .forms import MatrixForm
from .models import Matrix


def index(request):
    return render(request, 'methods/index.html')


def tables(request):
    if request.method == 'POST':
        text = request.POST.get('data_table')
        method_name = request.POST.get('method_name')
        write_to_file(method_name, 'methods/static/methods/json/method_name.json')
        write_table_to_file(text, 'methods/static/methods/json/tables.json')
        return redirect(reverse('result'))
    else:
        matrix = Matrix.objects.all()
        context = {'matrix': matrix}
        return render(request, 'methods/tables.html', context)


def result(request):
    methodName = read_from_file('methods/static/methods/json/method_name.json')
    array = read_table_to_array('methods/static/methods/json/tables.json')
    baseMethod = BaseMethods(array)
    dataMethod = baseMethod.method_selection(methodName)
    method_name = baseMethod.method_name_selection(methodName)
    exp_table = baseMethod.table['expDat']
    json_dump = json.dumps({'data': dataMethod},
                           cls=NumpyEncoder)
    json_exp = json.dumps({'exp_table': exp_table},
                           cls=NumpyEncoder)
    matrix = Matrix.objects.all()
    context={
        'data': dataMethod,
        'json_dump': json_dump,
        'json_exp': json_exp,
        'matrix': matrix,
        'method_name': method_name,
    }
    return render (request, 'methods/result.html', context)


class MatrixCreateView(CreateView):
    template_name = 'methods/create.html'
    form_class = MatrixForm
    success_url = reverse_lazy('tables')

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        table_param = Matrix.objects.all()
        table_param.delete()
        return context
