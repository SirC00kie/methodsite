from django.urls import path

from .views import index, MatrixCreateView, tables, result


urlpatterns =[
    path('result/', result, name='result'),
    path('tables/', tables, name='tables'),
    path('add/', MatrixCreateView.as_view(), name='add'),
    path('', index, name='index'),
]