from django.db import models

class Matrix(models.Model):
    stages = models.IntegerField(verbose_name="Количество стадий", null = False, blank = False)
    components = models.IntegerField(verbose_name="Количество компонентов", null=False, blank=False)
    experients = models.IntegerField(verbose_name="Количество опытов", null=False, blank=False)
    step = models.FloatField(verbose_name="Шаг", null=False, blank=False)
    class Meta:
        verbose_name = "Матрицы"
        verbose_name_plural = "Матрицы"