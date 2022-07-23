from django.forms import ModelForm

from .models import Matrix

class MatrixForm(ModelForm):
    class Meta:
        model = Matrix
        fields = {'stages', 'components','experients', 'step'}