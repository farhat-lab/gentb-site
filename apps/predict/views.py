"""
Views for predict app
"""
from __future__ import print_function

import logging

LOGGER = logging.getLogger(__name__)

from django.views.generic import (
    DetailView, ListView, CreateView, UpdateView, FormView,
    TemplateView,
)
from django.views.generic.detail import SingleObjectMixin
from django.urls import reverse
from django.http.response import JsonResponse
from django.http import Http404

from apps.maps.mixins import JsonView
from apps.mutations.models import GeneLocus

from .models import PredictDataset, PredictResult, PredictResultLocus, PredictDatasetNote
from .mixins import PredictMixin
from .forms import *


class Datasets(PredictMixin, ListView):
    title = "My Datasets"

    @classmethod
    def get_absolute_url(cls):
        return reverse('predict:view_my_datasets')


class DatasetView(PredictMixin, DetailView):
    parent = Datasets

class DatasetViewProcessing(DatasetView):
    template_name = 'predict/predictdataset_processing.html'

class DatasetViewOutput(DatasetView):
    template_name = 'predict/predictdataset_outputdata.html'

class DatasetViewPredict(DatasetView):
    template_name = 'predict/predictdataset_prediction.html'

class DatasetViewLineages(DatasetView):
    template_name = 'predict/predictdataset_lineages.html'


class Heatmap(PredictMixin, DetailView):
    queryset = PredictDataset.objects.all()
    template_name = 'predict/heatmap.html'

class UploadChoices(PredictMixin, TemplateView):
    template_name = 'predict/predictdataset_upload.html'
    forms = UploadForm.all_forms()
    title = "Create Prediction"
    parent = Datasets

    @classmethod
    def get_absolute_url(cls):
        return reverse('predict:upload')


class UploadView(PredictMixin, CreateView):
    model = PredictDataset
    parent = UploadChoices

    def get_title(self):
        return self.form_class.doc_title

    def get_template_names(self):
        default = super(UploadView, self).get_template_names()
        default = ['predict/predictdataset_%s.html' % self.kwargs['type']] + default
        return default

    @property
    def form_class(self):
        try:
            return UploadForm.get_form(self.kwargs['type'])
        except KeyError:
            raise Http404("No input type: %s" % self.kwargs['type'])

    def get_initial(self):
        return {
          'user': self.request.user,
          'file_type': self.kwargs['type'],
        }


class AddNote(PredictMixin, CreateView):
    model = PredictDatasetNote
    fields = ('note',)

    def form_invalid(self, form):
        return JsonResponse({'status': 'INVALID'})

    def form_valid(self, form):
        obj = form.save(commit=False)
        obj.title = str(self.request.user)
        obj.dataset = self.get_queryset().get(md5=self.kwargs['slug'])
        obj.save()
        return JsonResponse({
          'status': 'OK',
          'title': obj.title,
          'note': obj.note,
        })


class ScatterPlot(JsonView, SingleObjectMixin):
    model = PredictResult

    def build_loci_list(self, drug, dataset):
        """
        Compile a list of all loci for this drug used in the dataset.
        """
        # First compile a list used by the drug
        loci = list(drug.loci.all())
        # Next find all results for this drug for the ENTIRE dataset
        all_results = PredictResultLocus.objects.filter(result__drug=drug, result__strain__dataset=dataset)
        # Don't use it here! send the SQL back to the database to find the loci used
        for locus in GeneLocus.objects.filter(pk__in=all_results.values('locus_id')):
            if locus not in loci:
                loci.append(locus)
        return loci

    def get_context_data(self, **kwargs):
        result = self.get_object()
        loci = self.build_loci_list(result.drug, result.strain.dataset)

        ret = []
        for label, color in (
            ("Important", "255, 0, 0, 0.8"),
            ("Other", "0, 0, 255, 0.17"),
            ("New", "255, 255, 0, 0.5"),
            ("Lineage SNPs", "0, 255, 255, 0.5"),
        ):
            ret.append({
                "yAxis": "1",
                "cols": [str(l) for l in loci],
                "key": label,
                "color": "rgba(%s)" % color,
                "values": [
                    {"x": x, "y": 0, "size": 5, "tip": [f"No {l} mutations"]}
                    for x, l in enumerate(loci)
                ],
            })

        for item in result.loci.all():
            plot = ret[item.category - 1]
            locus_index = loci.index(item.locus)
            row = plot["values"][locus_index]
            row["size"] = 9
            row["tip"] = item.mutations.split("\n")
            row["y"] = len(row["tip"])

        return {'data': ret}
