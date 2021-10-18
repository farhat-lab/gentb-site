"""
A collection of useful report classes for pipelines.
"""


class PipelineProgramReport:
    completed = property(lambda self: sum([bool(run.is_complete) for run in self.runs] + [0]))
    not_submitted = property(lambda self: sum([not bool(run.is_submitted) for run in self.runs] + [0]))
    not_started = property(lambda self: sum([not bool(run.is_started) for run in self.runs] + [0]))
    no_errors = property(lambda self: sum([not bool(run.is_error) for run in self.runs] + [0]))

    def __init__(self, version, runs):
        self.version = version
        self.runs = list(runs.order_by('-pk')[:30])
        self.count = len(self.runs)

    @property
    def percent_complete(self):
        return "%.1f" % (self.no_errors / self.completed)
