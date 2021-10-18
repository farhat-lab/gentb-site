"""
A collection of useful report classes for pipelines.
"""


class PipelineProgramReport:
    completed = property(lambda self: sum([bool(run.is_submitted and run.is_started and run.is_complete and not run.is_error) for run in self.runs] + [0]))
    not_submitted = property(lambda self: sum([not bool(run.is_submitted) for run in self.runs] + [0]))
    not_started = property(lambda self: sum([bool(run.is_submitted and not run.is_started) for run in self.runs] + [0]))
    has_errors = property(lambda self: sum([bool(run.is_error) for run in self.runs] + [0]))

    def __init__(self, version, runs):
        self.version = version
        self.runs = list(runs.order_by('-pk')[:30])
        self.count = len(self.runs)

    @property
    def percent_complete(self):
        total = self.count - self.not_started
        if total:
            return "%.1f%%" % ((self.completed / total) * 100)
        return "-"
