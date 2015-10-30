from apps.predict.models import PredictDatasetStatus, PredictDataset, PredictDatasetNote, DatasetScriptRun, ScriptToRun
from apps.script_helper.script_runner_basic import run_script

class ErrMsg:
    def __init__(self, title, note):
        self.title = title
        self.note = note


def run_script_on_dataset(dataset):
    assert dataset is not None, 'The dataset cannot be None'

    # (1) get the script
    try:
        chosen_script = ScriptToRun.objects.get(is_chosen_script=True)
    except ScriptToRun.DoesNotExist:
        print("No ScriptToRun found!  Add a ScriptToRun through the admin. See apps/predict/script_run_helper.py")
        # Run failed, set status
        dataset.set_status_processing_failed()

        err_title='No chosen script'
        err_note = """Please go into the admin control panel and add a 'ScriptToRun'.  Talk to your administrator for details."""

        note = PredictDatasetNote(dataset=dataset,
                              title=err_title,
                              note=err_note)
        note.save()
        return (False, ErrMsg(err_title, err_note))

    # (2) Create a run object
    #
    dsr = DatasetScriptRun(dataset=dataset)
    dsr.save()

    # (2) Format actual command line args
    #
    json_args = dataset.get_script_args_json(run_md5=dsr.md5, as_list=True)

    # (3) Combine script + dataset args
    #
    full_args = chosen_script.get_script_args_as_list() + json_args
    print ('full_args', full_args)

    # (4) Save the command being run
    #
    dsr.notes='%s' % ' '.join(full_args)
    dsr.save()

    # (5) Update the dataset status to 'in process'
    #
    dataset.set_status_processing_started()

    # (6) Run the script -- not waiting for output
    #
    run_script(full_args)


    return (True, dsr)
