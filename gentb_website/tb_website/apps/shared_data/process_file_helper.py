from apps.script_helper.rscript_runner import RScriptRunner



def get_process_file_results(fullpath_to_file):
    
    runner = RScriptRunner(fullpath_to_file)
    runner.run_command_on_file()
    
    if runner.err_found:
        return (False, runner.get_err_msgs())
        
    return (True, runner.get_formatted_response())
    

