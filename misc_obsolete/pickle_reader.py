import pandas as pd

def read_pickle_file(file):
    pickle_data = pd.read_pickle("/project/berglandlab/DEST/dest_mapped/pipeline_output/PA_li_09_fall/PA_li_09_fall.cov")
    return pickle_data
