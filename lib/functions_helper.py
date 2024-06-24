import glob
import os
import time

def get_most_recent_log_time(
) -> float:
    """
    Parses the name of the most recently modified snakemake log file
    in order to retrieve its creation time.

    Returns
    -------
    float
        The time indicated in the most recent log file in epoch seconds
    """

    # Get the most recently modified log file
    last_log = max(glob.glob(os.path.join('.snakemake', 'log', '*.log')),
        key=os.path.getmtime)
    # Parse the name to get the datetime string
    last_time = last_log.split('/')[-1].split('.')[0]
    # Convert the string to epoch time
    last_time_f = time.mktime(time.strptime(last_time, '%Y-%m-%dT%H%M%S'))

    return last_time_f