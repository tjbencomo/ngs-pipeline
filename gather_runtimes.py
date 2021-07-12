"""
Parse log files to determine average runtime for each pipeline step
Instructions:
python gather_runtimes.py [log directory]
"""

import os
import re
import sys
from datetime import datetime as dt
import pandas as pd

def parse_log(fp):
    f = open(fp, 'r')
    lines = f.readlines()
    for i, line in enumerate(lines):
        if 'rule ' in line:
            rule_line = i
        elif 'Finished job' in line:
            finish_line = i
    start_time = lines[rule_line - 1].rstrip()
    start_time = start_time[1:len(start_time)-1]
    end_time = lines[finish_line - 1].rstrip()
    end_time = end_time[1:len(end_time)-1]
    rule = re.split(' |:', lines[rule_line])[1]
    start_time = dt.strptime(start_time, '%a %b %d %H:%M:%S %Y')
    end_time = dt.strptime(end_time, '%a %b %d %H:%M:%S %Y')
    d = end_time - start_time
    return rule, d.total_seconds()

def main():
    logdir = sys.argv[1]
    if not os.path.isdir(logdir):
        raise ValueError(f"{logdir} doesn't exist!")
    files = os.listdir(logdir)
    print(f"Found {len(files)} log files")
    rules = []
    seconds = []
    for fn in files:
        rule, time = parse_log(os.path.join(logdir, fn))
        rules.append(rule)
        seconds.append(time)

    df = pd.DataFrame({'rule' : rules, 'seconds' : seconds})
    seconds_per_hour = 60 ** 2
    df['hours'] = df['seconds'] / seconds_per_hour
    mean_times = df.groupby('rule').agg({'seconds' : 'max'})
    mean_times['hours'] = mean_times['seconds'] / seconds_per_hour
    mean_times = mean_times.sort_values(by=['hours'])
    print(mean_times)
    mean_times.to_csv('max_rule_times.csv')

if __name__ == '__main__':
    main()
