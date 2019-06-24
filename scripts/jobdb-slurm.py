#!/usr/bin/env python

import os
import sys
import time
import datetime
import sqlite3 as sqlt

user, start, end = sys.argv[1:4] # start/end  in YY-MM-DD
os.environ["SLURM_TIME_FORMAT"] = "standard"
sep = ","
csvfile = "jobs.csv"
dbfile = "jobs.db"
stmt = """
-- post-processing {0} table
-- assumes submission via xenon-cli
UPDATE {0} SET jobid=REPLACE(jobid, '.batch', '');
UPDATE {0} SET jobname=REPLACE(jobname, 'batch', '');

CREATE TABLE {1} AS
SELECT
    CAST(jobid AS INTEGER) AS jobid,
    GROUP_CONCAT(jobname, '') AS jobname,
    submit,
    start,
    end,
    to_sec(totalcpu) AS cputime, -- in sec
    MAX(CAST(REPLACE(maxvmsize, 'K', '') AS INTEGER) / 1024) AS mem_mb,
    MIN(CAST(reqcpus AS INTEGER)) AS reqcpus,
    nodelist,
    state
FROM {0}
GROUP BY jobid;

CREATE INDEX idx_{1}_jobid ON {1}(jobid);
CREATE INDEX idx_{1}_jobname ON {1}(jobname);

CREATE VIEW {2} AS
SELECT
    jobid,
    jobname,
    strftime('%Y-%m-%d %H:%M:%S', submit) AS subtime,
    strftime('%Y-%m-%d %H:%M:%S', start) AS stime,
    strftime('%Y-%m-%d %H:%M:%S', end) AS etime,
    strftime('%s', start) - strftime('%s', submit) AS qtime, -- in sec
    strftime('%s', end) - strftime('%s', start) AS runtime,  -- in sec
    cputime, -- in sec
    strftime('%Y-%m-%d %H:00:00', start) AS stime_bin,
    strftime('%Y-%m-%d %H:00:00', end) AS etime_bin,
    mem_mb,
    reqcpus AS n_cores,
    nodelist AS hostname,
    state AS status
FROM {1};

CREATE VIEW {3} AS
SELECT
    jobname,
    COUNT(jobid) AS n_jobs,
    status,
    MIN(runtime) AS min_runtime,
    MAX(runtime) AS max_runtime,
    CAST(ROUND(AVG(runtime)) AS INTEGER) AS mean_runtime,
    MIN(mem_mb) AS min_mem,
    MAX(mem_mb) AS max_mem,
    CAST(ROUND(AVG(mem_mb)) AS INTEGER) AS mean_mem
FROM {2}
GROUP BY jobname, status;

CREATE VIEW {4} AS
SELECT
    etime_bin,
    status,
    COUNT(DISTINCT hostname) AS n_hosts,
    SUM(n_cores) AS n_cores,
    COUNT(*) AS n_jobs
FROM {3}
GROUP BY etime_bin, status ORDER BY 1;

CREATE VIEW {5} AS
SELECT T1.etime_bin, T1.status, T1.n_jobs, SUM(T2.n_jobs) AS cum_jobs
FROM {3} AS T1, {3} AS T2
WHERE T2.rowid <= T1.rowid AND T2.status = T1.status
GROUP BY T1.status, T1.etime_bin, T1.n_jobs
ORDER BY T1.status, T1.rowid;
""".format('TMP', 'JOB', 'V_JOB', 'VV_JOB', 'V_JOB_E', 'VV_JOB_E')

sacct = """
sacct -u {} -S {} -E {} -P --delimiter={} \
-o jobid,jobname,submit,start,end,totalcpu,maxvmsize,reqcpus,nodelist,state > {}
""".format(user, start, end, sep, csvfile)
os.system(sacct)

sqlite = "sqlite3 {} -csv '.import {} TMP'".format(dbfile, csvfile)
os.system(sqlite)

def to_sec(timestr):
    t, sec, days, hms, ms = (None, 0, 0, '00:00:00', '00:00.000')
    if '-' in timestr:
        days, hms = timestr.split('-')
        days = int(days)
        t = time.strptime(hms, "%H:%M:%S")
    elif '.' in timestr:
        t = time.strptime(timestr, "%M:%S.%f")
    else:
        t = time.strptime(timestr, "%H:%M:%S")
    sec = int(datetime.timedelta(days=days,
                                 hours=t.tm_hour,
                                 minutes=t.tm_min,
                                 seconds=t.tm_sec).total_seconds())
    return int(sec)

with sqlt.connect(dbfile) as conn:
   conn.row_factory = sqlt.Row # enable column access by name: row['colnm']
   conn.create_function('to_sec', 1, to_sec)
   cur = conn.cursor()
   cur.executescript(stmt)
