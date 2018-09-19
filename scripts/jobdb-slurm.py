#!/usr/bin/env python

import sys, os


period = sys.argv[1:]

if len(period) == 0:
    period = ['2018-03-18', '2018-04-11']
sep = ","
csvfile = "jobs.csv"
sqlfile = "create_db.sql"
dbfile = "jobs.db"
stmts = """
CREATE TABLE {0} (
    jobid,
    jobname,
    submit,
    start,
    end,
    mem,
    reqcpus,
    nodelist,
    state
);

.separator {1}
.import {2} {0}

DELETE FROM {0} WHERE jobid = 'JobID';
-- post-processing {0} table (assumes submission via xenon-cli)
UPDATE {0} SET jobid=REPLACE(jobid, '.batch', '');
UPDATE {0} SET jobname=REPLACE(jobname, 'batch', '');

CREATE TABLE {3} AS
SELECT
    jobid,
    GROUP_CONCAT(jobname, '') AS jobname,
    start,
    end,
    (CAST(REPLACE(mem, 'K', '') AS INTEGER) / 1024) AS mem_mb,
    reqcpus,
    nodelist,
    state
FROM {0}
GROUP BY jobid;

CREATE VIEW {4} AS
SELECT
    jobid,
    jobname,
    strftime('%Y-%m-%d %H:%M:%D', start) AS stime,
    strftime('%Y-%m-%d %H:%M:%D', end) AS etime,
    strftime('%s',end) - strftime('%s', start) AS runtime, -- in sec
    strftime('%Y-%m-%d %H:00:00', start) AS stime_bin,
    strftime('%Y-%m-%d %H:00:00', end) AS etime_bin,
    mem_mb,
    reqcpus AS n_cores,
    nodelist AS hostname,
    state AS status
FROM {3};

CREATE VIEW {5} AS
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
FROM {4}
GROUP BY jobname, status;

CREATE VIEW {6} AS
SELECT
    etime_bin,
    status,
    COUNT(DISTINCT hostname) AS n_hosts,
    SUM(n_cores) AS n_cores,
    COUNT(*) AS n_jobs
FROM {4}
GROUP BY etime_bin, status ORDER BY 1;

CREATE VIEW {7} AS
SELECT T1.etime_bin, T1.status, T1.n_jobs, SUM(T2.n_jobs) AS cum_jobs
FROM {6} AS T1, {6} AS T2
WHERE T2.rowid <= T1.rowid AND T2.status = T1.status
GROUP BY T1.status, T1.etime_bin, T1.n_jobs
ORDER BY T1.status, T1.rowid;
""".format('TMP', sep, csvfile, 'JOB', 'V_JOB', 'VV_JOB', 'V_JOB_E', 'VV_JOB_E')

with open(sqlfile, "w") as fout:
    fout.write(stmts)

cmd = "sacct -u akuzniar -r normal -S {} -E {} -P --delimiter={} \
      -o jobid,jobname,submit,start,end,maxvmsize,reqcpus,nodelist,state > {} \
      && sqlite3 {} < {}".format(period[0], period[1], sep, csvfile, dbfile, sqlfile)
os.system(cmd)
