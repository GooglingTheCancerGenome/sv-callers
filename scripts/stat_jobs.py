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
CREATE TABLE {0} (jobid, jobname, start, end, mem, state);
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
    state
FROM {0}
GROUP BY jobid;

CREATE VIEW {4} AS
SELECT
    jobid,
    jobname,
    strftime('%Y-%m-%d %H:%M:%D', start) AS start,
    strftime('%Y-%m-%d %H:%M:%D', end) AS end,
    strftime('%s',end) - strftime('%s', start) AS rtime_s,
    strftime('%Y-%m-%d %H:00:00', start) AS start_bin,
    strftime('%Y-%m-%d %H:00:00', end) AS end_bin,
    mem_mb,
    state
FROM {3};

CREATE VIEW {5} AS
SELECT
    jobname,
    COUNT(jobid) AS n_jobs,
    state,
    MIN(rtime_s) AS min_rtime,
    MAX(rtime_s) AS max_rtime,
    CAST(ROUND(AVG(rtime_s)) AS INT) AS mean_rtime,
    MIN(mem_mb) AS min_mem,
    MAX(mem_mb) AS max_mem,
    CAST(ROUND(AVG(mem_mb)) AS INT) AS mean_mem
FROM {4}
GROUP BY jobname, state;

CREATE VIEW {6} AS
SELECT end_bin, state, COUNT(*) AS n_jobs
FROM {4}
GROUP BY end_bin, state ORDER BY 1;

CREATE VIEW {7} AS
SELECT T1.end_bin, "end", T1.state, T1.n_jobs, SUM(T2.n_jobs) AS cum_jobs
FROM {6} AS T1, V_JOB_E AS T2
WHERE T2.rowid <= T1.rowid AND T2.state = T1.state
GROUP BY T1.state, T1.end_bin, T1.n_jobs
ORDER BY T1.state, T1.rowid;
""".format('TMP', sep, csvfile, 'JOB', 'V_JOB', 'VV_JOB', 'V_JOB_E', 'VV_JOB_E')

with open(sqlfile, "w") as fout:
    fout.write(stmts)

cmd = "sacct -u akuzniar -r normal -S {} -E {} -P --delimiter={} \
      -o jobid,jobname,start,end,maxvmsize,state > {} && \
      sqlite3 {} < {}".format(period[0], period[1], sep, csvfile, dbfile, sqlfile)
os.system(cmd)
