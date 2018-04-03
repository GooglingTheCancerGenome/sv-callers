#!/usr/bin/env python

import os

table = "JOB"
sep = ","
csvfile = "jobs.csv"
sqlfile = "create_db.sql"
dbfile = "jobs.db"
stmts = """
CREATE TABLE {0} (`jobid`, `jobname`, `start`, `end`, `state`);
.separator {1}
.import {2} {0}
DELETE FROM {0} WHERE jobid = 'JobID';
""".format(table, sep, csvfile)

with open(sqlfile, "w") as fout:
    fout.write(stmts)

dt = '2018-03-18'
cmd = "sacct -u akuzniar -r normal -s CD -S {} -P --delimiter={} \
      -X -o jobid,jobname,start,end,state > {} && \
      sqlite3 {} < {}".format(dt, sep, csvfile, dbfile, sqlfile)
os.system(cmd)

"""
CREATE VIEW V_JOB AS
SELECT
    `jobid`,
    `jobname`,
    `start`,
    `end`,
    `state`,
    strftime('%s',`end`) - strftime('%s', `start`) AS elapsed,
    strftime("%Y-%m-%d %H:00:00", `start`) AS start_bin,
    strftime("%Y-%m-%d %H:00:00", `end`) AS end_bin
FROM JOB;

SELECT jobname, COUNT(jobid), min(elapsed), max(elapsed), avg(elapsed) FROM V_JOB GROUP BY jobname;

CREATE VIEW V_JOB_S AS
SELECT start_bin, COUNT(*) AS n_jobs
FROM V_JOB
GROUP BY start_bin ORDER BY 1;

CREATE VIEW V_JOB_E AS
SELECT end_bin, COUNT(*) AS n_jobs
FROM V_JOB WHERE state = 'COMPLETED'
GROUP BY end_bin ORDER BY 1;

CREATE VIEW VV_JOB_S AS
SELECT T1.start_bin, "start", T1.n_jobs, SUM(T2.n_jobs)
FROM V_JOB_S AS T1, V_JOB_S AS T2
WHERE T2.rowid <= T1.rowid
GROUP BY T1.start_bin, T1.n_jobs
ORDER BY T1.rowid;

CREATE VIEW VV_JOB_E AS
SELECT T1.end_bin, "end", T1.n_jobs, SUM(T2.n_jobs)
FROM V_JOB_E AS T1, V_JOB_E AS T2
WHERE T2.rowid <= T1.rowid
GROUP BY T1.end_bin, T1.n_jobs
ORDER BY T1.rowid;

SELECT * FROM VV_JOB_S
UNION
SELECT * FROM VV_JOB_E;
"""
