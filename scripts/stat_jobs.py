#!/usr/bin/env python

import os

table = "JOB"
sep = ","
csvfile = "jobs.csv"
sqlfile = "create_db.sql"
dbfile = "jobs.db"
stmts = """
CREATE TABLE {} (`jobid`, `jobname`, `start`, `end`);
.separator {}
.import {} {}
""".format(table, sep, csvfile, table)

with open(sqlfile, "w") as fout:
    fout.write(stmts)

cmd = "sacct -u akuzniar -r normal -s CD -S 2018-03-18 -P --delimiter={} \
      -X -o jobid,jobname,start,end > {} && \
      sqlite3 {} < {}".format(sep, csvfile, dbfile, sqlfile)
os.system(cmd)
