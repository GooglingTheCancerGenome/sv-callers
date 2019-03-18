#!/usr/bin/env python

from __future__ import print_function
import sys
import math
import numpy as np
import sqlite3 as sqlt
import datetime

# test user input
if len(sys.argv) == 1:
    print("Usage: %s [dbfile]" % sys.argv[0])
    sys.exit(2)

dbfile = sys.argv[1]

# user-defined function (UDF)
class Median:
    def __init__(self):
        self.arr = []

    def step(self, value):
        self.arr.append(value)

    def finalize(self):
        return np.median(np.array(self.arr))

# connect to db
with sqlt.connect(dbfile) as conn:
    conn.row_factory = sqlt.Row  # enable column access by name: row['colnm']
    conn.create_aggregate('median', 1, Median)
    cur = conn.cursor()
    sql_stmt = """
    SELECT
      jobname,
      COUNT(jobid) AS n,
      ROUND(MIN(runtime), 0) AS min_rt, -- in sec
      ROUND(MAX(runtime), 0) AS max_rt,
      ROUND(AVG(runtime), 0) AS mean_rt,
      ROUND(MEDIAN(runtime), 0) AS median_rt,
      ROUND(MIN(cputime), 0) AS min_ct, -- in sec
      ROUND(MAX(cputime), 0) AS max_ct,
      ROUND(AVG(cputime), 0) AS mean_ct,
      ROUND(MEDIAN(cputime), 0) AS median_ct
    FROM V_JOB
    WHERE status = 'COMPLETED'
    GROUP BY jobname"""

    cur.execute(sql_stmt)
    for row in cur.fetchall():
        for k in row.keys():
            v = row[k]
            if k in ('jobname'):
                print("%s: %s" % (k, v))
            else:
                print("%s: %d [%s]" % (k, v, datetime.timedelta(seconds=v)))
        print("\n")
