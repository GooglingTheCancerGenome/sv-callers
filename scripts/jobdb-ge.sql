
CREATE TABLE TMP (
  `qname`,
  `hostname`,
  `group`,
  `owner`,
  `jobname`,
  `jobnumber`,
  `account`,
  `priority`,
  `qsub_time`,
  `start_time`,
  `end_time`,
  `failed`,
  `exit_status`,
  `ru_wallclock`, -- = end_time - start_time
  `ru_utime`,
  `ru_stime`,
  `ru_maxrss`,
  `ru_ixrss`,
  `ru_ismrss`,
  `ru_idrss`,
  `ru_isrss`,
  `ru_minflt`,
  `ru_majflt`,
  `ru_nswap`,
  `ru_inblock`,
  `ru_oublock`,
  `ru_msgsnd`,
  `ru_msgrcv`,
  `ru_nsignals`,
  `ru_nvcsw`,
  `ru_nivcsw`,
  `project`,
  `department`,
  `granted_pe`,
  `slots`,
  `taskid`,
  `cpu`, -- = ru_utime + ru_stime
  `mem`,
  `io`,
  `category`,
  `iow`,
  `pe_taskid`,
  `maxvmem`,
  `arid`,
  `ar_submission_time`
);

.separator :
-- $SGE_ROOT/$SGE_CELL/common/accounting
.import accounting TMP

CREATE INDEX idx_TMP_owner ON TMP(`owner`);
CREATE INDEX idx_TMP_jobname ON TMP(`jobname`);
CREATE INDEX idx_TMP_jobid ON TMP(`jobnumber`);

DELETE FROM TMP WHERE owner != 'akuzniar';

CREATE VIEW V_JOB AS
SELECT
    CAST(`jobnumber` AS INTEGER) AS `jobid`,
    `jobname`,
    `hostname`,
    CAST(`slots` AS INTEGER) AS `slots`,
    strftime('%Y-%m-%d %H:%M:%S', `qsub_time`, 'unixepoch') AS `subtime`,
    strftime('%Y-%m-%d %H:%M:%S', `start_time`, 'unixepoch') AS `stime`,
    strftime('%Y-%m-%d %H:%M:%S', `end_time`, 'unixepoch') AS `etime`,
    CAST(`ru_wallclock` AS INTEGER) AS `runtime`, -- in sec
    CAST(`cpu` AS INTEGER) AS `cputime`, -- in sec
    strftime('%Y-%m-%d %H:00:00', `start_time`, 'unixepoch') AS `stime_bin`,
    strftime('%Y-%m-%d %H:00:00', `end_time`, 'unixepoch') AS `etime_bin`,
    (CAST(`maxvmem` AS INTEGER) / 1024  / 1024) AS `mem_mb`,
    (CASE WHEN CAST(`failed` AS INTEGER) = 0 AND CAST(`exit_status` AS INTEGER) = 0 THEN 'COMPLETED' ELSE 'FAILED' END) AS `status`
FROM TMP;

CREATE VIEW VV_JOB AS
SELECT
    `jobname`,
    COUNT(`jobid`) AS `n_jobs`,
    `status`,
    MIN(`runtime`) AS `min_runtime`,
    MAX(`runtime`) AS `max_runtime`,
    CAST(ROUND(AVG(`runtime`)) AS INTEGER) AS `mean_runtime`,
    MIN(`mem_mb`) AS `min_mem`,
    MAX(`mem_mb`) AS `max_mem`,
    CAST(ROUND(AVG(`mem_mb`)) AS INTEGER) AS `mean_mem`
FROM V_JOB
GROUP BY `jobname`, `status`;

CREATE VIEW V_JOB_E AS
SELECT
  `etime_bin`,
  `status`,
  COUNT(DISTINCT `hostname`) AS `n_hosts`,
  SUM(`slots`) AS `n_cores`,
  COUNT(*) AS `n_jobs`
FROM V_JOB
GROUP BY `etime_bin`, `status` ORDER BY 1;

CREATE VIEW VV_JOB_E AS
SELECT
  T1.etime_bin,
  T1.status,
  T1.n_jobs,
  SUM(T2.n_jobs) AS cum_jobs
FROM V_JOB_E AS T1, V_JOB_E AS T2
WHERE T2.rowid <= T1.rowid AND T2.status = T1.status
GROUP BY T1.status, T1.etime_bin, T1.n_jobs
ORDER BY T1.status, T1.rowid;

-- .header on
-- .mode csv
-- .output jobs.csv
-- SELECT etime_bin, cum_jobs FROM VV_JOB_E WHERE status = 'COMPLETED';
