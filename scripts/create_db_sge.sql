
CREATE TABLE TMP (`qname`,`hostname`,`group`,`owner`,`jobname`,`jobnumber`,`account`,`priority`,`qsub_time`,`start_time`,`end_time`,`failed`,`exit_status`,`ru_wallclock`,`ru_utime`,`ru_stime`,`ru_maxrss`,`ru_ixrss`,`ru_ismrss`,`ru_idrss`,`ru_isrss`,`ru_minflt`,`ru_majflt`,`ru_nswap`,`ru_inblock`,`ru_oublock`,`ru_msgsnd`,`ru_msgrcv`,`ru_nsignals`,`ru_nvcsw`,`ru_nivcsw`,`project`,`department`,`granted_pe`,`slots`,`taskid`,`cpu`,`mem`,`io`,`category`,`iow`,`pe_taskid`,`maxvmem`,`arid`,`ar_submission_time`);

.separator :
.import accounting TMP -- $SGE_ROOT/$SGE_CELL/common/accounting

CREATE INDEX idx_TMP_owner ON TMP(owner);
CREATE INDEX idx_TMP_jobname ON TMP(jobname);
CREATE INDEX idx_TMP_jobid ON TMP(jobnumber);

DELETE FROM TMP WHERE owner != 'akuzniar';

CREATE VIEW V_JOB AS
SELECT
    CAST(jobnumber AS INTEGER) AS jobid,
    jobname,
    strftime('%Y-%m-%d %H:%M:%S', start_time, 'unixepoch') AS start,
    strftime('%Y-%m-%d %H:%M:%S', end_time, 'unixepoch') AS end,
    CAST(ru_wallclock AS INTEGER) AS rtime_s,
    strftime('%Y-%m-%d %H:00:00', start_time, 'unixepoch') AS start_bin,
    strftime('%Y-%m-%d %H:00:00', end_time, 'unixepoch') AS end_bin,
    (CAST(maxvmem AS INTEGER) / 1024  / 1024) AS mem_mb,
    (CASE WHEN CAST(failed AS INTEGER) = 0 AND CAST(exit_status AS INTEGER) = 0 THEN 'COMPLETED' ELSE 'FAILED' END) AS state
FROM TMP;

CREATE VIEW VV_JOB AS
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
FROM V_JOB
GROUP BY jobname, state
