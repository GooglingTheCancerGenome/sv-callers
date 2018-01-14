rule delly:
    threads:
        get_nthreads("delly")
    resources:
        mem_mb=get_maxmem("delly")

include: "delly-BND.smk"
include: "delly-DEL.smk"
include: "delly-DUP.smk"
include: "delly-INS.smk"
include: "delly-INV.smk"
