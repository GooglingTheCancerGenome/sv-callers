suppressPackageStartupMessages(require(StructuralVariantAnnotation))

get_file <- function(sv_caller) {
    paste('benchmark/out/3/S3/', sv_caller, '_out/', sv_caller, '.vcf', sep='')
}

# SV type assignment based on
# https://github.com/PapenfussLab/gridss/blob/7b1fedfed32af9e03ed5c6863d368a821a4c699f/example/simple-event-annotation.R#L9
get_svtype <- function(gr) {
    return(ifelse(seqnames(gr) != seqnames(partner(gr)), "ITX", # LS: BP or CTX? (#74)
        ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS",
            ifelse(strand(gr) == strand(partner(gr)), "INV",
                ifelse(xor(start(gr) < start(partner(gr)), strand(gr) == "-"), "DEL", "DUP")))))
}

# convert (non-conformant) BED to BEDPE file
# LS: b37 reference used?
bed_file <- 'Personalis_1000_Genomes_deduplicated_deletions.bed'
bedpe_file <- 'Personalis_1000_Genomes_deduplicated_deletions.bedpe'
cmd <- paste('cat', bed_file , '|', "awk 'BEGIN {a=0; OFS=\"\t\"} NR>1 {print $1,$2,$2+1,$1,$3,$3+1, \"DEL_\" a,0,\"+\",\"+\",\"DEL\"; a+=1}'", '>', bedpe_file)
system(cmd)

# import truth set (DELs) from BEDPE file
bedpe <- rtracklayer::import(bedpe_file)
true_gr <- pairs2breakpointgr(bedpe)
print(length(true_gr))

# read SV callers' output in VCF
sv_callers <- c('manta', 'delly', 'lumpy', 'gridss')
for (c in sv_callers) {
    vcf_file <- get_file(c)
    vcf <- VariantAnnotation::readVcf(vcf_file, 'b37')
    #vcf
    gr <- breakpointRanges(vcf)
    print(paste(c, length(gr)))
    # select DELs only
    gr <- gr[gr$svtype == 'DEL']
    print(paste(c, length(gr))) # LS: 0 returned for gridss
}
q()
