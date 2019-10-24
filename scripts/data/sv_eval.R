#
# This script evaluates SV callsets obtained by Manta, DELLY, LUMPY and/or GRIDSS
# given a truth set (DELs only):
#  1. Personalis + 1k Genomes Project (Parikh et al., 2016)
#  2. PacBio/Moleculo (Layer et al., 2014)
#

suppressPackageStartupMessages(require(tools))
suppressPackageStartupMessages(require(StructuralVariantAnnotation))

min.supp <- 3  # min. number of callers supporting an SV
data.dir <- file.path('benchmark', 'out', '3', 'S3')  # with SV callsets

# "truth" sets with BED/BEDPE files
truth.sets <- list()
truth.sets$Personalis1kGP <- list(
    sv.file='Personalis_1000_Genomes_deduplicated_deletions.bed',
    excl.file='ENCFF001TDO.bed')
truth.sets$PacbioMoleculo <- list(
    sv.file='3717462611446476_add4.bedpe',
    excl.file='ceph18.b37.lumpy.exclude.2014-01-15.bed')
truth.sets <- lapply(truth.sets, lapply,
                     function(fn){file.path('benchmark', 'in', fn)})
sel.idx <- 1  # select either 1=Personalis1kGP or 2=PacbioMoleculo

# get VCF file path given a caller
getVcf <- function(caller) {
    return(file.path(data.dir, paste0(caller, '_out'), paste0(caller, '.vcf')))
}

# assign SV types
# https://github.com/PapenfussLab/gridss/blob/7b1fedfed32af9e03ed5c6863d368a821a4c699f/example/simple-event-annotation.R#L9
getSvType <- function(gr) {
    return(ifelse(seqnames(gr) != seqnames(partner(gr)), 'CTX',
        ifelse(gr$insLen >= abs(gr$svLen) * 0.7, 'INS',
            ifelse(strand(gr) == strand(partner(gr)), 'INV',
                ifelse(xor(start(gr) < start(partner(gr)), strand(gr) == '-'),
                    'DEL', 'DUP')))))
}

# compute performance metrics for a callset
getPerfMetrics <- function(callset, hits, n.true) {
    n <- length(hits)
    tp <- sum(hits)
    fp <- n - tp
    fn <- n.true - tp
    prec <- round(tp * 100 / n, digits=1)
    rec <- round(tp * 100 / n.true, digits=1)
    return(list(callset=callset, n=n, tp=tp, fp=fp, fn=fn, precision=prec,
                recall=rec))
}

# exclude genomic regions
excludeRegions <- function(query.gr, subject.gr) {
    return(query.gr[!(overlapsAny(query.gr, subject.gr) |
                      overlapsAny(partner(query.gr), subject.gr)), ])
}

# count breakpoint overlaps (hits)
getHits <- function(query.gr, subject.gr) {
    return(countBreakpointOverlaps(query.gr, subject.gr, maxgap=200,
                                   sizemargin=0.25, ignore.strand=TRUE,
                                   restrictMarginToSizeMultiple=0.5,
                                   countOnlyBest=TRUE))
}

# convert BED file to BEDPE if applicable
bed.file <- truth.sets$Personalis1kGP$sv.file
bedpe.file <- paste0(file_path_sans_ext(bed.file), '.bedpe')
if(!file.exists(bedpe.file)) {
    cmd <- paste("awk 'BEGIN {a=0; OFS=\"\t\"} NR>1 {print $1,$2,$2+1,$1,$3,\
                 $3+1,\"DEL_\" a,-1,\"+\",\"+\",\"DEL\"; a+=1}'", bed.file, '>',
                 bedpe.file)
    system(cmd)
    message("DEBUG:", cmd)
} else {
    truth.sets$Personalis1kGP$sv.file <- bedpe.file
}

# import deletions from BEDPE
bedpe.file <- truth.sets[[sel.idx]]$sv.file
message("DEBUG:", bedpe.file)
true.gr <- pairs2breakpointgr(rtracklayer::import(bedpe.file))
seqlevelsStyle(true.gr) <- 'NCBI'  # ensure chr[X] -> [X]
min.svLen <- min(abs(end(partner(true.gr)) - start(true.gr)) + 1)
message('### Truth set ###')
message('input = ', bedpe.file)
message('n = ', length(true.gr))
message('min.svLen = ', min.svLen)

# import the ENCODE's exclusion list
bed.excl.file <- truth.sets[[sel.idx]]$excl.file
excl.gr <- rtracklayer::import(bed.excl.file)
seqlevelsStyle(excl.gr) <- 'NCBI'  # chr[X] -> X
message('\n### Exclusion list ###')
message('input = ', bed.excl.file)
print(seqnames(excl.gr))

true.gr <- excludeRegions(true.gr, excl.gr)
min.svlen <- min(abs(end(partner(true.gr)) - start(true.gr)) + 1)
n.true <- length(true.gr)
message('\n### Truth set filtered by the exclusion list ###')
message('n = ', n.true)
message('min.svLen = ', min.svLen)

# import SV callsets from VCF files
callers <- c('manta', 'delly', 'lumpy', 'gridss')
hits.df <- data.frame(callset=character(), n=numeric(), tp=numeric(),
                      fp=numeric(), precision=numeric(), recall=numeric())
for (c in callers) {
  vcf.file <- getVcf(c)
  vcf <- VariantAnnotation::readVcf(vcf.file)
  # select only DELs
  gr <- breakpointRanges(vcf)
  gr$svtype <- getSvType(gr)
  gr <- gr[gr$svtype == 'DEL']
  message('\n### All DELs ###')
  message('# ', c)
  message('n = ', length(gr))
  print(summary(gr$svLen))

  message('\n### DELs svLen != NA ###')
  gr <- gr[!is.na(gr$svLen)]
  message('# ', c)
  message('n = ', length(gr))
  print(summary(gr$svLen))

  gr <- gr[abs(gr$svLen) >= min.svLen]
  message('\n### DELs svLen >= min.svLen ###')
  message('# ', c)
  message('n = ', length(gr))
  print(summary(gr$svLen))

  gr <- excludeRegions(gr, excl.gr)
  message('\n### DELs >= min.svlen AND filtered by the exclusion list ###')
  message('# ', c)
  message('n = ', length(gr))
  print(summary(gr$svLen))

  hits <- getHits(gr, true.gr)
  pm <- getPerfMetrics(c, hits, n.true)
  hits.df <- rbind(hits.df, data.frame(pm))
}

# import SURVIVOR merge callset from VCF file
vcf.infile <- file.path(data.dir, 'all.vcf')
vcf.outfile <- 'merge.vcf'  # fixed file: replace ':' by '_' in ID & SAMPLE fields
cmd <- paste("awk '{if ($1 ~ /^#/){print} else {id=$3; gsub(\":\",\"_\",$3);\
             gsub(id,$3,$10); print}}'", vcf.infile, '>', vcf.outfile)
system(cmd)

vcf <- VariantAnnotation::readVcf(vcf.outfile)
info(header(vcf))$Type[1:2] <- c("Integer", "Integer")  # fix INFO/CI{END,POS}
                                                        # types: String->Integer
vcf <- vcf[which(info(vcf)$SVTYPE == 'DEL')]  # keep only deletions
vcf <- vcf[which(as.integer(info(vcf)$SUPP) >= min.supp, TRUE)]  # filter calls
                                                                 # by support
#vcf <- vcf[which(info(vcf)$SUPP_VEC == '1101', TRUE)]
info(vcf)$CIPOS <- IntegerList(info(vcf)$CIPOS)  # fix type:
info(vcf)$CIEND <- IntegerList(info(vcf)$CIEND)  # CharacterList->IntegerList
gr <- breakpointRanges(vcf)
gr <- excludeRegions(gr, excl.gr)
hits <- getHits(gr, true.gr)
callset <- file_path_sans_ext(vcf.outfile)
pm <- getPerfMetrics(callset, hits, n.true)
hits.df <- rbind(hits.df, data.frame(pm))
message('\n### Performance metrics ###')
message('min.supp = ', min.supp, '\n')
hits.df
csv.file <- 'sv_eval.csv'
write.table(hits.df, file=csv.file, row.names=FALSE, col.names=TRUE,
            quote=FALSE, sep=',')
q()
