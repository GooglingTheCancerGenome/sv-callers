suppressPackageStartupMessages(require(StructuralVariantAnnotation))
suppressPackageStartupMessages(require(dplyr))

get_file <- function(sv_caller) {
  paste('benchmark/out/3/S3/',
        sv_caller,
        '_out/',
        sv_caller,
        '.vcf',
        sep = '')
}

# SV type assignment based on
# https://github.com/PapenfussLab/gridss/blob/7b1fedfed32af9e03ed5c6863d368a821a4c699f/example/simple-event-annotation.R#L9
get_svtype <- function(gr) {
  return(ifelse(
    seqnames(gr) != seqnames(partner(gr)),
    "ITX",
    # LS: BP or CTX? (#74)
    # AK: These are inter-chromosomal translocations (usually referred as CTX). However, in the sv_benchmark they are referred as 'BP'.
    ifelse(
      gr$insLen >= abs(gr$svLen) * 0.7,
      "INS",
      ifelse(
        strand(gr) == strand(partner(gr)),
        "INV",
        ifelse(xor(
          start(gr) < start(partner(gr)), strand(gr) == "-"
        ), "DEL", "DUP")
      )
    )
  ))
}

get_overlap <- function(gr1, gr2){
  
  countBreakpointOverlaps(
    gr1,
    gr2,
    # read pair based callers make imprecise calls.
    # A margin around the call position is required when matching with the truth set
    maxgap = 200,
    # Since we added a maxgap, we also need to restrict the mismatch between the
    # size of the events. We don't want to match a 100bp deletion with a
    # 5bp duplication. This will happen if we have a 100bp margin but don't also
    # require an approximate size match as well
    sizemargin = 0.25,
    ignore.strand = TRUE,
    # We also don't want to match a 20bp deletion with a 20bp deletion 80bp away
    # by restricting the margin based on the size of the event, we can make sure
    # that simple events actually do overlap
    restrictMarginToSizeMultiple =
      0.5,
    # Some callers make duplicate calls and will sometimes report a variant multiple
    # times with slightly different bounds. countOnlyBest prevents these being
    # double-counted as multiple true positives.
    countOnlyBest = TRUE
  )
  
}

# convert (non-conformant) BED to BEDPE file
# LS: b37 reference used?
# AK: the reference b37 is based on the reference version GRCh37 and so is the file Personalis_1000_Genomes_deduplicated_deletions.bed
bed_file <- 'Personalis_1000_Genomes_deduplicated_deletions.bed'
bedpe_file <- 'Personalis_1000_Genomes_deduplicated_deletions.bedpe'
cmd <-
  paste(
    'cat',
    bed_file ,
    '|',
    "awk 'BEGIN {a=0; OFS=\"\t\"} NR>1 {print $1,$2,$2+1,$1,$3,$3+1, \"DEL_\" a,0,\"+\",\"+\",\"DEL\"; a+=1}'",
    '>',
    bedpe_file
  )
system(cmd)

# import truth set (DELs) from BEDPE file
bedpe <- rtracklayer::import(bedpe_file)
true_gr <- pairs2breakpointgr(bedpe)
print(paste('bedpe', length(true_gr)))

# read SV callers' output in VCF
gr <- list()
sv_callers <- c('manta', 'delly', 'lumpy', 'gridss')
for (c in sv_callers) {
  vcf_file <- get_file(c)
  vcf <- VariantAnnotation::readVcf(vcf_file, 'b37')
  #vcf
  gr[[c]] <- breakpointRanges(vcf)
  print(paste(c, length(gr[[c]])))
  # select DELs only
  gr[[c]] <- gr[[c]][get_svtype(gr[[c]]) %in% c('DEL')]
  # select only DELs greater than 50 bp
  gr[[c]] <- gr[[c]][gr[[c]]$svLen <= (-50)]
  print(paste(c, 'DEL', length(gr[[c]])))
}

for (c in sv_callers) {
  gr[[c]]$caller <- c
}

#Adapted from the StructuralVariantAnnotation vignette:
#https://bioconductor.org/packages/devel/bioc/vignettes/StructuralVariantAnnotation/inst/doc/vignettes.html
for (c in sv_callers)
{
  gr[[c]]$truth_matches <- get_overlap(gr[[c]], true_gr)
}


#Generate all combinations of a binary vector of length 4
mat_comb <- expand.grid(replicate(4, 0:1, simplify = FALSE))[-1,]
colnames(mat_comb) <- sv_callers

gr_2way <- list()
#Generate 2-way overlap
for (c1 in sv_callers)
{
  for (c2 in sv_callers)
  {
    if (c1 != c2)
    {
      caller_name <- paste(c1, c2, sep = '_')
      gr_2way[[caller_name]] <- gr[[c1]]
      
      #Rename SV entries
      names(gr_2way[[caller_name]]) <- paste(caller_name,1:length(gr_2way[[caller_name]]),sep="_")
      
      gr_2way[[caller_name]]$caller <- caller_name
      gr_2way[[caller_name]]$truth_matches <- get_overlap(gr[[c1]], gr[[c2]])
      
    }
  }
}

gr_3way <- list()
#Generate 3-way overlap
for (c1 in names(gr_2way))
{
  for (c2 in sv_callers[-which(sv_callers%in%strsplit(c1, "_")[[1]])])
  {
    caller_name <- paste(c1, c2, sep = '_')
    gr_3way[[caller_name]] <- gr_2way[[c1]]
    
    #Rename SV entries
    names(gr_3way[[caller_name]]) <- paste(caller_name,1:length(gr_3way[[caller_name]]),sep="_")
    
    gr_3way[[caller_name]]$caller <- caller_name
    gr_3way[[caller_name]]$truth_matches <- get_overlap(gr_2way[[c1]], gr[[c2]])
      
  }
}

gr_4way <- list()
#Generate 3-way overlap
for (c1 in names(gr_3way))
{
  for (c2 in sv_callers[-which(sv_callers%in%strsplit(c1, "_")[[1]])])
  {
    caller_name <- paste(c1, c2, sep = '_')
    gr_4way[[caller_name]] <- gr_3way[[c1]]
    
    #Rename SV entries
    names(gr_4way[[caller_name]]) <- paste(caller_name,1:length(gr_4way[[caller_name]]),sep="_")
    
    gr_4way[[caller_name]]$caller <- caller_name
    gr_4way[[caller_name]]$truth_matches <- get_overlap(gr_3way[[c1]], gr[[c2]])
  }
}

#Generate full data.frame
svgr <- c(gr[['manta']],
          gr[['delly']],
          gr[['lumpy']],
          gr[['gridss']])

for(c in names(gr_2way))
{
  svgr <- c(svgr, gr_2way[[c]])
}
for(c in names(gr_3way))
{
  svgr <- c(svgr, gr_3way[[c]])
}
for(c in names(gr_4way))
{
  svgr <- c(svgr, gr_4way[[c]])
}
print(length(svgr))

#Calculate Precision, Recall and F1-score for the overlaps
res.df <- as.data.frame(svgr) %>%
  dplyr::select(caller, truth_matches) %>%
  dplyr::group_by(caller) %>%
  dplyr::summarise(calls = n(),
                   tp = sum(truth_matches > 0)) %>%
  dplyr::group_by(caller) %>%
  dplyr::mutate(
    n = cumsum(calls),
    fp = n - tp,
    precision = round(tp / n, digits = 1),
    recall = round(tp / length(true_gr), digits = 1)
  )
res.df$F1 = with(res.df, 2 * (precision * recall) / (precision + recall))

res.df$precision <- make_percent(res.df$precision)
res.df$recall <- make_percent(res.df$recall)
res.df$F1 <- make_percent(res.df$F1)

print(res.df)

q()
