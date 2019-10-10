suppressPackageStartupMessages(require(StructuralVariantAnnotation))

get_input_path <- function(sv_caller) {
    paste('benchmark/out/3/S3/', sv_caller, '_out/', sv_caller, '.vcf', sep='')
}

#SV type assignment based on
# https://github.com/PapenfussLab/gridss/blob/7b1fedfed32af9e03ed5c6863d368a821a4c699f/example/simple-event-annotation.R#L9
apply_svtype <- function(gr)
{
  gr$svtype <-
    ifelse(
      seqnames(gr) != seqnames(partner(gr)),
      "BP",
      ifelse(
        gr$insLen >= abs(gr$svLen) * 0.7,
        "INS",
        ifelse(
          strand(gr) == strand(partner(gr)),
          "INV",
          ifelse(xor(
            start(gr) < start(partner(gr)), strand(gr) == "-"
          ), "DEL",
          "DUP")
        )
      )
    )
  gr
}

load_sv_caller_vcf <- function(vcf_file, confidence_regions_gr, sample, sv_caller) {
    # vcf_file <- truth_set_file[[sample]]
    sv_callset_vcf <- VariantAnnotation::readVcf(vcf_file)

    if (sv_caller == 'lumpy') {
        # Read evidence support as a proxy for QUAL
      support <- unlist(info(sv_callset_vcf)$SU)
      fixed(sv_callset_vcf)$QUAL <- support
    } else if (sv_caller == 'delly') {
      # Split-read support plus Paired-end read support as a proxy for QUAL
      sr_support <- info(sv_callset_vcf)$SR
      sr_support[is.na(sr_support)] <- 0
      fixed(sv_callset_vcf)$QUAL <-
        sr_support + info(sv_callset_vcf)$PE
    }

    bpgr <- breakpointRanges(sv_callset_vcf)
    begr <- breakendRanges(sv_callset_vcf)
    gr <- sort(c(bpgr, begr))
    if (sv_caller %in% c('gridss', 'manta')) {
      gr <- apply_svtype(gr)
    }
    # Select DEL
    gr <- gr[which(gr$svtype == "DEL")]
    # Select DEL >= 50 bp
    gr <- gr[gr$svLen <= (-50)]
    gr
  }

load_truth_set_vcf <-
  function(vcf_file, sample)
  {
    sv_callset_vcf <-
      VariantAnnotation::readVcf(vcf_file)
    bpgr <- breakpointRanges(sv_callset_vcf)
    begr <- breakendRanges(sv_callset_vcf)
    gr <- sort(c(bpgr, begr))

    gr <- gr[which(gr$svtype == "DEL")]
    gr <- gr[gr$svLen <= (-50)]

    gr <- remove_blacklist(gr, confidence_regions_gr, sample)
    gr
  }

load_bedpe <- function(bedpe_file,
                       sample)
{
  sv_callset_bedpe <- rtracklayer::import(bedpe_file)
  bpgr <- pairs2breakpointgr(sv_callset_bedpe)
  gr <- sort(bpgr)
  #gr <- remove_blacklist(gr, confidence_regions_gr, sample)
  gr
}

make_percent <- function(x){
  round(x*100,digits = 1)
}

### main ###
test_sample <- c('NA12878')
truth_set_file <- list()
truth_set_file[[test_sample]] <-
  'Personalis_1000_Genomes_deduplicated_deletions.bed' # .bedpe
truth_set_file
gr <- list()
gr[[test_sample]] <- list()
truth_set <- list()
truth_set[[test_sample]] <- load_bedpe(truth_set_file[[test_sample]], test_sample)

#Specify output directory for results
outDir <- ""

sv_caller_list <- c('gridss', 'manta', 'lumpy', 'delly')

# Load sv_callers results
  for (sv_caller in sv_caller_list)
  {
    print(paste('Loading', sv_caller))

    vcf_file <- get_input_path(sv_caller)

    print(paste('Loading ',vcf_file,sep=''))
    gr[[test_sample]][[sv_caller]] <- load_sv_caller_vcf(vcf_file, test_sample, sv_caller)
  }

  # Only consider Chr1 to ChrX
  for (sv_caller in sv_caller_list)
  {
    print(sv_caller)
    print(length(gr[[test_sample]][[sv_caller]]))
    print(length(gr[[test_sample]][[sv_caller]][seqnames(gr[[test_sample]][[sv_caller]]) !=
                                             'Y']))
  }

  truth_svgr <- truth_set[[test_sample]]
  for (sv_caller in sv_caller_list)
  {
    gr[[test_sample]][[sv_caller]]$caller <- sv_caller
  }

  for (svcaller in sv_caller_list)
  {
    gr[[test_sample]][[svcaller]]$truth_matches <-
      countBreakpointOverlaps(
        gr[[test_sample]][[svcaller]],
        truth_svgr,
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

    svgr <- c(gr[[test_sample]][['gridss']],
              gr[[test_sample]][['manta']],
              gr[[test_sample]][['lumpy']],
              gr[[test_sample]][['delly']])

  res.df <- as.data.frame(svgr) %>%
    dplyr::select(caller, truth_matches) %>%
    dplyr::group_by(caller) %>%
    dplyr::summarise(calls = n(),
                     tp = sum(truth_matches > 0)) %>%
    dplyr::group_by(caller) %>%
    dplyr::mutate(
      cum_tp = cumsum(tp),
      cum_n = cumsum(calls),
      cum_fp = cum_n - cum_tp,
      precision = round(cum_tp / cum_n,digits = 1),
      recall = round(cum_tp / length(truth_svgr), digits = 1)
    )
  res.df$F1 = with(res.df, 2 * (precision * recall) / (precision + recall))

  res.df$precision <- make_percent(res.df$precision)
  res.df$recall <- make_percent(res.df$recall)
  res.df$F1 <- make_percent(res.df$F1)

  write.table(res.df, file=paste(outDir, 'performance_results.csv',sep=''), quote=F, row.names = F)

}
