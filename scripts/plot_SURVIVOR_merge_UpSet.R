#@author Luca Santuari
#@license Apache License, Version 2.0

require(ggplot2)
require(vcfR)
require(UpSetR)

int_folders <- "patients/Patient1"
#int_folders <- "trio/NA12878"

#Plot germline SVs?
#work_dir <-
#	paste("/Users/lsantuari/Documents/Data/germline/",int_folders,"/SV/Filtered/", sep="")
#out_dir <-
#	paste(work_dir, "SURVIVOR/", sep="")

work_dir <- getwd()
#output in working directory
out_dir <- work_dir

#Use prefix
#prefix <- "manta_10runs"
#prefix <- "run1"
#vcf_filepath <-
#	paste(work_dir, prefix, '_survivor_merge.vcf', sep = '')

vcf_filepath <-
paste(work_dir, 'survivor_merge.vcf', sep = '')
vcf_rda_filepath <- paste(vcf_filepath, '.rda', sep = '')

vcf <- read.vcfR(vcf_filepath, verbose = FALSE)
save(file = vcf_rda_filepath, vcf, compress = T)

#Extract INFO->SUPP_VEC, separate each 0/1
mysplit_list <-
strsplit(sapply(strsplit(sapply(strsplit(getFIX(vcf, getINFO = T)[, 'INFO'], ';'), function(x) {
    x[2]
}), '='), function(x) {
    x[2]
}), '')

#Create a 0/1 matrix
mat <-
matrix(0,
ncol = length(mysplit_list[[1]]),
nrow = length(mysplit_list))
for (i in 1 : nrow(mat))
{
    mat[i,] <- as.integer(mysplit_list[[i]])
}
#Get SVTYPE vector
svtype_vec <- getFIX(vcf)[, 'ALT']

#Name the columns, use proper pattern to extract the sample names
#pattern <- '/|_'
pattern <- '\\.'
colnames(mat) <-
sapply(strsplit(colnames(vcf@gt)[- 1], pattern), function(x) {
    x[[1]]
})

#Use the specified names?
#colnames(mat) <- c('Delly', 'GRIDSS', 'last_nanosv', 'Lumpy', 'Manta')
#head(mat)

if (prefix == 'run1')
{
    colnames(mat) <-
    c('Delly', 'GRIDSS', 'Lumpy', 'Manta', 'last_nanosv')
}

sv <- data.frame(ID = getFIX(vcf)[, 'ID'], mat, SVTYPE = svtype_vec)

sv.freq <- as.data.frame(table(apply(sv[which(sv$SVTYPE == "<INS>"),
2 : (ncol(sv) - 1)], 1, paste, collapse = '')))
head(sv.freq[order(sv.freq$Freq, decreasing = T),])
sum(sv.freq$Freq)

#How many SVs of a type are there?
#unique(sv$SVTYPE)
#table(sv$SVTYPE)

#Plot SVTYPES
for (svname in c("DEL", "INS", "DUP", "INV", "TRA"))
{
    idx <- which(sv$SVTYPE == paste("<", svname, ">", sep = ""))
    if (length(idx) > 0)
    {
        print(svname)
        jpeg(
        filename = paste(out_dir,
        #prefix, '_',
        svname,
        ".jpeg",
        sep = ""),
        width = 2000,
        height = 1000
        )
        upset(
        sv[idx,],
        sets = names(sv)[- c(1, ncol(sv))],
        sets.bar.color = "#56B4E9",
        order.by = "freq",
        #empty.intersections = "on",
        text.scale = 3,
        set_size.angles = 90,
        number.angles = 0
        )
        dev.off()
    }
}
