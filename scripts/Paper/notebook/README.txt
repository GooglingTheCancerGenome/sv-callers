=== NA12878 ===

# filter after merge

SURVIVOR merge input.txt 100 1 0 0 0 0 all.vcf
merging entries: 15081
merging entries: 11797
merging entries: 334
merging entries: 116891

SURVIVOR filter all.vcf ENCFF001TDO.bed -1 -1 0 -1 all.filtered.vcf
SVs ignored: 3753

all.filtered.vcf:113202
all.vcf:116955
delly.vcf:15081
gridss.vcf:11797
lumpy.vcf:334
manta.vcf:116891


# filter before merge

for f in $(ls *.vcf); do n=$(basename $f .vcf); SURVIVOR filter $f ENCFF001TDO.bed -1 -1 0 -1 $n.filtered.vcf; done
SVs ignored: 6
SVs ignored: 0
SVs ignored: 0
SVs ignored: 972

SURVIVOR merge input.filtered.txt 100 1 0 0 0 0 all.filtered.vcf
merging entries: 15075
merging entries: 11797
merging entries: 334
merging entries: 115919

grep -cv "#" *.vcf
all.filtered.vcf:116036
delly.filtered.vcf:15075
delly.vcf:15081
gridss.filtered.vcf:11797
gridss.vcf:11797
lumpy.filtered.vcf:334
lumpy.vcf:334
manta.filtered.vcf:115919
manta.vcf:116891

=== COLO829 ===

# filter after merge

SURVIVOR merge input.txt 100 1 0 0 0 0 all.vcf
merging entries: 31
merging entries: 164
merging entries: 1587
merging entries: 142

SURVIVOR filter all.vcf ENCFF001TDO.bed -1 -1 0 -1 all.filtered.vcf
SVs ignored: 135

grep -cv "#" *.vcf
all.filtered.vcf:1609
all.vcf:1744
delly.vcf:31
gridss.vcf:164
lumpy.vcf:1587
manta.vcf:142


# filter before merge

for f in $(ls *.vcf); do n=$(basename $f .vcf); SURVIVOR filter $f ENCFF001TDO.bed -1 -1 0 -1 $n.filtered.vcf; done
SVs ignored: 0
SVs ignored: 0
SVs ignored: 0
SVs ignored: 0

SURVIVOR merge input.filtered.txt 100 1 0 0 0 0 all.filtered.vcf
merging entries: 31
merging entries: 164
merging entries: 1587
merging entries: 142

grep -cv "#" *.vcf
all.filtered.vcf:1744
delly.filtered.vcf:31
delly.vcf:31
gridss.filtered.vcf:164
gridss.vcf:164
lumpy.filtered.vcf:1587
lumpy.vcf:1587
manta.filtered.vcf:142
manta.vcf:142

