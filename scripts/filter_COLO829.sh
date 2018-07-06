#!/bin/bash -x

#This bash script fixes Lumpy VCF header and filters the VCF files of the COLO829 runs by quality

BASEDIR="COLO829"

CMD="filter_VCF.py"

#Lumpy VCF do not contain contig information. This loop adds contig info to the header of Lumpy VCF file
for RUN in `seq 1 10`; do
    for SV_CALLER in Lumpy; do
        echo "Fixing header for "$SV_CALLER": run "$RUN
        SV_CALLER_LC=`echo "$SV_CALLER" | awk '{print tolower($0)}'`

        VCF_IN=$BASEDIR$RUN"/T"$RUN"--N"$RUN"/"$SV_CALLER_LC"_out/"$SV_CALLER_LC".vcf"
        cp $VCF_IN $VCF_IN".nocontigs"

        awk '/^#CHROM/ { printf("##contig=<ID=1,length=249250621>\n##contig=<ID=2,length=243199373>\n##contig=<ID=3,length=198022430>\n##contig=<ID=4,length=191154276>\n##contig=<ID=5,length=180915260>\n##contig=<ID=6,length=171115067>\n##contig=<ID=7,length=159138663>\n##contig=<ID=8,length=146364022>\n##contig=<ID=9,length=141213431>\n##contig=<ID=10,length=135534747>\n##contig=<ID=11,length=135006516>\n##contig=<ID=12,length=133851895>\n##contig=<ID=13,length=115169878>\n##contig=<ID=14,length=107349540>\n##contig=<ID=15,length=102531392>\n##contig=<ID=16,length=90354753>\n##contig=<ID=17,length=81195210>\n##contig=<ID=18,length=78077248>\n##contig=<ID=19,length=59128983>\n##contig=<ID=20,length=63025520>\n##contig=<ID=21,length=48129895>\n##contig=<ID=22,length=51304566>\n##contig=<ID=X,length=155270560>\n##contig=<ID=Y,length=59373566>\n##contig=<ID=MT,length=16569>\n");} {print;}' $VCF_IN".nocontigs" > $VCF_IN
    done
done

#Filter VCF files from the 10 runs
for RUN in `seq 1 10`; do
    for SV_CALLER in Delly Lumpy Manta GRIDSS; do
        echo "Running filtering for "$SV_CALLER": run "$RUN
        SV_CALLER_LC=`echo "$SV_CALLER" | awk '{print tolower($0)}'`

        OUT_DIR=$BASEDIR"../processed/"$RUN"/T"$RUN"--N"$RUN"/"$SV_CALLER_LC"_out/"
        [ ! -d "$OUT_DIR" ] && mkdir -p "$OUT_DIR"

        VCF_IN=$BASEDIR$RUN"/T"$RUN"--N"$RUN"/"$SV_CALLER_LC"_out/"$SV_CALLER_LC".vcf"
        VCF_OUT=$OUT_DIR$SV_CALLER_LC".filtered.vcf"

        echo "python "$CMD" -i "$VCF_IN" -o "$VCF_OUT" -c "$SV_CALLER | sh
    done
done
