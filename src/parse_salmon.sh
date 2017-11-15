#!/usr/bin/sh

mkdir -p data/salmon_output/tmp/length/
mkdir -p data/salmon_output/tmp/effective_length/
mkdir -p data/salmon_output/tmp/TPM/
mkdir -p data/salmon_output/tmp/counts/

cut -f 1 data/salmon_output/single_end/D15_P1373_Run1/Run00261_L5_YFVNM_D15_G11_P1373_1087_transcripts_quant/quant.sf > data/salmon_output/tmp/genes.csv

find data/salmon_output/ -name "quant.sf" | \
perl -pe 's/(.*)/\1 \1/g' | \
perl -pe 's/(.*) [^\s]*_(P\d{4}_[\d]+)_.*\/quant.sf/\1 \2/g' | \
awk '{system("sort -k 1 "$1" | cut -f 2 > data/salmon_output/tmp/length/"$2"_length.csv")}'

find data/salmon_output/ -name "quant.sf" | \
perl -pe 's/(.*)/\1 \1/g' | \
perl -pe 's/(.*) [^\s]*_(P\d{4}_[\d]+)_.*\/quant.sf/\1 \2/g' | \
awk '{system("sort -k 1 "$1" | cut -f 3 "$1" > data/salmon_output/tmp/effective_length/"$2"_effective_length.csv")}'

find data/salmon_output/ -name "quant.sf" | \
perl -pe 's/(.*)/\1 \1/g' | \
perl -pe 's/(.*) [^\s]*_(P\d{4}_[\d]+)_.*\/quant.sf/\1 \2/g' | \
awk '{system("sort -k 1 "$1" | cut -f 4 "$1" > data/salmon_output/tmp/TPM/"$2"_TPM.csv")}'

find data/salmon_output/ -name "quant.sf" | \
perl -pe 's/(.*)/\1 \1/g' | \
perl -pe 's/(.*) [^\s]*_(P\d{4}_[\d]+)_.*\/quant.sf/\1 \2/g' | \
awk '{system("sort -k 1 "$1" | cut -f 4 "$1" > data/salmon_output/tmp/"$2"_counts.csv")}'

find data/salmon_output/tmp/length/ -name "*.csv" | perl -pe 's/(.*(P\d{4}_[\d]+).*\.csv)/\1 \2/g' | awk '{system("echo "$2"; cat "$1)}' > data/salmon_output/tmp/length.csv

awk '
{
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' data/salmon_output/tmp/length.csv


less data/salmon_output/single_end/D15_P1373_Run1/Run00261_L5_YFVNM_D15_G11_P1373_1087_transcripts_quant/quant.sf
