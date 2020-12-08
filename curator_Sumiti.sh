#!/bin/bash -l


module load python3/3.6.5
module load umitools/1.0.0


#creating whitelist for barcodes

umi_tools whitelist --stdin SRR3879604_1_bc.fastq.gz \
--bc-pattern=CCCCCCCCCCCCCCCCCCCNNNNNN \
--knee-method=density \
--plot-prefix=knee_whitelist \
--log2stderr > whitelist_test_2_604.txt

umi_tools whitelist --stdin SRR3879605_1_bc.fastq.gz \
--bc-pattern=CCCCCCCCCCCCCCCCCCCNNNNNN \
--knee-method=density \
--plot-prefix=knee_whitelist \
--log2stderr > whitelist_test_605.txt

umi_tools whitelist --stdin SRR3879606_1_bc.fastq.gz \
--bc-pattern=CCCCCCCCCCCCCCCCCCCNNNNNN \
--knee-method=density \
--plot-prefix=knee_whitelist \
--log2stderr > whitelist_test_606.txt






#Extract the barcodes and filter the reads
umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCCCCNNNNNN \
--extract-method=regex\
--stdin SRR3879604_1.fastq.gz \
--stdout SRR3879604_1_extracted.fastq.gz \
--read2-in SRR3879604_2.fastq.gz \
--read2-out=SRR3879604_2_extracted.fastq.gz \
--filter-cell-barcode \
--whitelist=whitelist_test_2_604.txt

umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCCCCNNNNNN \
--extract-method=regex\
--stdin SRR3879605_1.fastq.gz \
--stdout SRR3879605_1_extracted.fastq.gz \
--read2-in SRR3879605_2.fastq.gz \
--read2-out=SRR3879605_2_extracted.fastq.gz \
--filter-cell-barcode \
--whitelist=whitelist_test_605.txt

umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCCCCNNNNNN \
--extract-method=regex\
--stdin SRR3879606_1.fastq.gz \
--stdout SRR3879606_1_extracted.fastq.gz \
--read2-in SRR3879606_2.fastq.gz \
--read2-out=SRR3879606_2_extracted.fastq.gz \
--filter-cell-barcode \
--whitelist=whitelist_test_606.txt

#salmon index


grep "^>" <(gunzip -c GRCh38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys.txt

cat gencode.v33.transcripts.fa.gz GRCh38.primary_assembly.genome.fa.gz > gentrome.fa.gz

salmon index -t gentrome.fa.gz -d decoys.txt -p 12 -i salmon_index --gencode


#transcript to gene map file


awk -F '|' '/ENST/ {print $1 "\t" $2}' <(gunzip -c GRCh38.primary_assembly.genome.fa.gz) > txp2gene_test.tsv
bioawk -c gff '$feature=="transcript" {print $group}' <(gunzip -c gencode.v31.primary_assembly.annotation.gtf.gz) | awk -F ' ' '{print substr($4,2,length($4)-3) "\t" substr($2,2,length($2)-3)}' - > txp2gene_test.tsv

salmon index -i index -k 31 --gencode -p 4 -t gencode.v33.transcripts.fa.gz

salmon alevin -lISR -1 SRR3879605_1.fastq.gz -2 SRR3879606_2.fastq.gz -i indexes -p 8 -o alevin_output --tgMap txp2gene.tsv --end 5 --barcodeLength 19 --umiLength 6 --whitelist SRR3879604_whitelist.csv

salmon alevin -lISR -1 SRR3879604_1_bc.fastq.gz -2 SRR3879604_2_extracted.fastq.gz -i indexes -p 16 --whitelist wl_umi_SRR3879604.tsv --dumpFeatures -o alevin_output_extracted_2 --tgMap txp2gene.tsv --end 5 --barcodeLength 19 --umiLength 6 --dumpMtx


