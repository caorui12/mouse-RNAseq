**Prepare work directory and QC**

```
cd /s3_d4/caorui/mouse_RNAseq/
mkdir second_stream_data/rawdata
cd second_stream_data/rawdata
ln -s ../../X201SC20123423-Z01-F004_01/raw_data/* .
ln -s ../../X201SC20123423-Z01-F004_02/raw_data/* .
ls > sampleList.txt ## prepare sampleList

for sample in `cat sampleList2.txt`
do
echo $sample
~/fqc.pl  adp_qual -p -f "$sample"/"$sample"_1.fq -r "$sample"/"$sample"_2.fq -o "$sample" 
done

```
**Build index**
```
nohup rsem-prepare-reference --gtf /s3_d4/caorui/mouse_RNAseq/genome/Mus_musculus.GRCm39.104.gtf /s3_d4/caorui/mouse_RNAseq/genome/Mus_musculus.GRCm39.dna.toplevel.fa --bowtie2 mouse_reference -p 40
```
**quantify
```
cd ..
mkdir quantify && cd quantify

for sample in `cat ../rawdata/sampleList.txt`
do      
echo $sample
rsem-calculate-expression --paired-end -no-bam-output --bowtie2 --append-names -p 20 \
../rawdata/"$sample"_1.fastq \
../rawdata/"$sample"_2.fastq \
../rawdata/mouse_reference \
$sample 
done
