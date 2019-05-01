# Variant Calling for RNA seq Data
## Installing STAR for alignning RNA reads

```javascript
conda install -c bioconda star 
conda install -c bioconda/label/cf201901 star 
```

## Generating genome index
## 1- STAR uses genome index files that must be saved in unique directories. (only ch22 will be used not the whole genome)

``` javascript
genomeDir=/home/ngs-01/workdir/Assignment2/chr22_with_ERCC92
mkdir $genomeDir
STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles /home/ngs-01/workdir/sample_data/chr22_with_ERCC92.fa  --runThreadN 1 --sjdbGTFfile /home/ngs-01/workdir/sample_data/chr22_with_ERCC92.gtf
```

## 2- Alignment jobs were executed as follows:

``` javascript
runDir=/home/ngs-01/workdir/Assignment2/1pass
mkdir $runDir
cd $runDir
STAR --genomeDir $genomeDir --readFilesIn /home/ngs-01/workdir/Assignment2/ngs2-assignment-data/SRR8797509_1.part_001.part_001.fastq /home/ngs-01/workdir/Assignment2/ngs2-assignment-data/SRR8797509_2.part_001.part_001.fastq --runThreadN 4
```

## 3- For the 2-pass STAR, a new index is then created using splice junction information contained in the file SJ.out.tab from the first pass:

```javascript
genomeDir=/home/ngs-01/workdir/Assignment2/chr22_with_ERCC92_2pass
mkdir $genomeDir
STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles /home/ngs-01/workdir/sample_data/chr22_with_ERCC92.fa    --sjdbFileChrStartEnd /home/ngs-01/workdir/Assignment2/1pass/SJ.out.tab --sjdbOverhang 75 --runThreadN 4
```


## 4- The resulting index is then used to produce the final alignments as follows:

```javascript 
runDir=/home/ngs-01/workdir/Assignment2/2pass
mkdir $runDir
cd $runDir
STAR --genomeDir $genomeDir --readFilesIn /home/ngs-01/workdir/Assignment2/ngs2-assignment-data/SRR8797509_1.part_001.part_001.fastq /home/ngs-01/workdir/Assignment2/ngs2-assignment-data/SRR8797509_2.part_001.part_001.fastq --runThreadN 4
```

##  Add read groups, sort, mark duplicates, and create index

```javascript
R1=/home/ngs-01/workdir/Assignment2/SRR8797509_1.part_001.part_001.fastq
R2=/home/ngs-01/workdir/Assignment2/SRR8797509_2.part_001.part_001.fastq

SM="SRR8797509"		               
LB="SRR8797509_2019329"				
PL="Illumina"		               
PU="HiSeqXTen"	                   
RGID=$(cat /home/ngs-01/workdir/Assignment2/ngs2-assignment-data/SRR8797509_*.part_001.part_001.fastq | head -n1 | sed 's/ /_/g' | cut -d "_" -f1)              

picard_path="/home/ngs-01/miniconda3/envs/ngs1/share/picard-2.19.2-0/"

java -jar $picard_path/picard.jar AddOrReplaceReadGroups I=/home/ngs-01/workdir/Assignment2/2pass/Aligned.out.sam O=rg_added_sorted.bam SO=coordinate RGID=$RGID RGLB=$LB RGPL=$PL RGPU=$PU RGSM=$SM

java -jar $picard_path/picard.jar MarkDuplicates I=rg_added_sorted.bam O=dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics 
```

# Split'N'Trim and reassign mapping qualities   
# conda install -c bioconda gatk4 

```javascript
samtools faidx /home/ngs-01/workdir/sample_data/chr22_with_ERCC92.fa
gatk CreateSequenceDictionary -R /home/ngs-01/workdir/sample_data/chr22_with_ERCC92.fa -O /home/ngs-01/workdir/sample_data/chr22_with_ERCC92.dict
gatk SplitNCigarReads -R /home/ngs-01/workdir/sample_data/chr22_with_ERCC92.fa -I dedupped.bam -O split.bam
```
#indexing of BAM
picard_path=$picard_path
java -Xmx2g -jar $picard_path/picard.jar BuildBamIndex VALIDATION_STRINGENCY=LENIENT INPUT=split.bam

# Download known varinats

```javascript
wget ftp://ftp.ensembl.org/pub/grch37/current/variation/vcf/homo_sapiens/homo_sapiens-chr22.vcf.gz -O chr22.vcf.gz```
gunzip chr22.vcf.gz
```

# Select variants on chr22 and correct chr name 

```javascript
grep "^#" chr22.vcf > chr_22.vcf
grep "^22" chr22.vcf | sed 's/^22/chr22/g' >> chr_22.vcf
gatk IndexFeatureFile -F chr_22.vcf
```
# Recalibrate Bases

```javascript
gatk --java-options "-Xmx2G" BaseRecalibrator -R /home/ngs-01/workdir/sample_data/chr22_with_ERCC92.fa -I split.bam --known-sites chr_22.vcf -O Recali_base.report
```


#git clone https://github.com/PhoebeMagdy/ngs2-assignment.git  
#echo "Phoebe Magdy Abd-El Massieh" > user_info.md
#echo "phoffa_m@hotmail.com" >> user_info.md
#git add user_info.md
#git commit -m "my info"
#git push
