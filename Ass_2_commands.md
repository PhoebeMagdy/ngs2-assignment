# Variant Calling for RNA seq Data
## Installing STAR for alignning RNA reads

```javascript
conda install -c bioconda star 
conda install -c bioconda/label/cf201901 star 
```

## Generating genome index

``` javascript
genomeDir=/home/ngs-01/workdir/Assignment2/chr22_with_ERCC92
mkdir $genomeDir
STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles /home/ngs-01/workdir/Assignment2/chr22_with_ERCC92.fa  --runThreadN 1 --sjdbGTFfile /home/ngs-01/workdir/sample_data/chr22_with_ERCC92.gtf
```

## Alignment jobs were executed as follows:

``` javascript
runDir=/home/ngs-01/workdir/Assignment2/1pass
mkdir $runDir
cd $runDir
STAR --genomeDir $genomeDir --readFilesIn /home/ngs-01/workdir/Assignment2/ngs2-assignment-data/SRR8797509_1.part_001.part_001.fastq /home/ngs-01/workdir/Assignment2/ngs2-assignment-data/SRR8797509_2.part_001.part_001.fastq --runThreadN 4
```

## For the 2-pass STAR, a new index is then created using splice junction information contained in the file sjdbList.out.tab from the first pass:

```javascript
genomeDir=/home/ngs-01/workdir/Assignment2/chr22_with_ERCC92_2pass
mkdir $genomeDir
STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles /home/ngs-01/workdir/Assignment2/chr22_with_ERCC92.fa    --sjdbFileChrStartEnd /home/ngs-01/workdir/Assignment2/chr22_with_ERCC92/sjdbList.out.tab --sjdbOverhang 75 --runThreadN 1
```



    
    
#git clone https://github.com/PhoebeMagdy/ngs2-assignment.git  
#echo "Phoebe Magdy Abd-El Massieh" > user_info.md
#echo "phoffa_m@hotmail.com" >> user_info.md
#git add user_info.md
#git commit -m "my info"
#git push