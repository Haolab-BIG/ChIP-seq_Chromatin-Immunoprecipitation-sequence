# ChIP-seq, Chromatin Immunoprecipitation sequence
For everyone who has a new face to HaoLab, we make this pipeline/protocol to standardize downstream analysis pipeline for ChIP-seq data.
## Part I Introduction
### Workflow
Here stands an throughout workflow of ChIP-seq data analysis.

![image](https://github.com/Haolab-BIG/ChIP-seq-Chromatin-Immunoprecipitation-sequence/blob/main/Figure/Firgue1_workflow.png)

As illustrated in the figure, \
(i) yellow circles represent the steps where commands need to be entered; \
(ii) pink dashed rectangular boxes represent the output results after processing at each step. 

We will proceed by structuring our workflow according to (ii).
### File Structure
Here stands an throughout file structure of ChIP-seq data analysis.
* *You can decide what your structure looks like, which makes it more efficient to work.* 

```
ChIP-seq 
├─ 1.raw_data 
│    └─ QC 
├─ 2.clean_data 
├─ 3.bam 
│    ├─ map_statistics
│    ├─ rmD_bam
│    └─ uniq_bam
├─ 4.fingerprint 
├─ 5.bw 
│    ├─ PCA 
│    └─ correlation 
└─ 6.peaks 
       └─ heatmap
```
### Conda Environment
You can configure a Conda environment named 'ChIP-seq' using the following code, which includes the essential software for ChIP-Seq analysis.

```
conda create -n ChIP-seq fastqc Trimmomatic bowtie2 samtools Picard deeptools MACS3
```

## Part II FASTQ2BAM
In this section, you will convert the raw FASTQ files into BAM files, which can be used for subsequent analysis.
### Raw Data Quality Check(QC)
You can perform quality check on the raw data to assess the sequencing quality.

```
fastqc -o /ChIP-seq/1.raw_data/QC/ /ChIP-seq/1.raw_data/raw_data_1.fq.gz /ChIP-seq/1.raw_data/raw_data_2.fq.gz
```

### Trimming after QC
Based on the quality control results, you can perform appropriate trimming on the raw data.

```
NO SPECIFIC CODE
```

### Mapping clean data to genome, filtering and remove duplicates
Align the cleaned data obtained in the previous step to the reference genome of the corresponding species. After filtering for high-quality sequences and removing duplicates, the resulting data will be prepared for further analysis.

#### Mapping, filtering
```
bowtie2 -t -k 1 --end-to-end --sensitive -p 20 --fr --no-mixed --no-discordant -X 1000 -x referrence_species_index -1 /ChIP-seq/2.clean_data/clean_data_1.fq.gz -2 /ChIP-seq/2.clean_data/clean_data_2.fq.gz 2> /ChIP-seq/3.bam/map_statistics/ | samtools view -q 255 -bS - | samtools sort - -o /ChIP-seq/3.bam/uniq_bam/uniq.sorted.bam
```
#### BAM indexing
```
samtools index /ChIP-seq/3.bam/uniq_bam/uniq.sorted.bam
```
#### Remove duplicates
```
java -jar picard.jar MarkDuplicates I=/ChIP-seq/3.bam/uniq_bam/uniq.sorted.bam O=/ChIP-seq/3.bam/rmD_bam/rmD.bam M=/ChIP-seq/3.bam/rmD_bam/marked_dup_metrics.txt REMOVE_DUPLICATES=true
```
## Part III
## Part IV
