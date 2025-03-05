# ChIP-seq, Chromatin Immunoprecipitation sequencing
For everyone who has a new face to HaoLab, we make this pipeline/protocol to standardize downstream analysis pipeline for ChIP-seq data.

## Part I Introduction
### i. Workflow
Here stands an throughout workflow of ChIP-seq data analysis.

![image](https://github.com/Haolab-BIG/ChIP-seq-Chromatin-Immunoprecipitation-sequence/blob/main/Figure/Firgue1_workflow.png)

As illustrated in the figure, \
(i) yellow circles represent the steps where commands need to be entered; \
(ii) pink dashed rectangular boxes represent the output results after processing at each step. 

We will proceed by structuring our workflow according to (ii).
### ii. File Structure
Here stands an throughout file structure of ChIP-seq data analysis.
* *You can decide what your structure looks like, which makes it more efficient to work.* 

```
ChIP-seq 
├─ 1.raw_data 
│    └─ QC 
├─ 2.clean_data 
├─ 3.bam 
│    ├─ map_statistics
│    ├─ uniq_bam
│    └─ rmD_bam
├─ 4.fingerprint 
├─ 5.bw
│    ├─ multiBigWigSummary
│    ├─ PCA 
│    └─ correlation 
└─ 6.peaks
       ├─ computeMatrix 
       └─ heatmap
```

### iii. Conda Environment
You can configure a Conda environment named 'ChIP-seq' using the following code, which includes the essential software for ChIP-Seq analysis. 

```
conda create -n ChIP-seq fastqc Trimmomatic bowtie2 samtools Picard deeptools MACS3
```

## Part II Generation of Data for Analysis: FASTQ2BAM
In this section, you will convert the raw FASTQ files into BAM files, which can be used for subsequent analysis.

### i. Raw Data Quality Check(QC)
You can perform quality check on the raw data to assess the sequencing quality.

```
fastqc -o /ChIP-seq/1.raw_data/QC/ /ChIP-seq/1.raw_data/raw_data_1.fq.gz /ChIP-seq/1.raw_data/raw_data_2.fq.gz
```

### ii. Trimming After QC
Based on the quality control results, you can perform appropriate trimming on the raw data.

```
java -jar trimmomatic-0.39.jar PE -threads 20 -phred33 \
ChIP-seq/1.raw_data/_raw_1.fq.gz ChIP-seq/1.raw_data/_raw_2.fq.gz \
ChIP-seq/2.clean_data/1.paired.fq.gz ChIP-seq/2.clean_data/1.unpaired.fq.gz \
ChIP-seq/2.clean_data/2.paired.fq.gz ChIP-seq/2.clean_data/2.unpaired.fq.gz \
ILLUMINACLIP:Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:36 HEADCROP:12
```

### iii. Mapping Clean Data to Genome, Filtering and Remove Duplicates
Align the cleaned data obtained in the previous step to the reference genome of the corresponding species. After filtering for high-quality sequences and removing duplicates, the resulting data will be prepared for further analysis.

#### 1. Mapping, Filtering
Align the quality-controlled FASTQ files to the reference genome, filter for high-quality aligned sequences, and perform sorting.

```
bowtie2 -t -k 1 --end-to-end --sensitive -p 20 --fr --no-mixed --no-discordant -X 1000 -x referrence_species_index -1 /ChIP-seq/2.clean_data/clean_data_1.fq.gz -2 /ChIP-seq/2.clean_data/clean_data_2.fq.gz 2> /ChIP-seq/3.bam/map_statistics/ | samtools view -q 255 -bS - | samtools sort - -o /ChIP-seq/3.bam/uniq_bam/uniq.sorted.bam
```

#### 2. BAM Indexing
Index the sorted files.

```
samtools index /ChIP-seq/3.bam/uniq_bam/uniq.sorted.bam
```

#### 3. Remove Duplicates
Remove duplicate sequences.

```
java -jar picard.jar MarkDuplicates I=/ChIP-seq/3.bam/uniq_bam/uniq.sorted.bam O=/ChIP-seq/3.bam/rmD_bam/.rmD.bam M=/ChIP-seq/3.bam/rmD_bam/marked_dup_metrics.txt REMOVE_DUPLICATES=true
```

### iv. Peakcalling
Identify regions of read enrichment based on the alignment results.

```
macs3 callpeak -t /ChIP-seq/3.bam/rmD_bam/Treated.rmD.bam -c /ChIP-seq/3.bam/rmD_bam/Ctrl.rmD.bam -f BAMPE -g species --outdir /ChIP-seq/6.peaks/ -n prefiix --broad --cutoff-analysis
```

#### 1. Peaks Quality Check
Assess the quality of the peaks.

```
computeMatrix scale-regions -S /ChIP-seq/5.bw/.bw -R /ChIP-seq/6.peaks/prefix_peaks_broadPeak.bed --beforeRegionStartLength 10000 --regionBodyLength 5000 --afterRegionStartLength 10000 --skipZeros -o /ChIP-seq/6.bw/commputeMatrix/.mat.gz -p 20
```

* *'.bw' file comes from Part (IV.i)* *

### OUTPUT 1
After completing the above steps, you will obtain a table similar to the one shown in 'example_report.ppt', slide (2), with columns including:

**DataName, FileName, RawData, FilterData(QC), CleanReads%, UniqueMapped%, Duplicate%, UnDuplicated, and Peaks.**


## Part III ChIP-seq Quality Check
you will obtain a plot similar to the one shown in 'example_report.ppt', slide (3,4).

```
plotFingerprint -b /ChIP-seq/3.bam/rmD_bam/.rmD.bam --labels label1 label2 --skipZeros --plotFile /ChIP-seq/4.fingerprint/.png
```

### OUTPUT 2

**By generating a fingerprint plot, you can identify the types of peaks, assess antibody strength to verify the effectiveness of immunoprecipitation (IP).**


## Part IV Visualization on Genome Browser & Check the Correlation in Samples
Convert the alignment results into files with lower resolution but smaller size for easier browsing across the genome; simultaneously, assess the correlation among the experimental samples.

### i. Visualization on Genome Browser: BAM2BW
You can visualize these results in the UCSC genome browser or IGV genome browser.

```
bamCoverage --numberOfProcessors 20 -v -b /ChIP-seq/3.bam/rmD_bam/.rmD.bam -o /ChIP-seq/5.bw/.bw --normalizeUsing CPM --binSize 10 --ignoreDuplicates --smoothLength 50 --skipNAs --ignoreForNormalization chrM
```

### ii. Check the Correlation in Samples
you will obtain a plot similar to the one shown in 'example_report.ppt', slide (5,6)

```
multiBigwigSummary bins -b /ChIP-seq/5.bw/.bw --labels label1 label2 label3 -out /ChIP-seq/5.bw/multiBigWigSummary/.npz
```

```
plotPCA -in /ChIP-seq/5.bw/multiBigWigSummary/.npz -o /ChIP-seq/5.bw/PCA/PCA.pdf
```

```
plotCorrelation  -in /ChIP-seq/5.bw/multiBigWigSummary/.npz --corMethod pearson --skipZeros -o /ChIP-seq/5.bw/PCA/correlation/cor.pdf --whatToPlot heatmap --colorMap RdYlBu --plotNumbers
```

### OUTPUT 3

**By generating Principal Component Analysis (PCA) plots and correlation heatmaps, you can assess biological reproducibility, experimental quality, and the extent of batch effects.**
