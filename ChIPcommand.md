```

#The location of rawdata
datadir=/mnt/haoyj/Project/WubinMa_translocation/Manuscript_20240415/Figure7/1.data/NACAPH2_HCT116_20140516
#The data2file the location: data2file.txt
#the QC of raw data
FastQCdir=/mnt/haoyj/Project/WubinMa_translocation/Manuscript_20240415/Figure7/3.bwconvert/1.qc

#Trimmed data and QC
trimedfadir=/mnt/haoyj/Project/WubinMa_translocation/Manuscript_20240415/Figure7/3.bwconvert/2.trim

####Note: The GGGG sequence not removed in the R2 files
#Mapped data to the human genome
bowtieoutdir=/mnt/haoyj/Project/WubinMa_translocation/Manuscript_20240415/Figure7/3.bwconvert/3.bowtie
hg38ref=/mnt/share/Index/hg38_bowtie_index
echo "The above mapping percentage depends on the R2 GGG number. We could use as single-end and remap to check the mapping quality"
echo "We try, but only increase 1% mapping rates"

#obtain bigwig files
bwoutdir=/mnt/haoyj/Project/WubinMa_translocation/Manuscript_20240415/Figure7/3.bwconvert/4.bw
echo "the previous command looks weird, we regenerated the bigwig files"

cd ${datadir}
#for i in `ls -1 *.fastq.gz |sed 's/_R1/\t/g'|sed 's/_R2/\t/g'|cut -f 1|sort|uniq`;
#for i in WM_240411_1_S16_L007 WM_240411_2_S17_L007 WM_240411_3_S18_L007 WM_240411_4_S19_L007 WM_240411_10_S25_L007 WM_240411_11_S26_L007 WM_240411_12_S27_L007 WM_240411_13_S28_L007 WM_240411_14_S29_L007 WM_240411_15_S30_L007;
#do
#     	bowtie2 -t -k 1 --end-to-end --sensitive -p 20 --fr --no-mixed --no-discordant -X 1000 -x ${hg38ref}/hg38_index -1 ${trimedfadir}/${i}_R1_001_val_1.fq.gz -2 ${trimedfadir}/${i}_R2_001_val_2.fq.gz 2> ${bowtieoutdir}/${i}_YJ.bowtie.stats | samtools view -q 255 -bS - | samtools sort - -o ${bowtieoutdir}/${i}.sort_YJ.bam
#	samtools index ${bowtieoutdir}/${i}.sort_YJ.bam
#	java -jar /mnt/share/software/picard/picard.jar MarkDuplicates I=${bowtieoutdir}/${i}.sort_YJ.bam O=${bowtieoutdir}/${i}.sort_YJ_rmD.bam M=${bowtieoutdir}/${i}.YJ_marked_dup_metrics.txt REMOVE_DUPLICATES=true
#	echo ${bowtieoutdir}/${i}.sort_YJ_rmD.bam
#	samtools index ${bowtieoutdir}/${i}.sort_YJ_rmD.bam
#   	bamCoverage --numberOfProcessors 20 -v -b ${bowtieoutdir}/${i}.no_dup.bam -o ${bwoutdir}/${i}_YJ.cpm.bin10.bw --normalizeUsing CPM --binSize 10 --ignoreDuplicates --smoothLength 50 --skipNAs --ignoreForNormalization chrM --extendReads --maxFragmentLength 1000
#done

echo "need rename the files"
#mv ${bwoutdir}/WM_240411_1_S16_L007_YJ.cpm.bin10.bw   ${bwoutdir}/WM_HCT116_Input_YJ.bw
#mv ${bwoutdir}/WM_240411_2_S17_L007_YJ.cpm.bin10.bw   ${bwoutdir}/WM_CAPH2AID_Input_YJ.bw
#mv ${bwoutdir}/WM_240411_3_S18_L007_YJ.cpm.bin10.bw   ${bwoutdir}/WM_HCT116_IgG_YJ.bw
#mv ${bwoutdir}/WM_240411_4_S19_L007_YJ.cpm.bin10.bw   ${bwoutdir}/WM_CAPH2IAA_HA_R1_YJ.bw
#mv ${bwoutdir}/WM_240411_5_S20_L007_YJ.cpm.bin10.bw   ${bwoutdir}/WM_CAPH2IAA_HA_R2_YJ.bw
#mv ${bwoutdir}/WM_240411_6_S21_L007_YJ.cpm.bin10.bw   ${bwoutdir}/WM_CAPH2IAA_HA_R3_YJ.bw
#mv ${bwoutdir}/WM_240411_7_S22_L007_YJ.cpm.bin10.bw   ${bwoutdir}/WM_HCT116_275A_R1_YJ.bw
#mv ${bwoutdir}/WM_240411_8_S23_L007_YJ.cpm.bin10.bw   ${bwoutdir}/WM_HCT116_275A_R2_YJ.bw
#mv ${bwoutdir}/WM_240411_9_S24_L007_YJ.cpm.bin10.bw   ${bwoutdir}/WM_HCT116_275A_R3_YJ.bw
#mv ${bwoutdir}/WM_240411_10_S25_L007_YJ.cpm.bin10.bw  ${bwoutdir}/WM-CAPH2AID_275A_R1_YJ.bw
#mv ${bwoutdir}/WM_240411_11_S26_L007_YJ.cpm.bin10.bw  ${bwoutdir}/WM-CAPH2AID_275A_R2_YJ.bw
#mv ${bwoutdir}/WM_240411_12_S27_L007_YJ.cpm.bin10.bw  ${bwoutdir}/WM-CAPH2AID_275A_R3_YJ.bw
#mv ${bwoutdir}/WM_240411_13_S28_L007_YJ.cpm.bin10.bw  ${bwoutdir}/WM-CAPH2IAA_275A_R1_YJ.bw
#mv ${bwoutdir}/WM_240411_14_S29_L007_YJ.cpm.bin10.bw  ${bwoutdir}/WM-CAPH2IAA_275A_R2_YJ.bw
#mv ${bwoutdir}/WM_240411_15_S30_L007_YJ.cpm.bin10.bw  ${bwoutdir}/WM-CAPH2IAA_275A_R3_YJ.bw

#check the correlation among replicates
FigureDir=/mnt/haoyj/Project/WubinMa_translocation/Manuscript_20240415/Figure7/3.bwconvert/6.deeptools
#check public HCT116 input datasets
HCT116input=/mnt/haoyj/Project/WubinMa_translocation/Manuscript_20240415/Figure7/3.bwconvert/PublicRefer

#generate bw files for HCT116
#bamCoverage --numberOfProcessors 20 -v -b ${HCT116input}/ENCFF040CVQ.bam -o ${HCT116input}/ENCFF040CVQ.bw --normalizeUsing CPM --binSize 10
#the public input is very different from ours, so we processed from raw fastq files
#bowtie2 -t -k 1 --end-to-end --sensitive -p 20 -x ${hg38ref}/hg38_index -U ${HCT116input}/ENCFF001HME.fastq 2> ${HCT116input}/HaoLab.bowtie.stats | samtools view -q 255 -bS - >${HCT116input}/HaoLab.Unsort.bam
#samtools sort -@ 20 -o ${HCT116input}/HaoLab.sort.bam ${HCT116input}/HaoLab.Unsort.bam
#samtools index ${HCT116input}/HaoLab.sort.bam
#bamCoverage --numberOfProcessors 20 -v -b ${HCT116input}/HaoLab.sort.bam -o ${HCT116input}/HaoLab.bw --normalizeUsing CPM --binSize 10 --ignoreDuplicates --smoothLength 50 --skipNAs --ignoreForNormalization chrM
#echo "when try our processing pipeline, the data looks good"

echo "check the correlation among replicates,includes the public data"
#multiBigwigSummary bins -b ${bwoutdir}/WM_HCT116_Input_YJ.bw ${bwoutdir}/WM_HCT116_IgG_YJ.bw ${bwoutdir}/WM_HCT116_275A_R1_YJ.bw ${bwoutdir}/WM_HCT116_275A_R2_YJ.bw ${bwoutdir}/WM_HCT116_275A_R3_YJ.bw ${HCT116input}/HaoLab_rmDup.bw --labels Input IgG Rep1 Rep2 Rep3 PublicInput -out ${FigureDir}/BW_compare_HCT116_WT.npz
#plotPCA -in ${FigureDir}/BW_compare_HCT116_WT.npz -o ${FigureDir}/BW_compare_HCT116_WT_PCA.pdf
#plotCorrelation  -in ${FigureDir}/BW_compare_HCT116_WT.npz --corMethod pearson --skipZeros -o ${FigureDir}/BW_compare_HCT116_WT_cor.pdf --whatToPlot heatmap --colorMap RdYlBu --plotNumbers

echo "Check the correlation among replicates in HCT116_AID cells"
multiBigwigSummary bins -b ${bwoutdir}/WM_CAPH2AID_Input_YJ.bw ${bwoutdir}/WM-CAPH2AID_275A_R1_YJ.bw ${bwoutdir}/WM-CAPH2AID_275A_R2_YJ.bw ${bwoutdir}/WM-CAPH2AID_275A_R3_YJ.bw ${bwoutdir}/WM-CAPH2IAA_275A_R1_YJ.bw ${bwoutdir}/WM-CAPH2IAA_275A_R2_YJ.bw ${bwoutdir}/WM-CAPH2IAA_275A_R3_YJ.bw ${bwoutdir}/WM_CAPH2IAA_HA_R1_YJ.bw ${bwoutdir}/WM_CAPH2IAA_HA_R2_YJ.bw ${bwoutdir}/WM_CAPH2IAA_HA_R3_YJ.bw --labels Input Ab_Rep1 Ab_Rep2 Ab_Rep3 IAA_Rep1 IAA_Rep2 IAA_Rep3 HA_Rep1 HA_Rep2 HA_Rep3 -out ${FigureDir}/BW_compare_HCT116_AID.npz -p 30
plotPCA -in ${FigureDir}/BW_compare_HCT116_AID.npz -o ${FigureDir}/BW_compare_HCT116_AID_PCA.pdf
plotCorrelation  -in ${FigureDir}/BW_compare_HCT116_AID.npz --corMethod pearson --skipZeros -o ${FigureDir}/BW_compare_HCT116_AID_cor.pdf --whatToPlot heatmap --colorMap RdYlBu --plotNumbers

echo "check the ChIP-seq quality"
#plotFingerprint -b ${bowtieoutdir}/WM_240411_1_S16_L007.sort_YJ_rmD.bam ${bowtieoutdir}/WM_240411_3_S18_L007.sort_YJ_rmD.bam ${bowtieoutdir}/WM_240411_7_S22_L007.sort_YJ_rmD.bam ${bowtieoutdir}/WM_240411_8_S23_L007.sort_YJ_rmD.bam ${bowtieoutdir}/WM_240411_9_S24_L007.sort_YJ_rmD.bam --labels Input IgG Rep1 Rep2 Rep3 --skipZeros --plotFile ${FigureDir}/fingerprints.png
#plotFingerprint -b ${bowtieoutdir}/WM_240411_2_S17_L007.sort_YJ_rmD.bam ${bowtieoutdir}/WM_240411_10_S25_L007.sort_YJ_rmD.bam ${bowtieoutdir}/WM_240411_11_S26_L007.sort_YJ_rmD.bam ${bowtieoutdir}/WM_240411_12_S27_L007.sort_YJ_rmD.bam ${bowtieoutdir}/WM_240411_13_S28_L007.sort_YJ_rmD.bam ${bowtieoutdir}/WM_240411_14_S29_L007.sort_YJ_rmD.bam ${bowtieoutdir}/WM_240411_15_S30_L007.sort_YJ_rmD.bam --labels Input Rep1 Rep2 Rep3 IAA_rep1 IAA_rep2 IAA_rep3 --skipZeros --plotFile ${FigureDir}/AIDcells_fingerprints.png
#plotFingerprint -b ${bowtieoutdir}/WM_240411_2_S17_L007.sort_YJ_rmD.bam ${bowtieoutdir}/WM_240411_4_S19_L007.sort_YJ_rmD.bam ${bowtieoutdir}/WM_240411_5_S20_L007.sort_YJ_rmD.bam ${bowtieoutdir}/WM_240411_6_S21_L007.sort_YJ_rmD.bam --labels Input HA_Rep1 HA_Rep2 HA_Rep3 --skipZeros --plotFile ${FigureDir}/AIDcells_HAab_fingerprints.png

echo "Call peaks"
pcoutdir=/mnt/haoyj/Project/WubinMa_translocation/Manuscript_20240415/Figure7/3.bwconvert/5.macs3
#/mnt/xuhj/software/anaconda3/bin/macs3 callpeak -t ${bowtieoutdir}/WM_240411_8_S23_L007.sort_YJ_rmD.bam -c ${bowtieoutdir}/WM_240411_1_S16_L007.sort_YJ_rmD.bam -f BAMPE -g hs --outdir ${pcoutdir} -n HCT116_WT_YJ --broad --cutoff-analysis 
#cut -f 1-6 ${pcoutdir}/HCT116_WT_YJ_peaks.broadPeak > ${pcoutdir}/HCT116_WT_YJ_peaks_broadPeak.bed

echo "Check peak quality"
#computeMatrix scale-regions -S ${bwoutdir}/WM_HCT116_275A_R2_YJ.bw ${bwoutdir}/WM_HCT116_Input_YJ.bw -R ${pcoutdir}/HCT116_WT_YJ_peaks_broadPeak.bed --beforeRegionStartLength 10000 --regionBodyLength 5000 --afterRegionStartLength 10000 --skipZeros -o ${FigureDir}/HCT116_WT.mat.gz
#plotHeatmap -m ${FigureDir}/HCT116_WT.mat.gz -out ${FigureDir}/HCT116_WT_peaks.pdf

echo "call peaks for HA-tagged files"
#/mnt/xuhj/software/anaconda3/bin/macs3 callpeak -t ${bowtieoutdir}/WM_240411_5_S20_L007.sort_YJ_rmD.bam ${bowtieoutdir}/WM_240411_6_S21_L007.sort_YJ_rmD.bam -c ${bowtieoutdir}/WM_240411_2_S17_L007.sort_YJ_rmD.bam -f BAMPE -g hs --outdir ${pcoutdir} -n HCT116_AID_HAab_YJ --broad --cutoff-analysis 
#cut -f 1-6 ${pcoutdir}/HCT116_AID_HAab_YJ_peaks.broadPeak > ${pcoutdir}/HCT116_AID_HAab_YJ_peaks_broadPeak.bed

echo "Check peak quality"
#computeMatrix scale-regions -S ${bwoutdir}/WM_CAPH2IAA_HA_R2_YJ.bw ${bwoutdir}/WM_CAPH2IAA_HA_R3_YJ.bw ${bwoutdir}/WM_CAPH2AID_Input_YJ.bw -R ${pcoutdir}/HCT116_AID_HAab_YJ_peaks_broadPeak.bed --beforeRegionStartLength 10000 --regionBodyLength 5000 --afterRegionStartLength 10000 --skipZeros -o ${FigureDir}/HCT116_AID_HAab_YJ.mat.gz
#plotHeatmap -m ${FigureDir}/HCT116_AID_HAab_YJ.mat.gz -out ${FigureDir}/HCT116_AID_HAab_YJ_peaks.pdf

echo "normalized signals"
#bigwigCompare --bigwig1 ${bwoutdir}/WM_HCT116_275A_R2_YJ.bw --bigwig2  ${bwoutdir}/WM_HCT116_Input_YJ.bw --skipZeroOverZero --operation ratio --skipNAs -p 20 -o ${bwoutdir}/WM_HCT116_WT_AbVSinput_ratio.bw -of bigwig 
#bigwigCompare --bigwig1 ${bwoutdir}/WM_HCT116_275A_R2_YJ.bw --bigwig2  ${bwoutdir}/WM_HCT116_Input_YJ.bw --skipZeroOverZero --operation subtract --skipNAs -p 20  -o ${bwoutdir}/WM_HCT116_WT_AbVSinput_subtract.bw -of bigwig

echo "Using normalized signal check peak quality"
#computeMatrix scale-regions -S ${bwoutdir}/WM_HCT116_WT_AbVSinput_ratio.bw -R ${pcoutdir}/HCT116_WT_YJ_peaks_broadPeak.bed --beforeRegionStartLength 10000 --regionBodyLength 5000 --afterRegionStartLength 10000 --skipZeros -o ${FigureDir}/HCT116_WT_AbVSinput_ratio.mat.gz -p 20
#plotHeatmap -m ${FigureDir}/HCT116_WT_AbVSinput_ratio.mat.gz -out ${FigureDir}/HCT116_WT_AbVSinput_ratio_peaks.pdf --colorMap inferno --heatmapHeight 15
#computeMatrix scale-regions -S ${bwoutdir}/WM_HCT116_WT_AbVSinput_subtract.bw -R ${pcoutdir}/HCT116_WT_YJ_peaks_broadPeak.bed --beforeRegionStartLength 10000 --regionBodyLength 5000 --afterRegionStartLength 10000 --skipZeros -o ${FigureDir}/HCT116_WT_AbVSinput_subtract.mat.gz -p 20
#plotHeatmap -m ${FigureDir}/HCT116_WT_AbVSinput_subtract.mat.gz -out ${FigureDir}/HCT116_WT_AbVSinput_subtract_peaks.pdf --colorMap inferno --heatmapHeight 20

echo "ratio normalization is better than subtract, so we chose this"
#bigwigCompare --bigwig1 ${bwoutdir}/WM_CAPH2IAA_HA_R2_YJ.bw --bigwig2 ${bwoutdir}/WM_CAPH2AID_Input_YJ.bw --skipZeroOverZero --operation ratio --skipNAs -p 20 -o ${bwoutdir}/WM_HCT116_AID_HAVSinput_R2_ratio.bw -of bigwig
#bigwigCompare --bigwig1 ${bwoutdir}/WM_CAPH2IAA_HA_R3_YJ.bw --bigwig2 ${bwoutdir}/WM_CAPH2AID_Input_YJ.bw --skipZeroOverZero --operation ratio --skipNAs -p 20 -o ${bwoutdir}/WM_HCT116_AID_HAVSinput_R3_ratio.bw -of bigwig
#computeMatrix scale-regions -S ${bwoutdir}/WM_HCT116_AID_HAVSinput_R2_ratio.bw ${bwoutdir}/WM_HCT116_AID_HAVSinput_R3_ratio.bw -R ${pcoutdir}/HCT116_AID_HAab_YJ_peaks_broadPeak.bed --beforeRegionStartLength 10000 --regionBodyLength 5000 --afterRegionStartLength 10000 --skipZeros -o ${FigureDir}/HCT116_AID_HAVSinput_ratio.mat.gz -p 20
#plotHeatmap -m ${FigureDir}/HCT116_AID_HAVSinput_ratio.mat.gz -out ${FigureDir}/HCT116_AID_HAVSinput_ratio_peaks.pdf --colorMap inferno --heatmapHeight 20


```
