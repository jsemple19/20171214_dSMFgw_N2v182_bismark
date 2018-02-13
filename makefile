#!/usr/bin/make -f
# mapping dSMF-gw sequences with bismark
## required software: fastqc, trimmomatic, bismark, bowtie2, samtools, picard, qualimap

###############################
########### VARIABLES #########
###############################

genomeDIR := ${HOME}/Documents/MeisterLab/GenomeVer/sequence
genomefile := ${genomeDIR}/c_elegans.PRJNA13758.WS250.genomic.fa
trimmomaticDIR := ${HOME}/Trimmomatic-0.36
bismarkDIR := ${HOME}/mySoftware/Bismark_v0.19.0
methIndGenomeFiles := $(addprefix ${genomeDIR}/Bisulfite_Genome/CT_conversion/, BS_CT.1.bt2	BS_CT.3.bt2	BS_CT.rev.1.bt2	genome_mfa.CT_conversion.fa BS_CT.2.bt2	BS_CT.4.bt2 BS_CT.rev.2.bt2) \
	$(addprefix ${genomeDIR}/Bisulfite_Genome/GA_conversion/, BS_GA.1.bt2 BS_GA.3.bt2 BS_GA.rev.1.bt2 genome_mfa.GA_conversion.fa BS_GA.2.bt2 BS_GA.4.bt2 BS_GA.rev.2.bt2)
bname :=  $(addprefix 180126_SNK268_A_L001_JIB-, 1L 2L 3L 4L)
longbname := $(addsuffix _R1, $(bname)) $(addsuffix _R2, $(bname))


#list of the final output files
objects :=  ${methIndGenomeFiles} \
	$(addsuffix .filt.bam, $(addprefix aln/, $(bname))) 

#list of the various reports and stats produced during analysis
statsObjects := $(addsuffix _fastqc.html, $(addprefix rawData/fastQC/, $(longbname))) \
	$(addsuffix _forward_paired_fastqc.html, $(addprefix trim/fastQC/, $(bname))) \
	$(addsuffix _forward_unpaired_fastqc.html, $(addprefix trim/fastQC/, $(bname))) \
	$(addsuffix _reverse_paired_fastqc.html, $(addprefix trim/fastQC/, $(bname))) \
	$(addsuffix _reverse_unpaired_fastqc.html, $(addprefix trim/fastQC/, $(bname))) \
	$(addprefix aln/prefilt/report_, $(addsuffix _flagstats.txt, $(bname))) \
	$(addprefix aln/prefilt/report_, $(addsuffix _stats.txt, $(bname))) \
	$(addprefix aln/prefilt/report_, $(addsuffix _flagstats.txt, $(bname))) \
	$(addprefix aln/prefilt/report_, $(addsuffix _picard_insert_size_metrics.txt, $(bname))) \
	$(addprefix aln/prefilt/report_, $(addsuffix _picard_insert_size_histogram.pdf, $(bname))) \
	$(addprefix aln/prefilt/report_, $(addsuffix _qualimap.pdf, $(bname))) \
	$(addprefix aln/postfilt/report_, $(addsuffix _picard.txt, $(bname))) \
	$(addprefix aln/postfilt/report_, $(addsuffix _flagstats.txt, $(bname))) \
	$(addprefix aln/postfilt/report_, $(addsuffix _stats.txt, $(bname))) \
	$(addprefix aln/postfilt/report_, $(addsuffix _flagstats.txt, $(bname))) \
	$(addprefix aln/postfilt/report_, $(addsuffix _picard_insert_size_metrics.txt, $(bname))) \
	$(addprefix aln/postfilt/report_, $(addsuffix _picard_insert_size_histogram.pdf, $(bname))) \
	$(addprefix aln/postfilt/report_, $(addsuffix _qualimap.pdf, $(bname))) 
	
#list of files to delete at the end of the analysis
intermediateFiles := $(addsuffix _forward_paired.fq.gz, $(addprefix trim/, $(bname))) \
	$(addsuffix _forward_unpaired.fq.gz, $(addprefix trim/, $(bname))) \
	$(addsuffix _reverse_paired.fq.gz, $(addprefix trim/, $(bname))) \
	$(addsuffix _reverse_unpaired.fq.gz, $(addprefix trim/, $(bname))) \
	$(addsuffix .bam, $(addprefix aln/, $(bname))) \
	$(addsuffix .dup.bam, $(addprefix aln/, $(bname)))

#list of secondary files to keep
secondaryFiles := $(addsuffix .sorted.bam, $(addprefix aln/, $(bname)))


###############################
########### RULES  ############
###############################

all: $(objects) $(statsObjects) $(secondaryFiles)

#use cleanall when you want to force rerunning all the analysis
cleanall:
	rm -f $(intermediateFiles)
	rm -f $(secondaryFiles)
	
#use clean if the intermediate files are not properly removed (should not be required)
clean:
	rm -f $(intermediateFiles)

cleanall4rerun:
	rm -f $(objects)
	rm -f $(statsObjects)
	rm -f $(secondaryFiles)

.PHONY: all clean cleanall cleanall4rerun
.INTERMEDIATE: $(intermediateFiles)
.SECONDARY: $(secondaryFiles)


#######################################################
## get initial read stats                            ##
#######################################################

#run fastqc on downloaded sequences
rawData/fastQC/%_fastqc.html: rawData/%.fastq.gz
	mkdir -p rawData/fastQC
	fastqc $^ -o rawData/fastQC 


#######################################################
## quality trim reads with Trimmomatic               ##
#######################################################

# use trimmomatic to trim
trim/%_forward_paired.fq.gz trim/%_forward_unpaired.fq.gz trim/%_reverse_paired.fq.gz trim/%_reverse_unpaired.fq.gz: rawData/%_R1.fastq.gz rawData/%_R2.fastq.gz
	mkdir -p trim
	java -jar ${trimmomaticDIR}/trimmomatic-0.36.jar PE rawData/$*_R1.fastq.gz rawData/$*_R2.fastq.gz trim/$*_forward_paired.fq.gz trim/$*_forward_unpaired.fq.gz trim/$*_reverse_paired.fq.gz trim/$*_reverse_unpaired.fq.gz ILLUMINACLIP:${trimmomaticDIR}/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> trim/report_$*_trimmomatic.txt
	# redo fastQC on trimmed reads
	mkdir -p trim/fastQC
	fastqc trim/$*_forward_paired.fq.gz trim/$*_forward_unpaired.fq.gz trim/$*_reverse_paired.fq.gz trim/$*_reverse_unpaired.fq.gz -o trim/fastQC


#######################################################
## align to genome with BWA-meth and convert to bam  ##
#######################################################
	
# convert and index genome file for bismark alignment
${methIndGenomeFiles}: ${genomefile}	
	${bismarkDIR}/bismark_genome_preparation --bowtie2 --verbose ${genomeDIR}

# align sequences to meth converted genome with bismark
aln/%.bam: trim/%_forward_paired.fq.gz trim/%_reverse_paired.fq.gz ${methIndGenomeFiles}
	mkdir -p aln
	${bismarkDIR}/bismark --bowtie2 -o aln -X 600 --genome ${genomeDIR} -1 trim/$*_forward_paired.fq.gz -2 trim/$*_reverse_paired.fq.gz
	mv aln/$*_forward_paired_bismark_bt2_pe.bam aln/$*.bam #simplify names

# use samtools to sort bam
aln/%.sorted.bam: aln/%.bam
	samtools view -u $^ | samtools sort -o $@


#######################################################
## Get alignment stats                               ##
#######################################################

# 	get alignment stats
aln/prefilt/report_%_flagstats.txt: aln/%.sorted.bam
	mkdir -p aln/prefilt
	samtools flagstat $^ > $@

aln/prefilt/report_%_stats.txt: aln/%.sorted.bam
	samtools stats $^ > $@

# Get insert size statistics and plots with picard and qualimap (set path to $QUALIMAP in .bash_profile)
aln/prefilt/report_%_picard_insert_size_metrics.txt aln/prefilt/report_%_picard_insert_size_histogram.pdf \
	aln/prefilt/report_%_qualimap.pdf: aln/%.sorted.bam ${PICARD} ${QUALIMAP}
	java -Xms1g -Xmx5g -jar ${PICARD} CollectInsertSizeMetrics I=aln/$*.sorted.bam \
	O=aln/prefilt/$*_picard_insert_size_metrics.txt \
    H=aln/prefilt/$*_picard_insert_size_histogram.pdf
	${QUALIMAP} bamqc -bam aln/$*.sorted.bam -c -outdir aln/prefilt -outfile $*_report_qualimap.pdf -outformat PDF


#######################################################
## Mark duplicates and filter reads. Then redo stats ##
#######################################################

# mark duplicates with picard (path to picard should be set in $PICARD variable in .bash_profile or in session)
aln/%.dup.bam aln/postfilt/report_%_picard.txt: aln/%.sorted.bam ${PICARD}
	mkdir -p aln/postfilt
	java -Xmx5g -jar ${PICARD} MarkDuplicates I=aln/$*.sorted.bam O=aln/$*.dup.bam M=aln/postfilt/report_$*_picard.txt

# 	remove mitochondrial reads
aln/%.filt.bam: aln/%.dup.bam
	samtools view -q 30 -F 1804 -b $^ > $@ 	
	# NOTE: sam flag 1804 means the following:
	# read unmapped
	# mate unmapped
	# not primary alignment
	# read fails platform/vendor quality checks
	# read is PCR or optical duplicate

# 	get alignment stats again post-filtering
aln/postfilt/report_%_flagstats.txt: aln/%.filt.bam
	samtools flagstat $^ > $@

aln/postfilt/report_%_stats.txt: aln/%.filt.bam
	samtools stats $^ > $@

# Get insert size statistics and plots with picard and qualimap post-filtering
aln/postfilt/report_%_picard_insert_size_metrics.txt aln/postfilt/report_%_picard_insert_size_histogram.pdf \
	aln/postfilt/report_%_qualimap.pdf: aln/%.filt.bam ${PICARD} ${QUALIMAP}
	java -Xms1g -Xmx5g -jar ${PICARD} CollectInsertSizeMetrics I=aln/$*.filt.bam \
	O=aln/postfilt/$*_picard_insert_size_metrics.txt \
    H=aln/postfilt/$*_picard_insert_size_histogram.pdf
	${QUALIMAP} bamqc -bam aln/$*.filt.bam -c -outdir aln/postfilt -outfile $*_report_qualimap.pdf -outformat PDF
