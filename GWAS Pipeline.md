# GWAS Pipeline

# Upload raw data files

scp -r /local/directory/path/to/files/ user@qb.loni.org:/path/to/folder

# Unzip files

gunzip -r /path/to/files/*gz.fq

# Concatenate files from the same sample 

cat samplename_lane_number samplename_lane_number > newsamplename.fq

##lane number will be different since reads were run on different lanes but the last number should be the same

# Run FastQC on all samples to check quality of reads

for i in /path/to/files/*.fq do /path/to/FastQC/fastqc $i; done

# FastQC creates fastqc.html files that need to be downloaded to computer to view 

scp user@qb.loni.org:/path/to/files/*.html /local/directory/
##scp securely copies files

# had high quality reads so didn't trim but would trim here 

# Quality filter reads using FastQ

module load fastx_toolkit/0.0.13.2/INTEL-14.0.2
for i in /path/to/files/*.fq; do fastq_quality_filter -Q33 -q20 -p 98 -i -o{i%.fq}.qual.fq; done

##path to the version of fastx toolkit installed on LONI 
##-Q sequence type, -q minimum quality score to keep, -p minimum percent of bases that have to have -q quality 
##output creates .qual.fq files (default is STDOUT)

# Index Reference Genome using BWA Index, Samtools and GATK

##Indexing the genome is like a book index, used to make alignment easier later 

# BWA Index
module load bwa/0.7.15/INTEL-14.0.2 
bwa index -a bwtsw /path/to/reference/.fa

##path to the version of BWA installed on LONI 
##-a algorithm to construct index from, bwtsw algorithm implemented in BWT-SW which is for large genomes, just need -a for small genomes 

# Samtools
source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate samtools-1.11
/path/to/conda-env/envs/samtools-1.11/bin/samtools faidx /path/to/reference/.fa

##Samtools path from download via conda
##Faidx will index input file and create a .fai file


# GATK
module load java/1.8.0 
/path/to/gatk-4.1.2.0/gatk --java-options "-Xmx2G" CreateSequenceDictionary -R /path/to/reference/.fa

##path to java version installed on LONI 
##-Xmx2G sets the initial and maximum heap size available to improve performance, CreateSequenceDictionary creates .dict file, -R reference


# Align reads to reference using BWA mem

module load bwa/0.7.15/INTEL-14.0.2
for i in /path/to/files/*.fq; do bwa mem -v 3 -M -P -a -t 10 
/path/to/reference/.fa 
$i ${i/1.fq/2.fq} > ${i%1.fq}PE.sam 2> ${i%1.fq}PE.mem.log; done

##path to the version of BWA installed on LONI
##bwa mem aligns 70bp-1Mbp sequences, -v verbose level (a value 0 for disabling all the output to stderr; 1 for outputting errors only; 2 for warnings and errors; 3 for all normal messages; 4 or higher for debugging), -M mark shorter split hits as secondary (for Picard compatibility), -P paired-end mode, -a output all found alignments, -t number of threads 
##output creates PE.sam and PE.mem.log files


# Create BAM files for downstream analyses using Samtools View

source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate samtools-1.11
/path/to/conda-env/envs/samtools-1.11/bin/samtools
for i in /path/to/files/*PE.sam; do samtools view -q 20 -bt /path/to/reference/.fa -o ${i%sam}bamq20 $i; done

##Samtools path from download via conda
##-q skip alignments with mapping quality less than specified number, -bt output in BAM format
##output creates bamq20 files

# Assign all reads to to read-group for GATK using Picard tools

cd $PBS_O_WORKDIR
mkdir -p tmp

module load  java/1.8.0

for i in /path/to/files/*.bamq20;
do
        java -Dpicard.useLegacyParser=false -Xmx2g -jar /path/to/picard/picard.jar \
        AddOrReplaceReadGroups -I $i -O ${i%PE.bamq20}.tag.bam -MAX_RECORDS_IN_RAM 1000000 -TMP_DIR $PWD/tmp \
        -SO coordinate -ID ${i%PE.bamq20} -LB 1 -PL illumina -PU 1 -SM ${i%PE.bamq20};
done

##path to java version installed on LONI
##-Dpicard.useLegacyParser=false for GATK 4.0 to run Picard tools from GATK, -Xmx2G sets the initial and maximum heap size available to improve performance,AddOrReplaceReadGroups assigns reads to read group, -I input file, -O output file, -MAX_RECORDS_IN_RAM specifies the number of records stored in RAM before storing in disk, TMP_DIR $PWD/tmp creates a temorary directory called tmp, -SO sort order to ouput in (options:unsorted,queryname,coordinate,duplicate,unknown), -ID read-group id, -LB read group library, -PL read group platform, -PU read group platform unit, -SM read group sample name
##creates .tag.bam files


# Mark PCR duplicates using Picard tools

module load  java/1.8.0
cd $PBS_O_WORKDIR
mkdir -p tmp

for i in /path/to/files/.tag.bam; 
do
        java -Xmx4g -jar /path/to/picard/picard.jar \
        MarkDuplicates INPUT=$i OUTPUT=${i%.tag.bam}.rmdup.bam MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=6000 MAX_RECORDS_IN_RAM=1000000 TMP_DIR=$PWD/tmp \
        METRICS_FILE=${i%.tag.bam}.rmdup.metrics ASSUME_SORTED=true;
done

##path to java version installed on LONI
##cd change directory, mkdir -d make directory and if neccessary, all parent directories, -Xmx2G sets the initial and maximum heap size available to improve performance, MarkDuplicates identifies duplicate reads, MAX_FILE_HANDLES_FOR_READ_ENDS_MAP maximum number of file handles to keep open when spilling read ends to disk (default is 8000), MAX_RECORDS_IN_RAM when writing files that need to be sorted this will specify the number of records stored in RAM before spilling to disk, TMP_DIR $PWD/tmp creates a temorary directory called tmp, METRICS_FILE File to write duplication metrics to, ASSUME_SORTED If true, assume that the input file is coordinate sorted 
##creates .rmdup.bam


# Index filtered files using Samtools

source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate samtools-1.11
/path/to/conda-env/envs/samtools-1.11/bin/samtools

for i in /path/to/files/*.rmdup.bam; do samtools index $i; done

##Samtools path from download via conda
##Samtools index indexes files for realignment 


# Create GVCF files using GATK

module load java/1.8.0

/path/to/gatk-4.1.2.0/gatk --java-options "-Xmx4g" HaplotypeCaller --ERC GVCF -R /path/to/reference/.fa  -I /path/to/file/file.rmdup.bam -O /path/to/file/file.g.vcf --min-base-quality-score 20

##path to java version installed on LONI
##-Xmx2G sets the initial and maximum heap size available to improve performance, HaplotypeCaller --ERC GVCF creates file with raw unfiltered SNP and indel calls to genotype, --min-base-quality-score minimum quality required to consider a base for calling


# Combine GVCF files using GATK

module load java/1.8.0

/path/to/gatk-4.1.2.0/gatk CombineGVCFs -R /path/to/reference/.fa  -V /path/to/file1/file1.g.vcf -V /path/to/file2/file2.g.vcf -O /path/to/file/combined.g.vcf 

##path to java version installed on LONI
##CombineGVCFs merge GVCF files into a single GVCF files, -R reference file, -V variant, -O output file


# Genotype GVCF file

module load java/1.8.0

/path/to/gatk-4.1.2.0/gatk --java-options "-Xmx4g" GenotypeGVCFs -R /path/to/reference/.fa  -V /path/to/file/combined.g.vcf -O /path/to/file/combined.vcf 

##path to java version installed on LONI
##-Xmx2G sets the initial and maximum heap size available to improve performance, GenotypeGVCFs performs joint genotyping of all samples within file, -R reference file, -V variant, -O output file
 
# Select SNP and MNP variants using GATK

module load java/1.8.0

/path/to/gatk-4.1.2.0/gatk SelectVariants --variant /path/to/file/combined.vcf -R /path/to/reference/.fa --output /path/to/file/file2.vcf --select-type-to-include SNP --select-type-to-include MNP --exclude-non-variants true --set-filtered-gt-to-nocall true 

##path to java version installed on LONI
##SelectVariants selects a subset of variants from VCF file, --variant VCF file with variants, -R reference file, --output output file --select-type-to-include options: INDEL, SNP, MIXED, MNP, SYMBOLIC, NO_VARIATION


# Quality filter variants and flag variants that don't meet criteria using GATK

module load java/1.8.0

/path/to/gatk-4.1.2.0/gatk VariantFiltration --variant /path/to/file/file2.vcf --output /path/to/file/file3.vcf -R /path/to/reference/.fa \
--filter-name "ReadPosRankSum_filter" \
--filter-expression "ReadPosRankSum < -8.0" \
--filter-name "MQRankSum_filter" \
--filter-expression "MQRankSum < -12.5" \
--filter-name "FS_filter" \
--filter-expression "FS > 60.0" \
--filter-name "QD_filter" \
--filter-expression "QD < 2.0" \
--genotype-filter-name "DP8filter" \
--genotype-filter-expression "DP < 8" 2>/dev/null 

##path to java version installed on LONI
##VariantFiltration filter varant calls, --variant VCF file with variants, --output output file, -R reference file, ReadPosRankSum compares whether positions of reference and alternate alleles are different within reads, MQRankSum compares mapping qualities of reads supporting the reference allele and alternate allele, FS Phred-scaled probability that there is strand bias at the site. FS value will be close to 0 when there is little to no strand bias at the site, QD is the variant confidence divided by the unfiltered depth of samples. This normalizes the variant quality in order to avoid inflation caused when there is deep coverage, DP is genotype depth of coverage

##I didn't know this until after I ran this code but: For filtering purposes it is better to use QD than either QUAL or DP directly (maybe the DP filter should be removed?)


# Remove flagged variants using GATK

module load java/1.8.0

/path/to/gatk-4.1.2.0/gatk SelectVariants -R /path/to/reference/.fa --variant /path/to/file/file3.vcf --output /path/to/file/file4.vcf --set-filtered-gt-to-nocall true

##path to java version installed on LONI
##SelectVariants selects a subset of variants from VCF file, -R reference file,--variant VCF file with variants,--output output file   


# View general stats of VCF file using VCFtools

export PERL5LIB=/path/to/vcftools/src/perl
/path/to/vcftools/bin/vcf-stats /path/to/file/.vcf  > output.txt

##PERL5LIB allows you to use perl through vcftools


# Remove poor sequences 

export PERL5LIB=/path/to/vcftools/src/perl
/path/to/vcftools/bin/vcftools --vcf /path/to/file/file.vcf --missing-indv --out /path/to/file/file

##--vcf VCF file, --missing-indv creates an .imiss file which has 5 columns: N_MISS= number of sites the individual does not have data for, F_MISS= frequency of missing data for the individual, create a text file with individuals to remove based on missing data, --out output file 

export PERL5LIB=/path/to/vcftools/src/perl
/path/to/vcftools/bin/vcftools --vcf /path/to/file/file.vcf --remove file.txt --out /path/to/file.vcf

##--vcf VCF file, --remove text file with individuals to remove, --out output file

# Create chromosome mapping file for VCF to PLINK conversion 

##List all the chromosomes or contigs in your VCF file, one chromosome or contig per line. I had to edit my VCF file because it didn't like the characters and underscores in the names of my contigs

I used: 
sed 's/part_to_remove//g' 
to remove the part of the contig name that PLINK didn't like (saved as a new vcf file), then copied the new names into my chromosome mapping file 

e.g

000001
000002
000003
000004

# Convert VCF to PLINK file for analyses

export PERL5LIB=/path/to/vcftools/src/perl
/path/to/vcftools/bin/vcftools --vcf /path/to/file/file.vcf --chrom-map  /path/to/file/file.txt --out /path/to/file/output  --plink

##--vcf VCF file, --chrom-map chromosome mapping file, --out output file, --plink convert to plink file

# Create phenotype file for PLINK association test

##Text file format uses FULL NAME from VCF file. My VCF file named all my individuals the full path (pretty sure that's the default) to the file so, double check what your individual names are in the VCF file you will be using for your association test.

##3 column format with Family ID (FID), Individual ID (IID) and phenotype [1= unaffected (or can be control), 2= affected (or can be case), 0= missing data], tab delimited columns 

##If there is no Family ID you can copy the Individual ID into that column or remove it completely (I haven't tried that though)

path/to/file    path/to/file    1,2 or 0

# Run association test in PLINK 

source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate plink-1.90

/path/to/miniconda3/bin/plink --file /path/to/file/file.plink --allow-extra-chr --allow-no-sex --pheno /path/to/file/file.txt --assoc --out /path/to/output

--file plink file, --allow-extra-chr allow extra chromosomes [permit unrecognized chromosome codes, treating their variants as unplaced], --allow-no-sex Do not force ambiguous-sex phenotypes to missing, --pheno phenotype file, --assoc basic association test, --out output  

# Sort association by p-value

sort -k 9,9 file.assoc > new_file

##p-values are in the 9th column 

# Convert PLINK files to .bed 

source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate plink-1.90

/home/gabby297/miniconda3/bin/plink --file /path/to/file/file --allow-extra-chr --allow-no-sex --pheno /path/to/file/file.txt --chr-set 40 no-xy --make-bed --out output

--file plink file, --allow-extra-chr allow extra chromosomes [permit unrecognized chromosome codes, treating their variants as unplaced], --allow-no-sex Do not force ambiguous-sex phenotypes to missing, --pheno phenotype file, --chr-set 40 changes the chromosome set. The first parameter specifies the number of diploid autosome pairs if positive or haploid chromosomes if negative. Given diploid autosomes the remaining modifiers let you indicate the absence of specific non-autosomal chromosomes. When there are n autosome pairs the X chromosome is assigned numeric code n+1, Y is n+2, XY (pseudo-autosomal region of X) is n+3, and MT (mitochondria) is n+4,  --make-bed generate .bed binary file, --out output
 
# Run PCA in PLINK 

source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate plink-1.90

/path/to/miniconda3/bin/plink --bfile /path/to/file.bed --allow-extra-chr --allow-no-sex --pheno /path/to/file.txt --chr-set 40 no-xy --pca --out output

--bfile bed file, --allow-extra-chr allow extra chromosomes [permit unrecognized chromosome codes, treating their variants as unplaced], --allow-no-sex Do not force ambiguous-sex phenotypes to missing, --pheno phenotype file, --chr-set 40 changes the chromosome set. The first parameter specifies the number of diploid autosome pairs if positive or haploid chromosomes if negative. Given diploid autosomes the remaining modifiers let you indicate the absence of specific non-autosomal chromosomes. When there are n autosome pairs the X chromosome is assigned numeric code n+1, Y is n+2, XY (pseudo-autosomal region of X) is n+3, and MT (mitochondria) is n+4,  --pca extract principal components, --out output

# Visualizing GWAS in R

library("ggplot2")

gwas <- read.table("file.assoc",header=T)

gwas <- na.omit(gwas)

ggplot(gwas, aes(x= CHR,y= P))+
  xlab("Contig") + ylab("P") + geom_point() + scale_x_log10()
  
##or

plot(gwas.HAAM, x= CHR, y= P, -log10(x))

##CHR= chromosome and P= p-value

# Visualizing PCA in R 

pca <-read.table("file.eigenvec",sep=" ",header=F)

plot(data=pca, column_3~column_4)

pca$Phenotype <- c(I created a new column listing out the phenotypes correlated with the number codes [0,1,2])

library("ggplot2")

ggplot(pca, aes(x= column_3,y= column_4, color= Phenotype))  +
  xlab("PC1") + ylab("PC2") + geom_point(shape=1, size= 3)