# Decoding-the-impacts-of-Sodium-Propionate-on-A549-Cell-line-An-Analysis-from-RNA-Seq-to-ML

# README for Gene Expression Data Analysis

## Introduction
This project involves comprehensive analysis of gene expression data, including clustering, correlation analysis, regression modeling, and visualization. The dataset contains numerical and categorical columns representing gene expression levels under various conditions. The analyses were performed using Python and R, leveraging multiple machine learning techniques and statistical models.

Name: [Jaya sruthi koppada]
Programming language: UNIX , R
Description: This script read text file.

Required files:
The given data from (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA934530Links to an external site.).  
Required Packages: fastqc, star, subread

Step:1 Open Command Prompt, logging in using carbonate account
       sshsobselva@carbonate.uits.iu.edu

     Enter the Password and press 1 to get confirmation code in Duo Push
     Logged in to the Carbonate account,
    
Command: pwd
Output: N/u/*****

# Changing to scratch account
Command: /N/scratch/*****

Step: 2 Open the given data link and move on to ENA browser.
       and copying the link of the sample 

# Downloading the data
Command:  
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR233/066/SRR23362166/SRR23362166_1.fastq.gz
Output: SRR23362166_1.fastq.gz 
(downloaded)
(Using wget I have downloaded all thesamples)

Step:3

# making a directory to the datafiles (creating a separate folder for datafiles)
Command: mkdir datafiles

# Copy the datafiles 
(From /N/scratch/sobselva/ to /N/scratch/*****/datafiles)
Command: mv SRR2336* /N/scratch/*****/datafiles
Output: 
[*****sobselva@h1 datafiles]$ ls
SRR23362162_1.fastq.gz                           SRR23362168_1.fastq.gz  SRR23362174_1.fastq.gz  SRR23362180_1.fastq.gz  SRR23362186_1.fastq.gz
SRR23362162_2.fastq.gz                           SRR23362168_2.fastq.gz  SRR23362174_2.fastq.gz  SRR23362180_2.fastq.gz  SRR23362186_2.fastq.gz
SRR23362163_1.fastq.gz                           SRR23362169_1.fastq.gz  SRR23362175_1.fastq.gz  SRR23362181_1.fastq.gz  SRR23362187_1.fastq.gz
......................................(like these all files has been downloaded)

Step: 3 Checking the quality (fastqc was used) (for doing the quality control check, checking the module and loading it)

module avail
module load  fastqc/0.11.9

# creating directories for quality checking files' output
Command: mkdir qualitycheck

Command: SRR23362166_1.fastq.gz  -o qualitycheck 

Output: SRR23362166_1_fastqc.html
        SRR23362166_1_fastqc.zip
Checking the .html file by copying the file to the local disk
Command:  scp *****@carbonate.uits.iu.edu:/N/scratch/*****/datafiles/qualitycheck/SRR23362166_1_fastqc.html .

(Using the above code, we have done the quality control check for all the samples)

Step:4 # Aligning a data to the human reference genome (star was used)
# Performing gene indexing 
( Using the command wget downloading the  human genome gencode.v38.chr_patch_hapl_scaff.annotation.gtf.gz , GRCh38.p13.genome.fa.gz)

# Unzipping these two files for performing downstream analysis.
Command: gunzip gencode.v38.chr_patch_hapl_scaff.annotation.gtf.gz
         gunzip GRCh38.p13.genome.fa.gz

Output: gencode.v38.chr_patch_hapl_scaff.annotation.gtf.gz
        GRCh38.p13.genome.fa (These are the output of unzipped file)

# Creating a file for indexing script (using command nano)
Command: nano indexing.sh

#!/bin/bash
#SBATCH --job-name=gene_indexing
#SBATCH --output=gene_indexing.out
#SBATCH --error=gene_indexing.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=78gb
#SBATCH --time=05:00:00
#SBATCH -A students

# Loading necessary module (using the tool star to align the reads to the human reference genome) 
module load star/2.7.3a

# Creating a pathway 
Command: indexing="/N/scratch/*****/"

# Creating a new directory for indexing
Command: mkdir -p $indexing/myindexing

# Running star command

STAR --runThreadN 12 \
--runMode genomeGenerate \
--genomeDir $indexing/myindexing \
--genomeFastaFiles $indexing/GRCh38.p13.genome.fa \
--sjdbGTFfile $indexing/gencode.v38.chr_patch_hapl_scaff.annotation.gtf \
--sjdbOverhang 124

(This is the bash script I have used for gene indexing)


# Submitting this bash script  (submitting the bash script using the command sbatch)
 Command: sbatch indexing.sh

# Checking the job run time (To check whether my job is running or not, using the command squeue)
 Command: squeue -u *****

# Checking the start time of the job
Comamnd: squeue -u ***** --start

( Once job done, I could see that in gene_indexing.out, finished successfully)

# Checking the indexing directory
Command: cd /N/scratch/*****/myindexing
Output: 
[*****@h1 assignment2]$ cd myindexing/
[*****@h1 myindexing]$ ls
chrLength.txt      chrStart.txt      geneInfo.tab          SA            sjdbList.fromGTF.out.tab
chrNameLength.txt  exonGeTrInfo.tab  Genome                SAindex       sjdbList.out.tab
chrName.txt        exonInfo.tab      genomeParameters.txt  sjdbInfo.txt  transcriptInfo.tab


Step: 5 # Performing gene mapping (star was used)

# Creating a new file for mapping script
Command: nano mapping.sh


#!/bin/bash
#SBATCH --job-name=new_mapping
#SBATCH --output=new_mapping.out
#SBATCH --error=new_mapping.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=*****@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=78gb # Increase memory allocation
#SBATCH --time=11:59:00
#SBATCH -A students

# Loading star module (from module avail command , I have loaded the star)
module load star/2.7.3a

# Create new directory for output files ( giving path of the newly created directory)
OUTFILE="/N/scratch/*****/my_mapping/"

# Path for the genome directory (giving path of the myindexing directory)
GENOMEDIR="/N/scratch/*****/myindexing/"

# Path to the datafiles(fastq files) (This is the path of the downloaded raw RNA seq data)
FASTQ_FILES="/N/scratch/*****/datafiles/"

# Using for loop for the sample files to map to the reference genome at a time)
for ((i=62; i<=88; i++)); do
    SAMPLE=" SRR233621${i}"
    FILE1="${FASTQ_FILES}/${SAMPLE}_1.fastq.gz"
    FILE2="${FASTQ_FILES}/${SAMPLE}_2.fastq.gz"

# Running the star command for mapping
 echo "Executing STAR command: STAR --runThreadN 12 --genomeDir $GENOMEDIR --readFilesIn $FILE1 $FILE2 --readFilesCommand zcat --quantMode GeneCounts --outFileNamePrefix /N/scratch/*****/assignment2/my_mapping/${SAMPLE}"  
  STAR --runThreadN 12 --genomeDir $GENOMEDIR --readFilesIn $FILE1 $FILE2 \
    --readFilesCommand zcat \
    --quantMode GeneCounts   --outFileNamePrefix /N/scratch/*****/my_mapping/${SAMPLE} \
    --outSAMmapqUnique 60 --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx
done


( I have given the path to the output file, I gave a code for creating output for all file in the separate folder )


# Submitting the script (Submitting the bash bash script for aligning the reads)
Command: sbatch mapping.sh

# Checking the job 
Command: squeue -u *****


# Changing to the mapping output directory 
Command: cd /N/scratch/*****/my_mapping/
Output:
 SRR23362167  SRR23362173  SRR23362179  SRR23362185 SRR23362162 SRR23362163           
 ......................................................................( Like these all samples reads will be stored in separate folder)
                                      
# Checking this directory

Command: cd SRR23362162
Output: [*****@h2 my_mapping]$ cd   SRR23362162
[*****@h2 SRR23362162]$ ls
feature_counts.csv          SRR23362162_Aligned.sortedByCoord.out.bam                    SRR23362162_Log.out               SRR23362162_SJ.out.tab
feature_counts.txt          SRR23362162_Aligned.sortedByCoord.out.bam.featureCounts.bam  SRR23362162_Log.progress.out      SRR23362162_Unmapped.out.mate1
feature_counts.txt.summary  SRR23362162_Log.final.out                                    SRR23362162_ReadsPerGene.out.tab  SRR23362162_Unmapped.out.mate2
[
(Using the bam file generated, I have carried out the further analysis)


Step: 6 # Performing quantification (subread was used)

module avail
module load subread/2.0.6

Command: featureCounts  -T 4 -a /N/scratch/*****/gencode.v44.chr_patch_hapl_scaff.annotation.gtf -o counts.txt *.bam

(Output has been seen in counts.txt)

Output: counts.txt
        counts.txt.summary

# Open the counts.txt file

Command: nano counts.txt

Output:
Geneid  Chr     Start   End     Strand  Length  /N/scratch/*****/my_mapping/SRR23362162/SRR23362162_Aligned.sortedByCoord.out.bam (like this all samples number with its pathway available)
ENSG00000223972.5       chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1    11869;12010;12179;12613;12613;12975;13221;13221;13453   12227;12057;12227;12721;12697;1305$ENSG00000227232.5       chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1;chr1  14404;15005;15796;16607;16858;17233;17606;17915;18268;24738;29534       14501;1503$ENSG00000278267.1       chr1    17369   17436   -       68      1
ENSG00000243485.5       chr1;chr1;chr1;chr1;chr1        29554;30267;30564;30976;30976   30039;30667;30667;31109;31097   +;+;+;+;+       1021    0
ENSG00000284332.1       chr1    30366   30503   +       138     0
ENSG00000237613.2       chr1;chr1;chr1;chr1;chr1        34554;35245;35277;35721;35721   35174;35481;35481;36073;36081   -;-;-;-;-       1219    0
..................(like these I have donw feature counts for all the samples)

# Sorting out only the necessary columns (sorting out all the columns of geneId and samples)
Command: cut -f1, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,18,19,20,21,22,23,24,25,26,27 counts.txt > counts_data.txt


# opening new file
Command: nano counts_data.txt
Output: 
.....................................................(Only the geneids and sample counts are extracted)

Step: 7 # Copying this file to the local disk
Command: cd ~/dowloads
          scp sobselva@carbonate.uits.iu.edu:/N/scratch/*****/my_mapping/counts_data.txt .

(Asked for passcode and duo push) 
(Copied to the local machine downloads folder)

Step:8 Performing downstream analysis ( DSEQ2 was used)

# setting the work directory
setwd("C:/Users/mahev/Downloads")
getwd()

# Loading and installing necessary packages and libraries

# Install DESeq2 from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("edgeR")

Command: library(edgeR)


# Loading data (loading the count matrix data in R)
counts_data1 <- as.matrix(read.table("Counts_data.txt", header = TRUE , row.names=1 ))

Command:
head(counts_data1)

Output: 

 Ensembl_ID      Control_3D_R1 Control_3D_R2 Control_3D_R3 SP_3D_R1 SP_3D_R2
  
1 ENSG00000186827             1             1             1        0        0
2 ENSG00000186891             0             0             0        1        3
3 ENSG00000160072            16             9            10       17        8
4 ENSG00000260179             7             4             7        2        6


# Create a DGEList object
y <- DGEList(counts = count_matrix1)

# Determine the number of samples for each condition
n_samples_Control_3D <- sum(grepl("Control_3D", colnames(count_matrix1)))
n_samples_SP_3D <- sum(grepl("SP_3D", colnames(count_matrix1)))
n_samples_SP_12D <- sum(grepl("SP_12D", colnames(count_matrix1)))

# Add sample information directly to the DGEList object
y$group <- factor(c(rep("Control_3D", n_samples_Control_3D), rep("SP_3D", n_samples_SP_3D), rep("SP_12D", n_samples_SP_12D)))

# Perform normalization
y <- calcNormFactors(y)

# Estimate dispersion
y <- estimateDisp(y)

# Assuming following number of samples for each condition
n_samples_Control_3D <- 3
n_samples_SP_3D <- 3
n_samples_SP_12D <- 3

# Add sample information directly to the DGEList object
y$samples$group <- factor(rep(c("Control_3D", "SP_3D", "SP_12D"), each = 3))

# Perform differential expression analysis
fit <- glmFit(y, design = model.matrix(~y$samples$group))

# Likelihood ratio test
lrt <- glmLRT(fit)

# Extract results
results1 <- topTags(lrt)

print(results1)


# Adjust criteria to be less stringent
significant_genes1 <- rownames(results1$table)[results1$table$FDR < 0.05 & abs(results1$table$logFC) > 1.0]

# Print significant genes
print(significant_genes1)

(Like these we have done differential expression analysis of another splitted count matrix of the conditions control 24 hours, sample 24 hours and 3 hours)
############DESeq2
Description: This script answers the code for our project
Required files: FILE1-COUNTS1.txt, FILE2-COLDATA1.csv, FILE3-COUNTS2.csv, and FILE4-COLDATA2.csv

Required packages in R:
1.DESeq2 package
2.ggplot2 package
3.EnhancedVolcano package
4.apeglm package
5.devtools package


Execution:
1)Open Rstudio
2)Import required packages.
3)Extract the required files mentioned in the description to obtain the input and output data for our analysis.
4)After reading the data, performs quality check and normalises data.
6)Run DESeq2 and perform the Visualizations such as MA plots, dispersion plots, Principal Component Analysis.


Output Files: Extracted top 20 significant genes: 'top20_sigs1_norm.csv' and 'top20_sigs2_norm.csv'
Visualizations, including PCA plots, MA plots, and dispersion plots, are saved during script execution.


PYTHON:
*Open Jupyter Notebook
*Execute commands to get outputs
Description:
This project involves data analysis and machine learning tasks implemented in Python. The code includes steps for data cleaning, exploratory data analysis (EDA), machine learning model development, and statistical analysis.
Dependencies:
A list of required Python packages and their versions is provided in the requirements.txt file.
1.pandas
2.numpy
3.matplotlib
4.Seaborn
5.OS
6.sklearn
7.statsmodels
Data Loading and Cleaning:
Data is loaded from a CSV file ('LR-1.csv') using pandas. Missing values and duplicate rows are removed to ensure data quality.

Exploratory Data Analysis (EDA):
The script includes EDA tasks such as generating summary statistics, creating a correlation heatmap, and visualizing clusters using K-Means and DBSCAN.

Machine Learning Model Development:
Two machine learning models are developed: a RandomForestRegressor and a Lasso regression model. The models are evaluated using Mean Squared Error, R-squared, and precision metrics.

Linear Regression Plot:

A linear regression model is created and evaluated, with results visualized through scatter plots using both Seaborn and Matplotlib.

Analysis of Variance (ANOVA):

ANOVA is performed on the dataset to analyze the impact of different categories on a dependent variable, and results are visualized using box plots.

Statistical Analysis with Statsmodels:

The statsmodels library is used to perform ANOVA, and a summary of the results is printed.

Conclusion:
Outputs and visualisations were executed.














## File Structure
- **Raw Data**:
  - `rawcounts.HvsM.txt`: Raw counts data for Healthy vs Mild Cognitive Impairment.
  - `rawcounts.HvsAD.txt`: Raw counts data for Healthy vs Alzheimer's Disease.
  - `coldata.HvsM.txt`: Metadata for Healthy vs Mild Cognitive Impairment.
  - `coldata.HvsAD.txt`: Metadata for Healthy vs Alzheimer's Disease.
- **Scripts**:
  - `analysis.R`: R script for differential expression analysis.
  - `analysis.py`: Python script for clustering, regression, and visualization.
- **Outputs**:
  - `DEGs_HvsM_results.csv`: Differentially expressed genes for Healthy vs MCI.
  - `DEGs_HvsAD_results.csv`: Differentially expressed genes for Healthy vs AD.
  - `H_MCI_volcano.png`: Volcano plot for Healthy vs MCI.
  - `H_AD_volcano.png`: Volcano plot for Healthy vs AD.
  - `H_MCI_AD_logFC.png`: Heatmap of log fold change for selected genes.

## Prerequisites
### Python Dependencies
- Python 3.x
- Required Python packages:
  ```bash
  pip install pandas numpy matplotlib seaborn scikit-learn statsmodels
  ```

### R Dependencies
- R 4.x
- Required R packages:
  ```R
  install.packages(c("DESeq2", "dplyr", "reshape", "ggplot2", "ggrepel", "pheatmap", "RColorBrewer"))
  ```

## Workflow
### Step 1: Data Preparation
1. Load the raw count data and metadata for Healthy vs MCI and Healthy vs AD using R.
2. Perform data normalization and filtering.

### Step 2: Differential Expression Analysis (R)
1. **Healthy vs MCI**:
    ```R
    dds.1 <- DESeqDataSetFromMatrix(countData = cts.1, colData = coldata.1, design = ~ condition)
    dds.1 <- DESeq(dds.1)
    res_HvsM <- results(dds.1, contrast=c("condition","Healthy","Mild_cog"))
    write.csv(res_HvsM, file="DEGs_HvsM_results.csv")
    ```

2. **Healthy vs AD**:
    ```R
    dds.2 <- DESeqDataSetFromMatrix(countData = cts.2, colData = coldata.2, design = ~ condition)
    dds.2 <- DESeq(dds.2)
    res_HvsAD <- results(dds.2, contrast=c("condition","Healthy","Alzheimers"))
    write.csv(res_HvsAD, file="DEGs_HvsAD_results.csv")
    ```

### Step 3: Clustering and Visualization (Python)
1. **Correlation Heatmap**:
    ```python
    import seaborn as sns
    correlation_matrix = file.corr()
    sns.heatmap(correlation_matrix, annot=False, cmap='coolwarm')
    ```

2. **K-Means Clustering**:
    ```python
    from sklearn.cluster import KMeans
    kmeans = KMeans(n_clusters=4)
    kmeans.fit(numeric_file)
    ```

3. **DBSCAN Clustering**:
    ```python
    from sklearn.cluster import DBSCAN
    dbscan = DBSCAN(eps=0.5, min_samples=5)
    cluster_labels = dbscan.fit_predict(numeric_file)
    ```

4. **Volcano Plots**:
    ```R
    ggplot(res.HvsM_df_ordered) +
      geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
      ggtitle("Healthy vs MCI")
    ```

### Step 4: Regression Modeling
1. **Random Forest Regressor**:
    ```python
    from sklearn.ensemble import RandomForestRegressor
    pipeline = Pipeline(steps=[('preprocessor', preprocessor),
                               ('regressor', RandomForestRegressor())])
    pipeline.fit(X_train, y_train)
    ```

2. **LASSO Regression**:
    ```python
    from sklearn.linear_model import Lasso
    lasso = Lasso()
    lasso.fit(X_train, y_train)
    ```

### Step 5: Saving Results
1. Export filtered DEGs and clustering results to CSV files.
2. Save visualizations as PNG files.

## Results
- **Differentially Expressed Genes**:
  - Healthy vs MCI: `DEGs_HvsM_results.csv`
  - Healthy vs AD: `DEGs_HvsAD_results.csv`
- **Clustering**:
  - Heatmap and clustering figures: `H_MCI_AD_logFC.png`
- **Regression**:
  - Random Forest and LASSO regression results: Accuracy and R-squared metrics printed in Python script.

## Author
This analysis was conducted by Jaya Sruthi Koppada. Please feel free to reach out for questions sruthikoppada@gmail.com
.


