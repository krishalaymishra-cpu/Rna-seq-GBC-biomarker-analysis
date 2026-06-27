# DEG-Analysis-pipeline
RNA-seq Feature Selection and TWAS Validation (GSE138109 / SRP223526)

Overview
--------
This repository documents a complete RNA-seq analysis pipeline starting from raw GEO/SRA datasets (GSE138109 and SRP223526) through alignment, quantification, differential expression analysis, feature selection, and validation using TWAS Hub and TWAS Atlas.

The workflow integrates HISAT2, featureCounts, and R-based analysis (DESeq2 + Elastic Net) to identify high-confidence gene candidates supported by both transcriptomic and genetic association evidence.

Workflow Steps
--------------
1. Data Extraction
   - Download raw FASTQ files from GEO/SRA (GSE138109 / SRP223526).
   - Convert SRR IDs to GSM IDs for sample tracking.

2. Alignment
   - Build HISAT2 genome index using Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.
   - Align paired-end reads with HISAT2.
   - Generate sorted and indexed BAM files using samtools.

3. Quantification
   - Use featureCounts to generate raw count matrices across all samples.
   - Produce a unified counts.txt file for downstream analysis.

4. Differential Expression Analysis (DEG)
   - Conduct DEG analysis in R using DESeq2.
   - Identify significantly upregulated and downregulated genes.

5. Feature Selection
   - Apply Elastic Net regression for feature selection.
   - Select candidate genes separately for upregulated and downregulated sets.

6. TWAS Validation
   - Query feature-selected genes against TWAS Hub (GTEx-based reference weights).
   - Validate common DEGs using TWAS Atlas (aggregated TWAS results).
   - Retain genes supported by both transcriptomic and genetic association evidence.

Outputs
-------
- Feature count matrix (final_feature_count.txt)
- DEG results (upregulated and downregulated gene lists)
- Elastic Net-selected features
- TWAS-validated candidate genes (Hub + Atlas)

Tools & Packages
----------------
- Alignment: HISAT2, samtools
- Quantification: featureCounts
- DEG Analysis: DESeq2 ,edgeR and limma-voom (R)
- Feature Selection: glmnet (Elastic Net)
- Validation: TWAS Hub, TWAS Atlas

How to Reproduce
----------------
1. Build HISAT2 index:
   hisat2-build Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa genome_index

2. Align reads:
   hisat2 -x genome_index -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz | samtools sort -o sample.sorted.bam
   samtools index sample.sorted.bam

3. Generate counts:
   featureCounts -T 8 -p -a annotation.gtf -o counts.txt *.sorted.bam

4. Run DEG + Elastic Net in R (see scripts in /R_scripts).

5. Validate selected genes using TWAS Hub & TWAS Atlas.

Citation
--------
- HISAT2: Kim et al., Nature Methods (2015)
- featureCounts: Liao et al., Bioinformatics (2014)
- DESeq2: Love et al., Genome Biology (2014)
- TWAS Hub / Atlas: Gusev et al., Nature Genetics (2016); Zhang et al., Nucleic Acids Research (2020)
