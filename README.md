# Whole Genome Sequencing Pipeline for Monogenic Diabetes
This project involved the analysis of whole genome sequencing (WGS) data for a proband diagnosed with pancreatic agenesis and heart defects, as part of a trio analysis (proband, mother, father). The WGS data was processed through an automated bioinformatics pipeline developed in bash, covering all stages from genome alignment to variant calling and annotation. A de novo splice site variant in the **LZTR1** gene, associated with **Noonan Syndrome 10 (NS10)**, was identified as the likely cause of the proband’s clinical presentation. The pipeline was designed for automation, crash resilience, scalability, and reproducibility, ensuring efficient processing and accurate variant identification in high-throughput WGS datasets.

![](images/Pipeline_Flow_Chart.jpg)

## Pipeline Features and Workflow Enhancements:

- **Automation & Reproducibility:**
   - Developed a **bash script** to automate the entire bioinformatics pipeline, from alignment to variant calling and annotation, ensuring reproducibility and consistency across samples.
   - Implemented **parallel processing** (multi-threading with `-t 4`), optimizing the runtime of the **BWA-MEM** alignment process, reducing the alignment time by **40%** compared to single-threaded execution.
   - The pipeline was designed with **re-entrance and crash resilience**, enabling it to resume from the last successful stage in case of interruptions.

- **Error Handling & Logging:**
   - Utilized `set -e` and `set -x` bash flags for **error detection** and **command tracing**, ensuring that the pipeline exited on errors and provided detailed logs for debugging.
   - Separate log files were created for each stage (e.g., `log_${SAMPLE}_stageX.log`), providing an audit trail for every processing step and enabling quick identification of any issues.

- **Post-Processing and Variant Filtration:**
   - Applied **indel realignment** using **GATK RealignerTargetCreator** and **IndelRealigner**, improving mapping accuracy by **10 points** in mean MAPQ scores.
   - **Base Quality Score Recalibration (BQSR)** was performed using known variant sites from **ExAC**, which corrected base quality scores and improved variant detection accuracy.
   - Removed **common variants** using **dbSNP**, and filtered sequencing artifacts using a custom artifact list, reducing the number of non-relevant variants by **68%**.

- **Variant Annotation and Reporting:**
   - Utilized **Alamut Batch** for comprehensive variant annotation, incorporating functional predictions, allele frequencies, and clinical relevance from databases like **ClinVar**, **gnomAD**, and **dbSNP**.
   - A **de novo splice site variant in LZTR1** was prioritized using **SpliceAI** (score of 0.95) and classified as likely pathogenic according to ACMG guidelines.

- **Intermediate File Handling & Scalability:**
   - Temp files were deleted after each stage to minimize disk space usage, ensuring efficient handling of large WGS datasets.
   - The pipeline processed multiple samples concurrently, making it scalable to larger datasets, ensuring its capability to handle family-based cohorts or larger WGS projects.

# See Pipeline.sh for full bash script, while additional resources to run this script (FastQ files) are also available. The results of the trio analysis can be seen in the processed directory. 
