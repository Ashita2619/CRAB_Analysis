
# General Package for State Lab for Carbapenem-resistant Acinetobacter baumannii (CRAB) Detection
_______________________________________

## Quick Start Installation

  1. Install necessary conda ENV from [CRAB.yml](resources/CRAB.yml) file in resources 
  > Only needed for HomeBrew pipeline

```
  conda env create -f CRAB.yaml

  conda activate CRAB
```


<br />
_______________________________________

## Before Running update Variables

  1. Update necessary locations in [Pipeline_Variable.json](data/pipeline_variables.json) for your system

<br />
_______________________________________

## Running Analysis

### Running [PHoeNIx](https://github.com/CDCgov/phoenix)
  1. Run [CRAB Runner Script] (scripts/CRAB_Analysis_Runner.py) which takes path to fastq, run date (MMDDYY), pipeline type
  2. Make the third argument "CDC"

```python 
  python CRAB_Analysis_Runner.py PathToFastq RunDate CDC
```
 
 > See [Phoenix Wrapper](docs/WF_Phoenix.md) documentation for more info

### Running HomeBrew Pipeline
  1. Run [CRAB Runner Script](scripts/CRAB_Analysis_Runner.py) which takes path to fastq, run date (MMDDYY)
  2. Do not include a third argument if running Homebrew

```python 
  python CRAB_Analysis_Runner.py PathToFastq RunDate
```

### Running Phylogenetic Tree Builder
  1. Run [CRAB Runner Script](scripts/CRAB_Analysis_Runner.py) which takes path to fastq, rundate (MMDDYY)
  2. Make the third argument "tree"

```python 
  python CRAB_Analysis_Runner.py PathToFastq RunDate tree
```

<br />


_______________________________________

## The package contains the following workflows in their respective subdirectories:

<br />

### **Phoenix Launcher:** [Run CDC Phoenix](docs/WF_Phoenix.md)
 - Wrapper to run CDC's [PHoeNIx](https://github.com/CDCgov/phoenix)
 - Parses and uploads results to internal DB  

<br />
<br />

### **Workflow 0:** [Assembley](docs/WF_0_Assembler.md)
 - Organize sample reads
 - Trim and removed low-quality reads
 - Assembly reads using [SPAdes](https://github.com/ablab/spades)
 - Check assembly quality  

  > This step is required, assembly is needed for downstream analysis.<br>

<br />
<br />

### **Workflow 1:** [Annotate](docs/WF_1_Annotate.md)
 - Find the Multilocus sequence typing (MLST).
 - Annotate genome.  

  > This step is not required but is needed for the final report.<br>

<br />
<br />

### **Workflow 2:** [Find AMR Genes](docs/WF_2_FindAMR.md)
 - Finds antimicrobial resistance genes within the assembled genome
 - Parses output

  > This step is required, to find resistances.
  
<br />
<br />

### **Workflow 3:** [Database Push](docs/WF_3_DB.md)
 - Open HORIZON LIMS database (Oracle).
 - Join all demographics with sample ID.
 - Join demographics with assembly stats.
 - Create resistance gene table.
 - Push new demographics, assembly stats, resistance gene to CRAB DB SQL database.
 - Write demographical information for final result file

  > This step is only required if uploading to a SQL DB.
 
<br />
<br />

### **Workflow 3.5:** [Phylogenetic Tree and SNP Heat Map](docs/WF_3_5_SNP_Phylo.md)
 - Takes consensus sequences and creates a SNP heatmap and Phylogenetic Tree
 
<br />
<br />

### **Workflow 4:** [Create Report](docs/WF_4_final_report.md)
 - Creates final report for Epidemiologists
 - Summarizes all results
  
<br />
<br />



