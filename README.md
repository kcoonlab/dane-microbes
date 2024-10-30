## Overview 
**Raw data and scripts for:**
Historical mosquito colonization dynamics drive patterns of microbial community assembly in aboveground aquatic habitats

## Authors 
* Serena Y. Zhao
* Kerri L. Coon - kerri.coon@wisc.edu

## Analysis overview 
Scripts for each analysis are written in R. Each directory contains necessary files and code to recreate each figure, table, and associated statistical analyses reported in the manuscript. To repeat the analysis, clone the repository, and then run each script. Do not `cd` into the cloned repository. 

**Example**
* Create a new project in RStudio. To run the script to recreate Fig. 1 in the manuscript: 
	* Navigate to the terminal window and `git clone https://github.com/kcoonlab/catch-basin-microbes`.
	* Open the script `Fig1.R` from the files panel window.
	* Install required packages. 
	* `cmd enter` from line `1`.

**Before getting started**
* Run the script `phyloseq-object.R` (code from line `1` to `42`) to generate a phyloseq object from the appropriate qiime artifacts. This object will be referenced in some downstream analyses.
* **Note**: this script also contains the code used to estimate alpha diversity (ASV richness and Shannon's H index) in each sample (see lines `46` to `52`)

## Recreate the manuscript figures, tables, and associated statistical analyses
Once the repository has been cloned (above), recreate each figure/table/data analysis as follows: 

**Fig. 1: Microbiota diversity in aquatic habitats with variable historical mosquito productivity**
* Fig 1A - Script: `Fig1.R`: run code from line `1` to `17`
* Fig 1B - Script: `Fig1.R`: run code from line `1` to `29`
* Fig 1C - Script: `Fig1.R`: run code from line `1` to `29`
* Fig 1D - Script: `Fig1.R`: run code from line `1` to `29`

**Fig. 2: Community assembly processes in habitat microbiota biotypes**
* Fig 2A - Script: `Fig2.R`: run code from line `1` to `29`
* Fig 2B - Script: `Fig2.R`: run code from line `1` to `50`

**Table 1: Partitions of variation in βNTI accounted for by mosquito productivity measures and spatial factors**
* Table 1 - Script: `Table1.R`: run code from line `1` to `53`

#### Supplementary information 
* Fig. S1A - Script: `supp-info.R`: run code from line `1` to `31`
* Fig. S1B - Script: `supp-info.R`: run code from line `1` to `31`
* Fig. S1C - Script: `supp-info.R`: run code from line `1` to `31`
* Fig. S2 - Script: `supp-info.R`: run code from line `1` to `61`
* Fig. S3A - Script: `supp-info.R`: run code from line `1` to `109`
* Fig. S3B - Script: `supp-info.R`: run code from line `1` to `109`
* Fig. S3C - Script: `supp-info.R`: run code from line `1` to `109`
* Fig. S3D - Script: `supp-info.R`: run code from line `1` to `109`
* Fig. S4A - Script: `supp-info.R`: run code from line `1` to `128`
* Fig. S4B - Script: `supp-info.R`: run code from line `1` to `143`
* Fig. S5 - Script: `supp-info.R`: run code from line `1` to `185`

## Citation 
Zhao SY, Hausbeck J, Coon KL. Historical mosquito colonization dynamics drive patterns of microbial community assembly in aboveground aquatic habitats.

