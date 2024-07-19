# Unique genetic and risk-factor profiles in clusters of major depressive disorder-related multimorbidity trajectories [from a study of 1.2 million subjects]

This repository houses all the software and R scripts utilized in the clustering pipeline as detailed in our study. Within, users will find the resources and instructions needed to replicate our clustering methodology and conduct analyses using our approach.

## Background
The heterogeneity and complexity of symptom presentation, comorbidities and genetic factors pose challenges to the identification of biological mechanisms underlying complex diseases. Current approaches used to identify biological subtypes of major depressive disorder (MDD) mainly focus on clinical characteristics that cannot be linked to specific biological models. In our study, we examined multimorbidities to identify MDD subtypes with distinct genetic and environmental factors. We leveraged dynamic Bayesian network approaches to determine a minimal set of multimorbidities relevant to MDD and identified seven clusters of disease-burden trajectories throughout the lifespan among 1.2 million participants from transnational cohorts. 

## Citation
Gabriella Juhasz, Andras Gezsi, Sandra Van der Auwera et al. Unique genetic and risk-factor profiles in multimorbidity clusters of depression-related disease trajectories from a study of 1.2 million subjects, 01 August 2023, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-3199113/v1]

## System Requirements

### R Environment

The R scripts in this repository are designed to run within an R environment. Key requirements include:

- **R Version**: 3.5.0 or higher. We recommend using the latest version for optimal performance and compatibility.
- **Operating System Compatibility**: The R scripts are compatible with Windows, macOS, and Linux. However, the full pipeline, including the in-house developed software `bnmcmc`, is designed to operate exclusively under Linux.

### Additional Software

- **bnmcmc**: Our clustering pipeline incorporates the `bnmcmc` software that performs Bayesian Direct Multimorbidity Mapping (BDMM) calculations, working only on Unix/Linux systems. The software's necessary libraries are statically linked, ensuring ease of use across Unix/Linux environments without the need for additional installations.

### R Packages

- Each R script is self-contained and will automatically install all necessary R packages upon execution. This ensures that all dependencies are satisfied without manual package installation.

## Installation

1. To set up your R environment, please download and install R from [CRAN](https://cran.r-project.org/).
2. Clone this repository to your local machine with `git clone`, followed by the repository URL.
3. Ensure you are operating on a Linux system if intending to use the full pipeline with `bnmcmc`.

## Usage

- **R Scripts**: To run any R script, navigate to the directory containing the script and execute it in the terminal. Each script includes detailed usage instructions accessible by running the script with the `--help` argument, e.g., `Rscript script_name.R --help`.
- **bnmcmc Software**: When executing the full pipeline, the appropriate R script will create a parameter file and a bash script that can be used to run this software.

## Compute MDD-related Consensual Clusters on Your Own Dataset
For users interested in calculating the MDD-related clusters as outlined in our study without running the entire pipeline, you can simply use the `compute_consensual_clusters.R` script (which can be found within the `consensual_clusters` directory of this repository) with your dataset, which should include disease onset information. For guidance on the required format of the input dataset file, please refer to the "Demo Dataset" section below.

### Outputs
Running the `compute_consensual_clusters.R` script generates three files in the specified output directory, each containing key data from the clustering process:

- **Cluster Indices (`cluster_indices.txt`)**: This file lists the cluster assignments for each individual in the dataset, providing a straightforward way to identify which cluster each subject belongs to.

- **Cluster Probabilities (`cluster_probabilities.txt`)**: Contains the probability of each individual's assignment to their respective cluster, offering insights into the confidence of each clustering decision.

- **Cluster Probabilities (Log Odds) (`cluster_probabilities_logodds.txt`)**: Contains the posterior log odds of each individual's assignment to their respective cluster, offering insights into the confidence of each clustering decision. These output variables were used in our study as dependent variables in regression analyses and genome-wide association calculations.

- **Cluster Variables (`cluster_variables.txt`)**: Details the variables that were considered in defining the clusters, enriching the context for interpreting the clustering results.

### How to Run

Open your terminal and execute the following command for usage instructions:

```sh
Rscript compute_consensual_clusters.R --help
```

### Usage Instructions

```plaintext
Usage: compute_consensual_clusters.R [options]

Options:
  -i INPUT, --input=INPUT
      Input file name (csv).

  -o OUTPUT, --output=OUTPUT
      Output directory name [default: . (current directory)].

  --cluster-definitions=CLUSTER-DEFINITIONS
      Cluster definition file name [default: NA].

  --id=ID
      Name of the ID column (if set, it will be copied to the result files to identify subjects) [default: NA].

  --age=AGE
      Name of the variable 'age' [default: Age].

  --onset-intervals=ONSET-INTERVALS
      Onset interval boundaries in the format of comma-separated numerical values (age in years) without spaces
      [default: 20,40,60,70].

  -h, --help
      Show this help message and exit.
```
### Required Arguments
- `-i` argument specifying the disease onset dataset. Ensure your dataset is prepared according to the specified requirements.
- `--cluster-definitions` argument specifying the cluster definitions. For the MDD-related clusters derived in our study, the `cluster-definitions-all-7.xlsx` file should be used.

### Example Usage on Demo Data

Execute the following command within the `consensual_clusters` folder:

```sh
Rscript compute_consensual_clusters.R -i ../data/demo_dataset.csv --cluster-definitions=cluster-definitions-all-7.xlsx
```

### Expected Run Time
The execution time for the `compute_consensual_clusters.R` script can vary from 1 to 15 minutes, primarily dependent on the size of the dataset being analyzed.

# Full Pipeline
## Demo Dataset

A demo dataset is included in this repository that serves as an example of the input data format required by our clustering pipeline. The dataset features basic individual descriptors including Age (real number, in years), Sex (integer), and Income (integer), alongside disease onset information denoted by the age at first diagnosis for a range of diseases encoded with ICD-10 codes (as the remaining columns of the table, real numbers, in years). Missing items (indicating a diagnosis is not present) should be coded with `NA`. Importantly, our pipeline is designed with flexibility in mind and can accommodate disease onset information coded in formats other than ICD-10.

### Accessing the Demo Dataset

The demo dataset can be found within the `data` directory of this repository. To begin exploring this example, navigate to the directory, and you will find the dataset file named `demo_dataset.csv`. This file is structured to provide a clear representation of the type of data our clustering pipeline analyzes, demonstrating the format and type of information required.

#### Sample Snippet from Demo Dataset:

```csv
Age,Sex,Income,A01,A02,A03,A04,A05,A06,A07,A08,A09,A15,A16,...
73,1,1,NA,NA,NA,NA,NA,NA,NA,NA,NA,19.22,NA
70,1,1,NA,NA,NA,NA,NA,NA,NA,NA,57.06,NA,NA
63,1,2,NA,NA,NA,NA,NA,NA,NA,NA,40.38,NA,NA
54,1,2,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA
49,2,2,NA,NA,7.2,NA,NA,NA,NA,NA,NA,NA,NA
74,1,2,NA,NA,NA,73.3,NA,NA,NA,NA,NA,NA,NA
66,2,1,NA,NA,NA,NA,NA,NA,NA,NA,NA,0,NA
72,1,2,NA,NA,NA,53.63,NA,NA,NA,NA,NA,NA,NA
...
```

*Note: This is a simplified example. Please refer to the actual `demo_dataset.csv` file for a comprehensive view.*

To utilize this dataset with our clustering pipeline, simply follow the instructions provided in the `Usage` section, substituting your data paths with those leading to the `demo_dataset.csv` file within the `data` directory.

## Step 1: Dataset Check with `check_dataset.R`, collecting descriptive statistics, creating exploratory plots

The initial step in our full pipeline involves running the `pipeline/step1/check_dataset.R` script on your dataset. This script conducts a basic data sanity check to ensure the dataset is properly formatted and complete for further analysis.

### What it Checks

- Verifies the presence of essential variables: Age, Sex, and Household Income.
- Identifies missing values in Age and Sex variables.
- Searches for any negative values.
- Ensures Age is not lower than any disease onset.

### Outputs

Upon successful execution, the script:

- Computes descriptive statistics for disease onsets and their prevalences, saving the results in `disease-prevalences.xlsx`.
- Calculates pairwise associations between diseases, with findings documented in `disease-comorbidities.xlsx`.
- Generates several plots as PDF files, visualizing key data aspects.

### How to Run

Open your terminal and execute the following command for usage instructions in the `pipeline/step1` folder:

```sh
Rscript check_dataset.R --help
```

### Usage Instructions

```plaintext
Usage: check_dataset.R [options]

Options:
  -i INPUT, --input=INPUT
      Input file name (csv).

  --age=AGE
      Name of the variable 'age' [default: Age].

  --sex=SEX
      Name of the variable 'sex' [default: Sex].

  --household-income=HOUSEHOLD-INCOME
      Name of the variable 'household income' [default: Income].

  --depression=DEPRESSION
      Name of the variable 'depression' (for multiple variables, format as comma-separated strings without spaces) [default: F32,F33].

	--codetable=CODETABLE
		  Filename for code table containing disease names for each disease code (columns of the input file) [default: NA]

  -r REMOVE, --remove=REMOVE
      Names of variables to remove (format: comma-separated strings without spaces).

  --plot
      Create plots [default: TRUE].

  --no-plot
      Don't create plots.

  -h, --help
      Show this help message and exit.
```

### Expected Run Time
The execution time for the `check_dataset.R` script can last several hours, which can be reduced to several minutes by disabling the creation of plots (`--no-plot` argument).

### Optional Argument

The `--codetable` argument can be used to assign disease descriptions to the output prevalence tables. For ICD-10 codes, the provided `data/disease_ICD_10.xlsx` file can be used. Check this file if you have specific needs.

### Example Usage on Demo Data

Execute the following command within the `pipeline/step1` folder:

```sh
Rscript check_dataset.R -i ../../data/demo_dataset.csv --codetable=../../data/disease_ICD_10.xlsx
```

## Step 2: Data Transformation and Calculation of Relevance Scores Using Inhomogeneous Dynamic Bayesian Networks
This step is divided into two parts: the first for transforming raw onset data and the second for performing BDMM computations.

### Part 1: Data Transformation

Transform your raw onset data into a format suitable for BDMM analysis (utilizing inhomogeneous dynamic Bayesian networks, ihDBN) using `pipeline/step2/transform_onset_to_inhomogenDBN.R`.

```sh
Rscript transform_onset_to_inhomogenDBN.R -i disease_onsets.csv
```

Open your terminal and execute the following command for detailed usage instructions in the `pipeline/step2` folder:

```sh
Rscript transform_onset_to_inhomogenDBN.R --help
```

**Outputs:**
- Time slice datasets (`dataset_ihDBN_interval-*.csv`).
- Auxiliary files, including model files (`*.xml`).
- A bash script file (`run_BDMM.sh`) for executing BDMM computations.

**Expected Run Time**
The execution time for the `transform_onset_to_inhomogenDBN.R` script takes several minutes.

**Example Usage on Demo Data**
Execute the following command within the `pipeline/step2` folder:

```sh
Rscript transform_onset_to_inhomogenDBN.R -i ../../data/demo_dataset.csv
```

### Part 2: BDMM Computations

Execute the BDMM computations for each time slice dataset by running the bash script created in Part 1. This process is designed to run in parallel, significantly optimizing computation time.

```sh
bash run_BDMM.sh
```

**Note:**
Please ensure to grant executable permissions to the `bnmcmc.exe` file.

**Outputs:**
- Log files documenting the computation process.
- Output files for each time slice (`bdmm_output_*`).

**Expected Run Time**
The execution time for the `run_BDMM.sh` script can last several hours.

## Step 3: Patient Clustering, Creating Exploratory Plots

The final step involves executing `pipeline/step3/compute_clusters.R` to perform the actual patient clustering based on the processed data and BDMM analysis results.

## Outputs

Upon successful execution of `compute_clusters.R`, the script generates the following outputs, where `N` denotes the number of clusters specified:

- **Cluster Indices (`cluster_indices_N.txt`)**: Lists the cluster assignments for each subject, indicating which cluster each individual has been grouped into.

- **Cluster Probabilities (`cluster_probabilities_N.txt`)**: Provides the probability scores for each subject's assignment to their respective cluster, highlighting the confidence level of each clustering decision.

- **Cluster Variables (`cluster_variables_N.txt`)**: Details the variables considered in the formation of clusters, offering insights into the factors driving the clustering process.

- **Plots**: Includes several visual representations of the clustering analysis, facilitating an intuitive understanding of the data and clustering outcomes.

- **Excel Table for Disease Over/Underrepresentations**: An Excel table summarizing diseases that are overrepresented or underrepresented within each cluster, providing valuable insights into the disease composition and patterns across different clusters.

### How to Run

Open your terminal and execute the following command for usage instructions in the `pipeline/step3` folder::

```sh
Rscript compute_clusters.R --help
```

### Usage Instructions

```plaintext
Usage: compute_clusters.R [options]

Options:
  -i INPUT, --input=INPUT
      Input file name (csv).

  -o OUTPUT, --output=OUTPUT
      Output directory name prefix (will be appended with the number of clusters) [default: clusters].

  --cohort-profile=COHORT-PROFILE
      Cohort-specific comorbidity profile definition based on BDMM analysis results, name of directory [default: NA].

  --number-of-clusters=NUMBER-OF-CLUSTERS
      Number of clusters to generate [default: NA].

  --age=AGE
      Name of the variable 'age' [default: Age].

  --sex=SEX
      Name of the variable 'sex' [default: Sex].

  --household-income=HOUSEHOLD-INCOME
      Name of the variable 'household income' [default: Income].

  --depression=DEPRESSION
      Name of the variable 'depression' (for multiple variables, format as comma-separated strings without spaces) [default: F32,F33].

	--codetable=CODETABLE
		  Filename for code table containing disease names for each disease code (columns of the input file) [default: NA]

  --id=ID
      Name of ID column (will be included in result files to identify subjects) [default: NA].

  --onset-intervals=ONSET-INTERVALS
      Onset interval boundaries (format: comma-separated numerical values [of age in years] without spaces) [default: 20,40,60,70].
```

### Required Arguments
- `-i` argument specifying the disease onset dataset. Ensure your dataset is prepared according to the specified requirements. This file should be the one that was used in the earlier steps.
- `--cohort-profile` argument specifying the comorbidity profile definition based on BDMM analysis results. This could be the directory with the outputs of the BDMM analysis.
- `--number-of-clusters` argument specifying the desired number of clusters.

### Expected Run Time
The execution time for the `compute_clusters.R` script takes several minutes.

### Example Usage on Demo Data
Execute the following command within the `pipeline/step3` folder:

```sh
Rscript compute_clusters.R -i ../../data/demo_dataset.csv --cohort-profile=../step2 --number-of-clusters=5
```
