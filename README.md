# Predicting the bacterial host range of plasmid genomes using the language model-based one-class SVM algorithm

## Introduction
The prediction of the plasmid host range is crucial for investigating the dissemination of plasmids and the transfer of resistance and virulence genes mediated by plasmids.

Several machine learning-based tools have been developed to predict plasmid host ranges. These tools have been trained and tested based on the bacterial host records of plasmids in related databases. Typically, a plasmid genome in databases such as NCBI is annotated with only one or a few bacterial hosts, which does not encompass all possible hosts.

Consequently, existing methods may significantly underestimate the host ranges of mobilizable plasmids. In this work, we propose a novel method named HRPredict, which employs a word vector model to digitally represent the encoded proteins on plasmid genomes.

## Citations
Feng T, Chen X, Wu S, Zhou H and Fang Z. "Predicting the bacterial host range of plasmid genomes using the language model-based one-class SVM algorithm."

## Requires
+ Python >= 3.6.15
+ hmmer >= 3.3.2
+ numpy >= 1.19.5
+ perl >= 5.26.2
+ biopython >= 1.79

## Dependencies
+ jellyfish >= 2.2.10
+ prokka >= 1.14.6
+ r-base >= 4.2.2
+ r-e1071 >= 1.7_13
+ r-dplyr >= 1.1.0

## Preparation
HRPredict is developed using some dependencies, and we recommend using conda to install them.

### Step 1: add channels
```
% conda config --add channels defaults
% conda config --add channels conda-forge
% conda config --add channels bioconda
```

### Step 2: create conda environment

```
% conda env create -f environment.yaml  # create the environment using configured yaml file
```

```
% conda create -n HRPredict # Manually create the conda environment
% conda activate HRPredict
% conda install -c bioconda prokka
% conda install python
% conda install r-base
% conda install r-e1071
% conda install r-dplyr
% conda install biopython
```

## Installation
Clone this repository to your local linux PC.
```
% git clone https://github.com/FengTaoSMU/HRPredict.git
% cd HRPredict
```

## Using HRPredict to perform plasmid host range prediction

```
# input: Plasmid complete fasta file with unique id
% python3 HRPredict.py -i ./test/test.fasta -m ./model/ -r ./reference/ -o .
```

# Output file
| File | Description |
| ------------ | ------------ |
| HostRange_Family.tsv | This file contians the predicted host range of each plasmid at family level |
| HostRange_Genus.tsv | This file contians the predicted host range of each plasmid at genus level |
| HostRange_Species.tsv | This file contians the predicted host range of each plasmid at species level |

# Output file format
| field | Description |
| :---------: | :---------: | 
| id | Plasmid id |
| host_range | Predicted host ranges of each plasmid at different level |

## Contact
Tao Feng - fengtaosmu@foxmail.com

## License

HRPredict is distributed under a GPL-3.0 license.
