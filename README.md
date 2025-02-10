# Predicting the bacterial host range of plasmid genomes using the language model-based one-class SVM algorithm

## Change logs
2025.02.10 The README.md was updated with a few modification.

## Introduction
The prediction of the plasmid host range is crucial for investigating the dissemination of plasmids and the transfer of resistance and virulence genes mediated by plasmids.

Several machine learning-based tools have been developed to predict plasmid host ranges. These tools have been trained and tested based on the bacterial host records of plasmids in related databases. Typically, a plasmid genome in databases such as NCBI is annotated with only one or a few bacterial hosts, which does not encompass all possible hosts.

Consequently, existing methods may significantly underestimate the host ranges of mobilizable plasmids. In this work, we propose a novel method named HRPredict, which employs a word vector model to digitally represent the encoded proteins on plasmid genomes.

## Citations
Feng T, Chen X, Wu S, Zhou H and Fang Z. "Predicting the bacterial host range of plasmid genomes using the language model-based one-class SVM algorithm."

## Requires
+ biopython >= 1.78
+ prokka >= 1.14.6
+ python >= 3.9.19
+ r-base >= 4.3.3

## Dependencies
+ bioconductor-biocparallel >= 1.36.0
+ jellyfish >= 2.2.10
+ r-caret >= 6.0_94
+ r-dplyr >= 1.1.0
+ r-e1071 >= 1.7_13
+ r-rocr >= 1.0_11
+ r-stringr >= 1.5.1
+ pandas >= 2.2.3

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
# create the environment using configured yaml file

% conda env create -f environment.yaml
```

```
# Manually create the conda environment

% conda create -n HRPredict
% conda activate HRPredict
% conda install conda-forge::r-base
% conda install bioconda::bioconductor-biocparallel
% conda install conda-forge::r-caret
% conda install conda-forge::r-dplyr
% conda install conda-forge::r-e1071
% conda install conda-forge::r-stringr
% conda install r::r-rocr
% conda install bioconda::prokka
% conda install anaconda::python
% conda install biopython
% conda install bioconda::jellyfish
% conda install anaconda::pandas
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

## Output file
| File | Description |
| ------------ | ------------ |
| HostRange_Family.tsv | This file contians the predicted host range of each plasmid at family level |
| HostRange_Genus.tsv | This file contians the predicted host range of each plasmid at genus level |
| HostRange_Species.tsv | This file contians the predicted host range of each plasmid at species level |

## Output file format
| field | Description |
| :---------: | :---------: | 
| id | Plasmid id |
| host_range | Predicted host ranges of each plasmid at different level |

## Output example
| ID | Predicted host range with decision value (SVM) |
| :---------: | :---------: | 
| NZ_CP122477.1 | f__Enterobacteriaceae(1.01);g__Enterobacter(0.03);s__Enterobacter_hormaechei(-0.01)	f__Enterobacteriaceae(1.01);g__Salmonella(0.06);s__Salmonella_enterica(0.07)	f__Enterobacteriaceae(1.01);g__Shigella(0.02);s__Shigella_flexneri(0.01) |
| NC_005970.1 | No species was predicted |

## Contact
Tao Feng - fengtaosmu@foxmail.com

## License

HRPredict is distributed under a GPL-3.0 license.
