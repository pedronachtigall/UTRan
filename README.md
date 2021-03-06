# UTRan
**UTR an**notation of *de novo* assembled transcripts. It is an experimental approach designed to identify UTR sequences in transcripts obtained through *de novo* methods using the coding sequence annotation and the read coverage within transcripts.

It was not thoroughly tested and we give no warranty of confident results, please use it with caution.

## Pipeline
```
|☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬|
|>Parse CDSs from annotation file                                |
|>Map reads into contigs                                         |
|>Hard-filtering of mapped reads                                 |
|>Get coverage data                                              |
|>Identify changing point in 5'UTR                               |
|>Identify changing point and PAS in 3'UTR                       |
|>Output plots and UTR sequences                                 |
|☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬☬|
```
## Requirements

- [Python3](https://www.python.org/)
    - [Biopython](https://biopython.org/wiki/Download)
    - [NumPy](https://numpy.org/)
    - [Pandas](https://pandas.pydata.org/)
    - [Matplotlib](https://matplotlib.org/2.0.2/index.html)
    - [Scipy](https://scipy.org/)
    - [Pathos](https://pypi.org/project/pathos/)
- [BWA](http://bio-bwa.sourceforge.net/)
- [Samtools](http://samtools.sourceforge.net/)
- [Bedtools](https://bedtools.readthedocs.io/)
- [GATK4](https://github.com/broadinstitute/gatk)
- [Picard](https://broadinstitute.github.io/picard/)

### Installation

If the user has all requirements installed an properly working, you just need to do the following steps (please, change ```/.bash_profile``` to ```/.bashrc``` if it is your bash source):
```
git clone https://github.com/pedronachtigall/UTRan.git
echo "export PATH=$PATH:$(pwd)/UTRan/bin/" >> ~/.bash_profile
source ~/.bash_profile
```

Alternatively, the user can install all requirements through [conda](https://conda.io/) manager as follow (please, change ```/.bash_profile``` to ```/.bashrc``` if it is your bash source):
```
conda create -n utran_env python=3.7 biopython bwa samtools bedtools picard pandas matplotlib scipy gatk4 pathos
git clone https://github.com/pedronachtigall/UTRan.git
echo "export PATH=$PATH:$(pwd)/UTRan/bin/" >> ~/.bash_profile
source ~/.bash_profile
conda activate utran_env
```

It may be needed to apply "execution permission" to all bin executables: ```chmod +x UTRan/bin/*```

## Usage

## Inputs

## Outputs

## Contact
:bug::sos::speech_balloon:

To report bugs, to ask for help and to give any feedback, please contact **Pedro G. Nachtigall**: pedronachtigall@gmail.com

## Cite

If you use UTRan in your analysis, please cite the following (which is the manuscript that UTRan was designed to be used):

Nachtigall et al. (2022) Differences in PLA2 Constitution Distinguish the Venom of Two Endemic Brazilian Mountain Lanceheads, Bothrops cotiara and Bothrops fonsecai. Toxins, 14(4), 237. DOI:[https://doi.org/10.3390/toxins14040237](https://doi.org/10.3390/toxins14040237)
