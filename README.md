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
```
usage: UTRan.py [-h] [-t TRANSCRIPTS] [-c CDS] [-r READS] [-o [OUTPUT]]
                [-cz [CHECKZEROS]] [-m [MISMATCHES]] [-cp [CHANGINGPOINT]]
                [-ws [WINDOWSIZE]] [-minSF [MINSIZEFIVE]]
                [-minST [MINSIZETHREE]] [-ts [TOOSHORT]] [-mq [MAPQUALITY]]
                [-minR [MINREAD]] [-maxR [MAXREAD]] [-T THREADS]

UTRan - UTR sequence annotator using de novo transcriptome assembly.

optional arguments:
  -h, --help            show this help message and exit
  -t TRANSCRIPTS, --transcripts TRANSCRIPTS
                        Mandatory - Fasta file of assembled contigs. Use the
                        full length transcripts assembled.
  -c CDS, --cds CDS     Mandatory - file with CDS annotation. This file can be
                        in BED, GTF or GFF format. It also accepts a FASTA
                        format with the CDSs, but it is important to note
                        that, in this case, the CDSs must be extracted
                        directly from the assembled contigs to ensure
                        identical residues between CDSs and contigs.
  -r READS, --reads READS
                        Mandatory - Fastq(.gz) file of UNPAIRED reads. Please
                        use merged reads if working with paired-end data. See
                        UTRan repository for more details.
  -o [OUTPUT], --output [OUTPUT]
                        Optional - path to output folder (e.g.,
                        /path/to/output/folder/). If not declared, it will be
                        created at the working directory.
                        [Default="UTRan_output"]
  -cz [CHECKZEROS], --CheckZeros [CHECKZEROS]
                        Optional - Turn on this parameter to perform an
                        additional check to detect zero coverage in the middle
                        of the annotated CDS. Zero coverage in the middle of
                        CDS may represent a putative chimeric transcript. To
                        perform this step, specify "True". [Default="False"]
  -m [MISMATCHES], --mismatches [MISMATCHES]
                        Optional - Number of allowable mismatches to keep in
                        alignment. [Default = 0]
  -cp [CHANGINGPOINT], --changingpoint [CHANGINGPOINT]
                        Optional - Value between 0 and 1 to be used as
                        threshold to detect changing points in the read
                        coverage of contigs to detect the UTRs. [Default =
                        0.6]
  -ws [WINDOWSIZE], --WindowSize [WINDOWSIZE]
                        Optional - Window size to be used to detect changing
                        points. [Default = 10]
  -minSF [MINSIZEFIVE], --minSizeFive [MINSIZEFIVE]
                        Optional - Minimum size of 5'UTR. [Default = 30]
  -minST [MINSIZETHREE], --minSizeThree [MINSIZETHREE]
                        Optional - Minimum size of 3'UTR. [Default = 30]
  -ts [TOOSHORT], --tooShort [TOOSHORT]
                        Optional - Minimum length of clipped reads. [Default =
                        100]
  -mq [MAPQUALITY], --mapQuality [MAPQUALITY]
                        Optional - Minimum mapping quality. Please, note that
                        reads with multiple mappings get assigned a quality of
                        0 by bwa. [Default = 0]
  -minR [MINREAD], --minRead [MINREAD]
                        Optional - Minimum read length. [Default = 100]
  -maxR [MAXREAD], --maxRead [MAXREAD]
                        Optional - Maximum read length. [Default = 1000]
  -T THREADS, --threads THREADS
                        Optional - Number of threads. [Default = 2]
```

## Outputs
```
UTRan_output/
├── 3utr_sequences.fa
├── 5utr_sequences.fa
├── coverage.txt
├── report.txt
└── Plots
    ├── transcript_1.pdf
    ├── ...
    └── transcript_N.pdf
```

## Contact
:bug::sos::speech_balloon:

To report bugs, to ask for help and to give any feedback, please contact **Pedro G. Nachtigall**: pedronachtigall@gmail.com

## Cite

If you use UTRan in your analysis, please cite the following (which is the manuscript that UTRan was designed to be used):

Nachtigall et al. (2022) Differences in PLA2 Constitution Distinguish the Venom of Two Endemic Brazilian Mountain Lanceheads, Bothrops cotiara and Bothrops fonsecai. Toxins, 14(4), 237. DOI:[https://doi.org/10.3390/toxins14040237](https://doi.org/10.3390/toxins14040237)
