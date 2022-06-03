#!/usr/bin/env python3
'''
UTRan - UTR sequence annotator
Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)
'''

##Import modules
import re
#from optparse import OptionParser
import argparse
import datetime as dt
import subprocess
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import pathos.multiprocessing as mp
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
except:
    print('''UTRan was not able to run due to the ERROR below:
    The biopython package is not properly installed or it is not active in your actual environment.
    Please, install biopython properly (check https://biopython.org/ for more details), or active the environment where biopython is installed!''')
    quit()
try:
    from Bio.Alphabet import generic_dna
except ImportError:
    generic_dna = None

##Functions
#generate the reverse complement of a sequence
def _ReverseComplement_(sequence):
    if generic_dna:
        DNA = Seq(sequence, generic_dna)
    else:
        DNA = Seq(sequence)
    RC = str(DNA.reverse_complement())
    return RC

#parse fasta file
def _ParseFasta_(fasta):
    final = {}
    for record in SeqIO.parse(fasta,"fasta"):
        ID = str(record.id)
        SEQ = str(record.seq)
        final[ID] = SEQ
    return final

#parse GTF/GFF/BED file
def _ParseGTF_(CDS):
    final = {}
    if CDS.endswith("bed"):
        a = open(CDS, "r")
        for line in a:
            line1 = line.rstrip().split("\t")
            final[line1[0]] = [line1[1], line1[2], line1[5]]
        a.close()
    else:
        a = open(CDS, "r")
        for line in a:
            if "\tCDS\t" in line:
                line1 = line.rstrip().split("\t")
                final[line1[0]] = [line1[3], line1[4], line1[6]]
        a.close()
    return final

#convert all contigs into "plus" in a temp file
def _Convert2Plus_(T, Annotation, outF):
    Aconverted = {}
    tmp = open(outF+"transcripts_temp.fa", "w")
    for k in Annotation.keys():
        if Annotation[k][2] == "+":
            tmp.write(">"+k+"\n"+T[k]+"\n")
            start = int(Annotation[k][0])
            end = int(Annotation[k][1])
            Aconverted[k] = [start, end, "+"]
        if Annotation[k][2] == "-":
            RC = _ReverseComplement_(T[k])
            tmp.write(">"+k+"\n"+RC+"\n")
            start = len(RC)-int(Annotation[k][1])
            end = len(RC)-int(Annotation[k][0])
            Aconverted[k] = [start, end, "+"]
    tmp.close()
    return Aconverted

#find PAS positions in the sequence mc[AATAAA, ATTAAA] & lc[TATAAA, AGTAAA, TGTAAA]
    #canonical PAS (AAUAAA) and its variants (AUUAAA,UAUAAA, AGUAAA, UUUAAA, CAUAAA, AAUACA, AAUGAA,AAUAUA, GAUAAA, and UGUAAA)
def _FindPAS_(contig):
    pas = []
    p = re.compile('[A][A|T][T][A][A][A]')
    for m in p.finditer(contig):
        pas.append((m.start(), m.group()))
    if pas == []:
        p = re.compile('[T|A][A|G][T][A][A][A]')
        for m in p.finditer(contig):
            pas.append((m.start(), m.group()))
    return pas

def _FinalChart_(x, transcript, Futr, fiveprime, st, end, Tutr, threeprime, pas, outF):

    x_data = list(x['pos'])
    y_data = list(x['cov'])

    cds = list(range(st, end))
    poscds = [0.8]*len(cds)

    PAS = []
    for i in pas:
        PAS.append(i[0])
    pospas = [0.75]*len(PAS)

    fig, axes = plt.subplots(2)
    fig.suptitle(transcript)

    axes[0].bar(x_data, y_data, color="blue", ec="blue")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("Coverage")
    axes[0].set_xlim(x_data[0], x_data[-1])
    axes[1].plot(cds, poscds, color='yellow', linewidth=6)
    axes[1].plot([Futr, st], [0.8, 0.8], 'k-', color='green', linewidth=5)
    axes[1].plot([end, Tutr], [0.8, 0.8], 'k-', color='magenta', linewidth=5)
    if threeprime != "":
        axes[1].text(end+10, 0.85, "3\'UTR", rotation=45)
    axes[1].text(st+150, 0.85, "CDS", rotation=45)
    if fiveprime != "":
        axes[1].text(Futr, 0.85, "5\'UTR", rotation=45)
    if PAS != []:
        axes[1].plot(PAS, pospas, 'o', color="red")
        axes[1].text(PAS[0]-300, 0.70, "Poly-A signal")
    axes[1].set_ylim(0.6, 1.0)
    axes[1].set_yticks([0.6, 1.0])
    axes[1].set_xlabel("Contig position")
    axes[1].set_ylabel("Annotation")
    axes[1].set_xlim(x_data[0], x_data[-1])

    plt.savefig(outF+"Plots/"+transcript+".pdf")
    plt.close()

#main function
def _UTRan_(transcripts, CDS, reads, outF, threads, CheckZeros, mismatches, ChangingPoint, WindowSize, minSizeFive, minSizeThree, tooShort, mapQuality, minRead, maxRead):
    ##step 1 - convert all contigs into the plus strand (based on annotation)
        #this annotation will return the putative 5' and 3' UTRs from contigs
    print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> Interpreting CDS annotation data")
    T = _ParseFasta_(transcripts)
    if CDS.split(".")[-1] in ["fa","fasta","fas"]:
        Annotation = {}
        C = _ParseFasta_(CDS)
        #if ids in CDS.fasta dn contigs.fasta are equal
        for k in C.keys():
            if k in T.keys():
                contig = T[k]
                if C[k].find(contig):
                    m = re.search(C[k], contig)
                    if m != None:
                        start = m.span()[0]+1
                        end = m.span()[1]
                        strand = "+"
                        Annotation[id] = [start, end, strand]
                if C[k].find(contig) == -1:
                    RC = _ReverseComplement_(contig)
                    if C[k] in RC:
                        m = re.search(C[k], RC)
                        if m != None:
                            start = m.span()[1]
                            end = m.span()[0]
                            strand = "-"
                            Annotation[id] = [start, end, strand]
            #else, find the seqs from CDS.fasta in the contigs.fasta
            if k not in T.keys():
                for id in T.keys():
                    contig = T[id]
                    if C[k].find(contig):
                        m = re.search(C[k], contig)
                        if m != None:
                            start = m.span()[0]
                            end = m.span()[1]
                            strand = "+"
                            Annotation[id] = [start, end, strand]
                    if C[k].find(contig) == -1:
                        RC = _ReverseComplement_(contig)
                        if C[k] in RC:
                            m = re.search(C[k], RC)
                            if m != None:
                                start = len(RC)-int(m.span()[1])
                                end = len(RC)-int(m.span()[0])
                                strand = "-"
                                Annotation[id] = [start, end, strand]
    if CDS.split(".")[-1] in ["bed","gff","gtf"]:
        Annotation = _ParseGTF_(CDS)

    #create a temporary file with all contigs converted
    Aconverted = _Convert2Plus_(T, Annotation, outF)
    Tconverted = _ParseFasta_(outF+"transcripts_temp.fa")

    ##step 2 - perform the alignment and post-processing
        #if there is changing points in the middle of CDS print a warning message and exclude this transcript from the analysis
    print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> Mapping and read filtering step")

    grepNM = "grep -E 'NM:i:[0-" + str(mismatches) + "][[:space:]]|^@'"
    # Build the bwa index
    command = "bwa index "+outF+"transcripts_temp.fa"
    subprocess.call(command,shell=True)
    # Generate the initial sam alignment
    command = "bwa mem -M -t "+str(threads) +" -R \'@RG\\tID:Sample\\tSM:Reads\' "+outF+"transcripts_temp.fa "+ reads +" | "+ grepNM  +" > "+outF+"tmp1.sam"
    subprocess.call(command,shell=True)
    # Create a sorted bam file
    command = "picard SortSam INPUT="+outF+"tmp1.sam OUTPUT="+outF+"tmp2.bam SORT_ORDER=coordinate USE_JDK_DEFLATER=true USE_JDK_INFLATER=true"
    subprocess.call(command,shell=True)
    command =  "picard BuildBamIndex INPUT="+outF+"tmp2.bam USE_JDK_DEFLATER=true USE_JDK_INFLATER=true"
    subprocess.call(command,shell=True)
    os.remove(outF+"tmp1.sam")
    # Remove overclipped reads
    command = "picard CreateSequenceDictionary REFERENCE="+outF+"transcripts_temp.fa OUTPUT="+outF+"transcripts_temp.dict USE_JDK_DEFLATER=true USE_JDK_INFLATER=true"
    subprocess.call(command,shell=True)
    command = "samtools faidx "+outF+"transcripts_temp.fa"
    subprocess.call(command,shell=True)
    command = "gatk PrintReads -R "+outF+"transcripts_temp.fa -I "+outF+"tmp2.bam -RF OverclippedReadFilter --filter-too-short "+str(tooShort)+" --dont-require-soft-clips-both-ends -RF MappingQualityReadFilter --minimum-mapping-quality "+str(mapQuality)+" -RF ReadLengthReadFilter --min-read-length "+str(minRead)+" --max-read-length "+str(maxRead)+" -O "+outF+"transcripts_temp.bam"
    subprocess.call(command,shell=True)
    # Calculate coverage
    command = "bedtools genomecov -d -ibam "+outF+"transcripts_temp.bam > "+outF+"coverage.txt"
    subprocess.call(command,shell=True)
    subprocess.call("rm tmp*[bs]a[im]",shell=True)
    os.remove(outF+"tmp2.bam")

    #read coverage data
    X = pd.read_csv(outF+"coverage.txt",sep='\t',names=['transcript','pos','cov'])
    print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> Identified coverage data for " + str(len(set(X['transcript']))) + " transcripts")
    print("\t-> Filtering transcripts with zero coverage in the middle of CDS (putative chimeras)")
    TR = list(set(X['transcript']))
    #check for any sites with no coverage in the CDS (putative chimeras) - if "-cz True"
    if CheckZeros != False:
        zeroBad = []
        for i in TR:
            CDSpos = Aconverted[i]
            X2 = X[X['transcript'] == i]
            X3 = X2[X2['pos'] >= CDSpos[0]]
            X4 = X3[X3['pos'] <= CDSpos[1]]
            zeros = list(X4[X4['transcript'] == i]['cov']).count(0)
            if zeros :
                zeroBad.append(i)
        #remove the contigs with zero coverage sites in the CDS
        TR = [x for x in TR if x not in zeroBad]
        print("\t-> "+str(len(TR))+" transcripts passed through the check step of CDSs with zero coverage")

    print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> Identifying changing points in read coverage of UTR regions")

    ##step 3 - annotate the pas signals in putative 3'UTRs
    PAS = {}
    for k in TR:
        end = Aconverted[k][1]
        contig = Tconverted[k]
        pas = _FindPAS_(contig)
        PAS.setdefault(k, [])
        for i in pas:
            if i[0] > end:
                PAS[k].append(i)
        #PAS[k] = pas

    ##step 4 - identify the changing points in read coverage of UTR regions
    #checking the changing points in read coverage in the UTR regions
    FIVE = {}
    THREE = {}
    def _WorkingCode_(i):
        FIVE = {}
        THREE = {}
        report = {}
        report.setdefault(i, [])
        Fout = open(outF+"5utr_sequences.fa","a+")
        Tout = open(outF+"3utr_sequences.fa","a+")
        CDSpos = Aconverted[i]
        X2 = X[X['transcript'] == i]

        contig = Tconverted[i]
        st = int(Aconverted[i][0])
        end = int(Aconverted[i][1])

        print("\t-> Analyzing transcript "+i+" [CDS:"+str(st)+"-"+str(end)+"]")
        report[i].append("\t-> CDS:"+str(st)+"-"+str(end))

        #work with 5'UTR
        coverages = []
        fiveend = st
        count = 1
        avg = 10.0
        #if the 5utr region is longer than 'minSizeFive' nts
        if st > minSizeFive:
            report[i].append("\t-> checked the changing point in read coverage to find the correct 5\'UTR")
            for n in reversed(range(st+1)[1:]):
                if count < minSizeFive+1:
                    coverages.append(X2.loc[X2['pos'] == n, 'cov'].iloc[0])
                    count += 1
                    if X2.loc[X2['pos'] == n, 'cov'].iloc[0] == 0:
                        fiveend = n
                        report[i].append("\t\t5\'UTR pos:"+str(fiveend)+"-"+str(st))
                        break
                if count >= minSizeFive+1:
                    count += 1
                    avg = np.mean(coverages[-WindowSize:])
                    if X2.loc[X2['pos'] == n, 'cov'].iloc[0] == 0:
                        fiveend = n
                        report[i].append("\t\t5\'UTR pos:"+str(fiveend)+"-"+str(st))
                        break
                    if X2.loc[X2['pos'] == n, 'cov'].iloc[0] <= avg*ChangingPoint:
                        fiveend = n #changing point -1?
                        report[i].append("\t\t5\'UTR pos:"+str(fiveend)+"-"+str(st))
                        break
                    if n == 1:
                        fiveend = 1
                        report[i].append("\t\t5\'UTR pos:"+str(fiveend)+"-"+str(st))
                        break
                    coverages.append(X2.loc[X2['pos'] == n, 'cov'].iloc[0])

        #if the putative 5utr region is shorter than 'minSizeFive' nts, it will get all nts (truncated)
        if st <= minSizeFive:
            fiveend = 1
            report[i].append("\t-> truncated 5\'UTR")
            report[i].append("\t\t5\'UTR pos:"+str(fiveend)+"-"+str(st))

        FIVE[i] = [fiveend, X2.loc[X2['pos'] == fiveend, 'cov'].iloc[0], avg]

        #work with 3'UTR
        coverages = []
        threeend = len(contig)
        count = 1
        avg = 10.0
        #if the 3utr region is longer than 'minSizeThree' nts
        if len(contig)-end > minSizeThree:
            report[i].append("\t-> checked the changing point in read coverage to find the correct 3\'UTR")
            for n in range(end, len(contig)):
                if count < minSizeThree+1:
                    coverages.append(X2.loc[X2['pos'] == n, 'cov'].iloc[0])
                    count += 1
                    if X2.loc[X2['pos'] == n, 'cov'].iloc[0] == 0:
                        threeend = n
                        report[i].append("\t\t3\'UTR pos:"+str(end)+"-"+str(threeend))
                        break
                if count >= minSizeThree+1:
                    count += 1
                    avg = np.mean(coverages[-WindowSize:])
                    if X2.loc[X2['pos'] == n, 'cov'].iloc[0] == 0:
                        threeend = n
                        report[i].append("\t\t3\'UTR pos:"+str(end)+"-"+str(threeend))
                        break
                    if X2.loc[X2['pos'] == n, 'cov'].iloc[0] <= avg*ChangingPoint:
                        threeend = n #changing point +1?
                        report[i].append("\t\t3\'UTR pos:"+str(end)+"-"+str(threeend))
                        break
                    if n == len(contig):
                        threeend = len(contig)
                        report[i].append("\t\t3\'UTR pos:"+str(end)+"-"+str(threeend))
                        break
                    coverages.append(X2.loc[X2['pos'] == n, 'cov'].iloc[0])
        #if the putative 3utr region is shorter than 'minSizeThree' nts, it will get all nts (truncated)
        if len(contig)-end <= minSizeThree:
            threeend = len(contig)
            report[i].append("\t-> truncated 3\'UTR")
            report[i].append("\t\t3\'UTR pos:"+str(end)+"-"+str(threeend))

        THREE[i] = [threeend, X2.loc[X2['pos'] == threeend, 'cov'].iloc[0], avg]
        if threeend == len(contig):
            report[i].append("\t\t3\'UTR pos:"+str(end)+"-"+str(threeend))

        for pas in PAS[i]:
            if pas[0] < threeend and pas[0] > threeend-40:
                report[i].append("\t\tchanging point near to the PAS [pos:"+str(pas[0])+","+pas[1]+"]")
        if PAS[i] != []:
            report[i].append("\t\tall PAS identified in the full contig: "+";".join([str(x) for x in PAS[i]]))

        fiveprime = contig[FIVE[i][0]:st]
        if fiveprime != "":
            Fout.write(">"+i+"\n"+fiveprime+"\n")
        threeprime = contig[end:THREE[i][0]] #+1? contig[end:THREE[i][0]+1]
        if threeprime != "":
            Tout.write(">"+i+"\n"+threeprime+"\n")

        _FinalChart_(X2, i, FIVE[i][0], fiveprime, st, end, THREE[i][0], threeprime, PAS[i], outF)
        Fout.close()
        Tout.close()

        Rout = open(outF+"report.txt","a+")
        for k in report.keys():
            Rout.write(">"+k+"\n"+"\n".join(report[k])+"\n")
        Rout.close()

    os.mkdir(outF+"Plots/")

    def _RunAnalysis_(TR):
        p = mp.ProcessPool(nodes=threads)
        p.map(_WorkingCode_, TR)
        p.clear()
        p.close()
        p.join()

    _RunAnalysis_(TR)

    ##step 5 - correlate the changing points with the pas signals to identify the correct 3'UTR and use the changing point to identify the correct 5'UTR
        #if the transcrpt ends, we will have a truncated UTR, which is not a problem
        #if the read coverage reach to 0, it is where the UTR end in the transcript (it will no be necessary to correlated with pas in the case of 3'UTR)

##files and arguments to test script
#mandatory
#transcripts = "/home/nachtigall/Desktop/posdoc_butantan/miRNA_seq/UTR/SB0265VGL_contigs_CLEAN.fasta"
#CDS = "/home/nachtigall/Desktop/posdoc_butantan/miRNA_seq/UTR/SB0265VGL_CDS_final.fasta"
#reads = "/media/nachtigall/413C3C5E0D0853F9/Linux/bothrops_datasets/merged_reads/SB0265VGL_merged_reads.assembled.fastq.gz"
#outF = "/home/nachtigall/Desktop/posdoc_butantan/miRNA_seq/UTR/SB0265_utran/"

#optional - here I am setting the default for the optional parameters
#threads = 4
#mismatches = 0
#ChangingPoint = 0.6
#minSizeFive = 30
#minSizeThree = 100
#tooShort = 150
#mapQuality = 0
#minRead = 150
#maxRead = 1000

#options.output_folder
#if os.path.isdir(outF) == False:
#    os.mkdir(outF)

#print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> Starting UTRan analysis!")
#_UTRan_(transcripts, CDS, reads, outF, threads, mismatches, ChangingPoint, minSizeFive, minSizeThree, tooShort, mapQuality, minRead, maxRead)
#print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> Finished!")

##Options
def _Main_():
    parser = argparse.ArgumentParser(description='UTRan - UTR sequence annotator using de novo transcriptome assembly.')
    #Mandatory parameters
    parser.add_argument("-t","--transcripts",
                        type=str,
                        default = None,
                        help="Mandatory - Fasta file of assembled contigs. Use the full length transcripts assembled.")
    parser.add_argument("-c","--cds",
                        type=str,
                        default = None,
                        help="Mandatory - file with CDS annotation. This file can be in BED, GTF or GFF format. It also accepts a FASTA format with the CDSs, but it is important to note that, in this case, the CDSs must be extracted directly from the assembled contigs to ensure identical residues between CDSs and contigs.")
    parser.add_argument("-r","--reads",
                        type=str,
                        default = None,
                        help="Mandatory - Fastq(.gz) file of UNPAIRED reads. Please use merged reads if working with paired-end data. See UTRan repository for more details.")
    #Optional parameters
    parser.add_argument("-o","--output",
                        nargs='?',
                        type=str,
                        default=None,
                        help="Optional - path to output folder (e.g., /path/to/output/folder/). If not declared, it will be created at the working directory. [Default=\"UTRan_output\"]")
    parser.add_argument("-cz","--CheckZeros",
                        nargs='?',
                        type=str,
                        default=False,
                        help="Optional - Turn on this parameter to perform an additional check to detect zero coverage in the middle of the annotated CDS. Zero coverage in the middle of CDS may represent a putative chimeric transcript. To perform this step, specify \"True\". [Default=\"False\"]")
    parser.add_argument("-m","--mismatches",
                        nargs='?',
                        type=int,
                        default=0,
                        help="Optional - Number of allowable mismatches to keep in alignment. [Default = 0]")
    parser.add_argument("-cp","--changingpoint",
                        nargs='?',
                        type=float,
                        default=0.6,
                        help="Optional - Value between 0 and 1 to be used as threshold to detect changing points in the read coverage of contigs to detect the UTRs. [Default = 0.6]")
    parser.add_argument("-ws","--WindowSize",
                        nargs='?',
                        type=int,
                        default=10,
                        help="Optional - Window size to be used to detect changing points. [Default = 10]")
    parser.add_argument("-minSF","--minSizeFive",
                        nargs='?',
                        type=int,
                        default=30,
                        help="Optional - Minimum size of 5\'UTR. [Default = 30]")
    parser.add_argument("-minST","--minSizeThree",
                        nargs='?',
                        type=int,
                        default=30,
                        help="Optional - Minimum size of 3\'UTR. [Default = 30]")
    parser.add_argument("-ts","--tooShort",
                        nargs='?',
                        type=int,
                        default=100,
                        help="Optional - Minimum length of clipped reads. [Default = 100]")
    parser.add_argument("-mq","--mapQuality",
                        nargs='?',
                        type=int,
                        default=0,
                        help="Optional - Minimum mapping quality. Please, note that reads with multiple mappings get assigned a quality of 0 by bwa. [Default = 0]")
    parser.add_argument("-minR","--minRead",
                        nargs='?',
                        type=int,
                        default=100,
                        help="Optional - Minimum read length. [Default = 100]")
    parser.add_argument("-maxR","--maxRead",
                        nargs='?',
                        type=int,
                        default=1000,
                        help="Optional - Maximum read length. [Default = 1000]")
    parser.add_argument("-T","--threads",
                        type=int,
                        default=2,
                        help="Optional - Number of threads. [Default = 2]")

    args = parser.parse_args()

    transcripts = args.transcripts
    CDS = args.cds
    reads = args.reads
    outF = args.output
    CheckZeros = args.CheckZeros
    mismatches = args.mismatches
    ChangingPoint = args.changingpoint
    WindowSize = args.WindowSize
    minSizeFive = args.minSizeFive
    minSizeThree = args.minSizeThree
    tooShort=args.tooShort
    mapQuality=args.mapQuality
    minRead=args.minRead
    maxRead=args.maxRead
    threads = args.threads

    if transcripts == None or CDS == None or reads == None:
        print("""
   __  ____________
  / / / /_  __/ __ \\____ _____
 / / / / / / / /_/ / __ `/ __ \\
/ /_/ / / / / _, _/ /_/ / / / /
\____/ /_/ /_/ |_|\__,_/_/ /_/

>>>> UTRan v1.0 April 2021 <<<<
****Use -h for help!****

BASIC USAGE (with mandatory arguments):
UTRan.py -t transcripts.fasta -c cds(.fasta or .bed/gtf/gff) -r reads.fastq(.gz)
""")
        quit()

    if outF != None:
        if not outF.endswith("/"):
            outF += "/"
    if outF == None:
        CWD = os.getcwd()
        outF = CWD+"/UTRan_output/"

    if os.path.isdir(outF) == False:
        os.mkdir(outF)

    if transcripts != None and os.path.isfile(transcripts) == False:
            print('''UTRan was not able to run due to the ERROR below:
        The transcripts file indicated is not a valid file. Please, indicate a valid fasta file to the \"-i\" option.''')
            quit()

    if CDS != None and os.path.isfile(CDS) == False:
            print('''UTRan was not able to run due to the ERROR below:
        The CDS file indicated is not a valid file. Please, indicate a valid fasta/annotation file to the \"-c\" option.''')
            quit()

    if reads != None and os.path.isfile(reads) == False:
            print('''UTRan was not able to run due to the ERROR below:
        The reads file indicated is not a valid file. Please, indicate a valid fastq file to the \"-r\" option.''')
            quit()

    if transcripts != None and CDS != None and reads != None:
        print("""
   __  ____________
  / / / /_  __/ __ \\____ _____
 / / / / / / / /_/ / __ `/ __ \\
/ /_/ / / / / _, _/ /_/ / / / /
\____/ /_/ /_/ |_|\__,_/_/ /_/

""")
        print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> Starting UTRan (v1.0 April 2021) analysis!")

        print("\ttranscript file ->",transcripts)
        print("\tCDS file ->",CDS)
        print("\treads ->",reads)
        print("\toutput folder ->",outF)
        print("\tcheck for zero coverage in CDS ->",CheckZeros)
        print("\tmismatches ->",mismatches)
        print("\tchanging point threshold ->",ChangingPoint)
        print("\twindows size of coverage to be analyzed ->",WindowSize)
        print("\tminimum size of 5\'UTR ->",minSizeFive)
        print("\tminimum size of 3\'UTR ->",minSizeThree)
        print("\tminimum length of clipped reads ->",tooShort)
        print("\tminimum mapping quality ->",mapQuality)
        print("\tminimum read length ->",minRead)
        print("\tmaximum read length ->",maxRead)
        print("\tnumber of threads ->",threads)

        _UTRan_(transcripts, CDS, reads, outF,
                threads, CheckZeros,
                mismatches, ChangingPoint, WindowSize,
                minSizeFive, minSizeThree,
                tooShort, mapQuality, minRead, maxRead)

        print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> Finished!")

_Main_()
#END
