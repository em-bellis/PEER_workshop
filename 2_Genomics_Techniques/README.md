# Session 2: Genomics Techniques
This session will provide an overview for understanding FASTQ file formats, interpreting read quality, mapping to a reference, and calling genomic variants. It is derived from the [Data Wrangling and Processing for Genomics Lesson](https://datacarpentry.org/wrangling-genomics/) from [Data Carpentry](https://datacarpentry.org/lessons/) which provides a much better and in-depth view of these topics. However, rather than running the analyses on a pre-imaged Amazon Web Service instance, we will try to facilitate analysis on your local machine, though it will require some patience with software installation. A smaller dataset (instead of the <i>E. coli</i> used for the Data Carpentry lession) is provided to help each step run more quickly. 

## Prior to the session: 
1. **(WINDOWS ONLY): Download and install Cygwin.**  Cygwin provides a Linux-like environment for Windows. It is built on top of Windows, but allows you to interact with your computer using Linux commands. Follow the instructions [here](http://www1.udel.edu/CIS/105/pcline/07J/useful-links/cygwin/) to install. Once the installation is complete, you can start using Cygwin by launching it using the desktop shortcut or from the start menu.

2. **Download and install FastQC.** Click the 'Download Now' on the [FastQC website](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and follow the instructions under 'Installation and setup instructions'. Note, FastQC will require 

3. **Download and install BWA, samtools, and bcftools.** Unfortunately, these programs are only available for MacOS and Linux systems. Participants are welcome to follow along during the workshop even if they are not able to run the analyses on their own machine. Alternatively, some users might be interested in trying out [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10).  You can choose which Linux distribution to install; Ubuntu 18.04 LTS is recommended. However, this workshop hasn't been extensively tested with WSL.

- **[BWA](http://bio-bwa.sourceforge.net):**
```
$ curl -OL http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2
$ tar jxvf bwa-0.7.17.tar.bz2
$ cd bwa-0.7.17
$ make
```

- **[samtools](http://samtools.sourceforge.net):**
```
$ curl -OkL https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
$ tar jxvf samtools-1.9.tar.bz2
$ cd samtools-1.9
$ make
```
- **[bcftools](http://samtools.github.io/bcftools/bcftools.html):**
```
$ curl -OkL https://github.com/samtools/bcftools/releases/download/1.8/bcftools-1.8.tar.bz2
$ tar jxvf bcftools-1.8.tar.bz2
$ cd bcftools-1.8
$ make
```
---

## 2a: The FASTQ format:
See [here](https://datacarpentry.org/wrangling-genomics/02-quality-control/index.html) for an in-depth description.

Download the two example files in this repository, x and x. 

Look at the top of each read file.

How many sequences are in each file?

What do you notice about the length of the sequences?

## 2b: Assessing read quality with FastQC:
Ensure that FastQC is installed properly.

Generate a FastQC report on the command-line.

Alternatively, FastQC can be run interactively using the graphical user interface.

## 2c. Mapping reads to the reference with BWA:
See (here)[https://datacarpentry.org/wrangling-genomics/04-variant_calling/index.html] for full details. Note, we will be skipping the quality trimming and filtering steps today; the reads in this repository have already been trimmed and filtered using [BBDuk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/). 

## 2d. Variant calling:
We will continue with the Data Carpentry tutorial using **bcftools** for variant calling. Note, there are many programs available to perform variant calling; **bcftools** is just one of them.

## 2e. Assess the alignment:

