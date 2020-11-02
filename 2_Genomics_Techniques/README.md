# Session 2: Genomics Techniques
This session will provide an overview for understanding FASTQ file formats, interpreting read quality, mapping to a reference, and calling genomic variants. It is derived from the [Data Wrangling and Processing for Genomics Lesson](https://datacarpentry.org/wrangling-genomics/) from [Data Carpentry](https://datacarpentry.org/lessons/) which provides a much better and in-depth view of these topics. However, rather than running the analyses on a pre-imaged Amazon Web Service instance, we will try to facilitate analysis on your local machine, though it will require some patience with software installation. A smaller dataset (instead of the <i>E. coli</i> used for the Data Carpentry lession) is also provided. 

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

## Overview: 
