# Session 2: Genomics Techniques
This session will provide an overview for understanding FASTQ file formats, interpreting read quality, mapping to a reference, and calling genomic variants. It is derived from the [Data Wrangling and Processing for Genomics Lesson](https://datacarpentry.org/wrangling-genomics/) from [Data Carpentry](https://datacarpentry.org/lessons/) which provides a much better and in-depth view of these topics. However, rather than running the analyses on a pre-imaged Amazon Web Service instance, we will try to facilitate analysis on your local machine, though it will require some patience with software installation. A smaller dataset (instead of the <i>E. coli</i> used for the Data Carpentry lession) is also provided. 

## Prior to the session: 
1. **(WINDOWS ONLY): Download and install Cygwin.**  Cygwin provides a Unix-like environment for Windows. It is built on top of Windows, but allows you to interact with your computer using Unix commands. Follow the instructions (here)[http://www1.udel.edu/CIS/105/pcline/07J/useful-links/cygwin/] to install. Once the installation is complete, you can start using Cygwin by launching it using the desktop shortcut or from the start menu.

Alternatively, some users might wish to consider installing (Windows Subsystem for Linux)[https://docs.microsoft.com/en-us/windows/wsl/install-win10].  You can choose which Linux distribution to install; Ubuntu 18.04 LTS is recommended.

2.  **Download and install FASTQC.**   

