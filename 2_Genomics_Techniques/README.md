# Session 2: Genomics Techniques
This session will provide an overview for understanding FASTQ file formats, interpreting read quality, mapping to a reference, and calling genomic variants. It is derived from the [Data Wrangling and Processing for Genomics Lesson](https://datacarpentry.org/wrangling-genomics/) from [Data Carpentry](https://datacarpentry.org/lessons/) which provides a much better and in-depth view of these topics. However, rather than running the analyses on a pre-imaged Amazon Web Service instance, we will try to facilitate analysis on your local machine, though it will require some more patience with software installation. A smaller dataset (instead of the <i>E. coli</i> used for the Data Carpentry lession) is provided to help each step run more quickly. 

## Prior to the session: 
1. **(WINDOWS ONLY): Download and install Cygwin.**  Cygwin provides a Linux-like environment for Windows. It is built on top of Windows, but allows you to interact with your computer using Linux commands. Follow the instructions [here](http://www1.udel.edu/CIS/105/pcline/07J/useful-links/cygwin/) to install. Once the installation is complete, you can start using Cygwin by launching it using the desktop shortcut or from the start menu.

2. **Download and install FastQC.** Click 'Download Now' on the [FastQC website](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and follow the instructions under 'Installation and setup instructions'.

3. **Download and install BWA, samtools, and bcftools.** Unfortunately, these programs are only available for MacOS and Linux systems. Participants using Windows machines are welcome to follow along during the workshop even if they are not able to run the analyses on their own machine. I have not tested it with these programs, but some users might be interested in trying out [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10).  You can choose which Linux distribution to install; Ubuntu 18.04 LTS should work well.

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
See Section 2 of the Data Carpentry lesson [here](https://datacarpentry.org/wrangling-genomics/02-quality-control/index.html) for an in-depth description.

1. Download the two example files in this repository, `SH009_R1.fastq` and `SH009_R2.fastq` to a new working directory. These files include a subset of genomic reads from a *Striga hermonthica* individual collected in a field of maize in Mumias, Kenya by Emily Bellis, Sylvia Mutinda, Calvins Odero, and Steven Runo in 2018.

2. Navigate to the directory on your computer where these files were downloaded. Inspect the first sequence in each read file using the `head` command. Is this sequence of good quality?
```
$ head -n 4 SH009_R1.fastq
```

3. How many sequences are in each file? We can use the command `grep` with the `-c` flag to find out.
```
$ grep -c '@A00755' SH009_R1.fastq
$ grep -c '@A00755' SH009_R2.fastq
```

## 2b: Assessing read quality with FastQC:
1. Ensure that FastQC is installed properly. We can run it with the `-h` flag to print the help screen. 
```
$ fastqc -h
```

2. Now, generate a FastQC report for each of the read files. When run on the command-line, FastQC will generate an `.html` file as output to display the results.  The `*` is a wildcard character, indicating to the program to perform the operation on any file that ends in `.fastq`. 
```
$ fastq *.fastq
```

3. Check out the FastQC output by clicking on the `.html` file to open it.

Alternatively, [FastQC can be run interactively](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt) using the graphical user interface. 

## 2c. Mapping reads to the reference with BWA (MacOS/Linux):
See Section 4 [here](https://datacarpentry.org/wrangling-genomics/04-variant_calling/index.html) for full details. Note, we will be skipping the quality trimming and filtering steps (Section 3) today; the reads we are using have already been trimmed and filtered using [BBDuk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/).

1. First, we must download the reference. It will be easiest if you download to the same directory where you have the read files. The reference we are using as an example is actually just 1 sequence, that of a transcript for a pectin methylesterase inhibitor gene from *Striga hermonthica*.  Typically though we would use a genome reference consisting of multiple chromosomes or contigs.

2. Inside the directory where you have the file, inspect the reference sequence. It is a `.fasta` file, with the sequence id preceded by a `>` on one line, followed by the sequence on the following lines. In this case the sequence is wrapped across several lines.  
```
$ cat Sther_PMEI.fasta 
```

3. Next, index the reference sequence with `bwa`. Unless you have added the bwa program to your PATH environment variable (outside the scope of this workshop), you may need to specify the full file path to where you downloaded the bwa program. 
```
$ bwa index Sther_PMEI.fasta
```

4. Then, map the forward and reverse read files for sample 'SH009' to the reference with the `bwa mem` algorithm. After the mapping has finished, we can inspect the resulting `.sam` file since it is a human readable format.
```
$ bwa mem SH009_R1.fastq SH009_R2.fastq > SH009.aligned.sam 
$ cat SH009.aligned.sam
```

5. This file is actually pretty small! Most `.sam` files can be much bigger. See how big it is with:
```
$ du -sh SH009.aligned.sam
```

6. We will create a compressed version (a `.bam` file) that is smaller and more efficient for processing by computers. We use `samtools view` to convert the file, and `samtools sort` to sort it. Again, you may have to specfify the full path to the `samtools` program.
```
$ samtools view -S -b SH009.aligned.sam > SH009.aligned.bam
$ samtools sort SH009.aligned.bam > SH009.aligned.sorted.bam
```

7. See how much smaller the `.bam` file is?
```
$ du -sh SH009.aligned.sam
```

8. Once we have the `.bam` file we can use `samtools flagstat` to check how many reads mapped to the reference.
```
$ samtools flagstat SH009.aligned.sorted.bam
```

9. We can also visualize the alignment with `samtools tview`. You can scroll with the arrow keys on your keyboard, or type `?` for a menu of options. 
```
$ samtools tview SH009.aligned.sorted.bam Sther_PMEI.fasta
```

## 2d. Variant calling:
There are many programs available to perform variant calling; for now we will keep following the Data Carpentry tutorial and use `bcftools`. 

1. `bcftools` first counts read coverage at each position. The flag `-O b` tells `bcftools` to generate a bcf format output file, `-o` specifies where to write the output file, and `-f` flags the path to the reference genome.
```
$ bcftools mpileup -O b -o SH009_raw.bcf -f Sther_PMEI.fasta SH009.aligned.sorted.bam
```

2. Next, call single nucleotide polymorphisms (SNPs).  Specify ploidy with `--ploidy`, `-m` allows for multiallelic and rare-variant calling, `-v` tells the program to output variant sites only (not every site in the genome), and `-o` specifies where to write the output file:
```
$ bcftools call --ploidy 2 -m -v -o SH009_variants.vcf SH009_raw.bcf 
```

3. Most pipelines also involve a filtering step after variant calling, which can vary according to the specific project. For now, let's inspect the `.vcf` file, one of the standard file format for genomic variants!
```
$ cat SH009_variants.vcf
```
