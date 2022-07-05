# Computational Biology - Salmonella Outbreak

## Table of contents
* [General info](#general-info)
* [Requirements](#technologies)
* [Setup](#setup)
* [Results](#results)


## General info
The tool takes as input two different FASTA files (the Wild strain and the Variantstrain) and give in output a list of SNPs, which are a germline substitutions of a single nucleotide at a specific position in the genome.

	
## Requirements
Project is written in Python and works with the following version:
*  python_version >= '3.6'
All the other requirements can be installed through the terminal writing:
```
$ pip install -r requirements.txt
```
Main libraries:
* levenshtein==0.16.0
* matplotlib==3.5.0
* scipy==1.7.3
* termcolor==1.1.0
* biopython==1.79
* numpy==1.21.4

## Setup
Before running the programme, reading the  [Report.pdf](Report.pdf) is highly recommended.

To run this project, after installing the requirements, open the terminal and launch it as:
```
$ python main.py file1.fasta file2.fasta -k 50
```
where:
* **file1.fasta** is the Wild Type file
* **file2.fasta** is the Variant Type file
* **k** is the lenght of the k-mer (default value 50)

Once launched, after the computation of the **file1** and **file2** the tool shows a plot of the aboundance of the k-mers. After each plot an integer input parameter **p** (sequencing error parameter) has to be inserted. The files used can be downloaded at: https://cloud-ljk.imag.fr/index.php/s/HkxDLozHRcqBcqz

## Results
Results are presented as:
```
SNP
    R: TTACATGCCAATACAATGTAGGCTGCTCTACACCTAGCTTCTGGGCGAGTTTACGGGTTGTTAAACCTTCGATTCCGACCTCATTAAGCAGCTCTAATGCG
    V: TTACATGCCAATACAATGTAGGCTGCTCTACACCTAGCTTCTGGGCGAGGGGACGGGTTGTTAAACCTTCGATTCCGACCTCATTAAGCAGCTCTAATGCG
   	 Distance: 3
------------------------------------------------------
SNP
    R: CGCATTAGAGCTGCTTAATGAGGTCGGAATCGAAGGTTTAACAACCCGTAAACTCGCCCAGAAGCTAGGTGTAGAGCAGCCTACATTGTATTGGCATGTAA
    V: CGCATTAGAGCTGCTTAATGAGGTCGGAATCGAAGGTTTAACAACCCGTCCCCTCGCCCAGAAGCTAGGTGTAGAGCAGCCTACATTGTATTGGCATGTAA
   	 Distance: 3
```



