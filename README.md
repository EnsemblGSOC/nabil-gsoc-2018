# nabil-gsoc-2018
## Transcript Comparisons
In Comparative Genomics we compare an unknown gene with some other known genes, for better inference of biological properties of that unknown gene. Identification of Gene Orthology Relation is the most important task of Comparative Genomics, as they tend to preserve similar molecular and biological functions. Thus if we can establish orthology relationships between two genes, we can obtain valuable evolutionary history of the two genes. However, with advanced sequencing depth and expansion in transcriptome data, genes are no longer the proper units for interrogation in functional conservation, evolutionary events, and expressional patterns, especially in the field of alternative splicing. As the accumulation of transcriptomic data, alternative splicing is taken into account in the assignments of gene orthologs and the orthology is suggested to be further considered at transcript level. Whether gene or transcript orthology, exons are the basic units that represent the whole gene structure; however, there is not much reported study on how to build exon level orthology in a whole genome scale. Therefore, it is essential to establish a transcription oriented gene orthology algorithm.


## Brief Description
Though Orthologous genes are very closely related to each other, some slight to moderate deviations are observed among the several transcripts of the same gene. This phenomenon is more visible for the Eukaryotes, where a plethora of protein isoforms are formed from the same gene due to alternative splicing. In this work we attempted to model the different variants of alternative splicing and based on it made a ranking of the orthologous transcripts of a transcript. The rankings are formulated with the objective that the most similar transcript will be ranked higher.

We have developed a complete software tool for comparing orthologous transcripts. Our goal is to enable biologists or bioinformaticians to study the various transcripts and identify the most related ones. In addition to providing apis, we have developed a very user friendly gui to benefit the researchers. 


## Algorithm

A manuscript is being prepared explaining the algorithm and the results.....


## Requirements

### Python

The entire system is built upon Python, more precisely Python3. The visualization tool is a local Flask app.  

The following modules were used.

* Flask
* BioPython
* tqdm

These are listed in the **requirements.txt** file.

The python modules are automatically installed when the **setup.sh** file is run.

### Perl

Perl is used for downloading the data through Ensembl api. 

BioPerl, Perl DBI, DBD modules and most importantly Ensembl api is required.

### C++ 
Some computationally expensive dynamic programming modules are implemented in C++, for faster computation. No external librarires were used.

The C++ codes are converted into executables when the **setup.sh** file is run.

### HTML, CSS
HTML and CSS are used in the front-end of the application. The styles are adopted from [MaterializeCSS](https://materializecss.com). The CSS files are bundled with the package, thus no more extra hustle.

### Javascript
Javascript controls the dynamicity of the application. [MaterializeCSS](https://materializecss.com) and [jQuery](https://jquery.com/) have been used, which are bundled with the package.

### MUSCLE
[MUSCLE](https://www.ebi.ac.uk/Tools/msa/muscle/) is used for sequence alignment. A platfrom dependent exe file is required.


## Installation

First of all make the shell files executable

```
$ chmod +x ./setup.sh
$ chmod +x ./run.sh
```

Almost all the requirements are dealt with when the **setup.sh** is run, which needs sudo privilege.

```
$ sudo ./setup.sh
```

To use MUSCLE, you would first need to download the appropriate exe file for your platform from [here](https://www.drive5.com/muscle/downloads.htm). 

Next you would need to specify the path to the exe file in the following text file : 
    exon_alignments/muscle_path.txt


## Data 

The data is obtained from Ensembl. As an user you can do one of the following:

1. Download the data yourself by using the perl scripts in the data_acquisition directory
2. Download the data from [here](https://drive.google.com/open?id=1X51bFX-WVbsp346byt36VW3kCUb6mPtb) and save it in the path data_acquisition/data/db.db 


## Running the app

Simply run the run.sh file
```
$ ./run.sh
```


## Tutorial 
A comprehensive step by step tutorial can be found [here](https://ensemblgsoc.github.io/nabil-gsoc-2018/)
