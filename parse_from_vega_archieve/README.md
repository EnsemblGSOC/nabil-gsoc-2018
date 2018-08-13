# Parse from VEGA Archieve

These python scripts downloads files from vega arhieve, extracts the information and stores them in a local database

## Instructions

Please create two directories named 'data' and 'vega_files'

Running the following script downloads the .gtf and .fa files from ensembl ftp

```
$ python3 download_from_ftp.py
```

Next, running the following script extracts information from the .gtf files and stores them in a local sqlite database

```
$ python3 parse_gtf.py
```

## Requirements 

The following python3 modules are used, which can be downloaded using a pip3 command

```
wget
biopython
tqdm
```
