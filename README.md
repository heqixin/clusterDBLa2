# clusterDBLalpha
A python script to generate an otu table from cleaned DBLa sequences.

## Installation
The script depends on python2.7. 
All the necessary programs are in the `third_party` folder.
The hmm file for dbla search is under `data` folder

## Usage
The fasta file is assumed to be either in the format used by Thomas Rask's pipeline 
or the updated pipeline which uses the Usearch format to assign reads to isolates.  
you can either only perform cluster or only perform upsTypesSearch by choosing the 
corresponding options.  
output from DBLa domain search is saved as the file called `_DBLaUpsType.csv`, which has four columns:

|Type|Domain|Grouping|Score|
|----|------|:-------:|:---:|
|M00123:163:000000000-B5J7V:1:2117:23936:9166;sample=S5MRS1424.MID54-54.P5.Jun17;size=1063;|DBLa0.22|BC|208.8|
|M00566:68:000000000-ARD58:1:1119:9923:17310;sample=S3MRS2115.MID89-89.P4.May16;size=148;|DBLa1.7|A|175.9|
|M00566:161:000000000-B5HBJ:1:1118:3983:14249;sample=RS4MRS2224.MID5-5.P4.Jun17;size=448;|DBLa1.1|A|203.6|

```
usage: clusterDBLa_2.py [-h] -o OUTPUTDIR -r READ [--perID PERID] [--cpu CPU]
                      [--minEval MINEVAL] [--cluster_only]
                      [--upsTypesSearch_only] [--verbose]

Cluster cleaned DBLalpha sequence tags, translate tags and assign to DBLa ups
types.

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUTDIR, --outputDir OUTPUTDIR
                        location of output directory. Will be created if it
                        doesn't exist
  -r READ, --read READ  location of fasta file containing sequences.
  --perID PERID         percentage ID threshold. (default=0.96)
  --cpu CPU             number of cpus to use. (default=1)
  --minEval MINEVAL     minimum E-value for Hmmer match. (default=1E-30)
  --cluster_only        whether to perform clustering only. (default=False)
  --upsTypesSearch_only
                        whether to perform translation and search upsTypes
                        only. If using upsTypesSearch only, -r option should
                        be your centroid fasta of type sequences
                        (default=False)
  --verbose             print verbose output (default=False)
```
