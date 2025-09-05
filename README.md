# GFFcompare

Parse and compare the content of 2 GFF files. This allows you to quickly determine the difference in content between 2 annotated genomes. This is particularly useful to compare to very close bacterial strains.

## Requirement

BcBio package is required for GFF parsing.

```
pip install bcbio-gff
```

## Usage
```
gffcomp.py -1 genome1.gff -2 genome2.gff
```