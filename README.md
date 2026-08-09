# DNA Inbred Estimate

Compute F_ROH from raw genetic DNA data kit file using [PLINK](https://github.com/chrchang/plink-ng).

F_ROH estimates provide an approximate representation of genomic inbreeding or shared parental ancestry:
 
* ~25.000% - Parent-child or siblings (example: King Charles II of Spain)
* ~12.500% - Half siblings, grandparent-grandchild, uncle–niece or aunt–nephew
* ~6.2500% - First cousins, half uncle–niece or half aunt–nephew
* ~3.1250% - First cousins once removed
* ~1.5625% - Second cousins

## Requirements

> A raw genetic DNA data kit file used for genealogical research.

Supported DNA kit file formats:
* Ancestry.com (array v2, converter v1)
* MyHeritage.com (format v1)
* 23andMe.com
* __TODO__: [test more formats](https://github.com/R055A/dna-inbred-estimate/issues/1)

## Install

```bash
pip install -r requirements.txt
```

```bash
python utils/install_plink.py
```

## Run

```bash
python app.py -file <file-name> [OPTIONS]
```
