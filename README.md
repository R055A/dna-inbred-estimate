# DNA Inbred Estimate

Estimate the F_ROH representation of inbreeding/shared-ancestry from genetic DNA data kit file.

* ~25.0% - Parent-child/Sibling inbreeding (example: King Charles the II of Hapsburg)
* ~12.5% - First-cousin inbreeding
* ~6.25% - Second-cousin inbreeding
* < ~6.0% - Distant parental relatedness, endogamy
* < ~3.0% - Ancestral endogamy, founder effect

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

## Run

```bash
python dna_inbred_estimate.py -file <file-name> [OPTIONS]
```
