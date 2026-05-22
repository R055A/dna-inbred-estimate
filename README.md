# DNA Inbred Estimate

Estimate the F_ROH representation of inbreeding/shared-ancestry from genetic DNA data kit file.

## Requirements

> A raw genetic DNA data kit file used for genealogical research.

Supported DNA kit file formats:
* Ancestry.com (array v2, converter v1)
* MyHeritage.com (format v1)
* __TODO__: [test more formats](https://github.com/R055A/dna-inbred-estimate/issues/1)

## Install

```bash
pip install -r requirements.txt
```

## Run

```bash
python dna_inbred_estimate.py -file <file-name>
```

## Options

```bash
python dna_inbred_estimate.py -help
```

## Contribute

Run the following code formatting before making a PR:

```bash
black .
```

```bash
ruff check .
```
