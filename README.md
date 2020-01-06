# differannotate
Differentiate genome annotations

Tools like [gffcompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml) and [ParsEval](https://aegean.readthedocs.io/en/v0.16.0/parseval.html) already exist for comparing genome annotations, but they only consider transcription features.
Differannotate was designed to perform a more generic analysis and evaluate as many features as possible.

For each feature category in the input GFF files, performance metrics and figures are generated at the region level and optionally the base level if a genome reference is included.
Metrics are also broken into unstranded (+/-), forward (+), and reverse (-) categories.
Regions are considered to be "true" if they reciprocally overlap with at least [90%] of a feature from the control annotation.

## Testing

```
git clone https://github.com/zyndagj/differannotate.git
cd differannotate
python setup.py test
```

> This produces output text during CLI tests, but should yield no errors.

## Installation

```
pip install git+https://github.com/zyndagj/differannotate.git
```

> Add the `--user` parameter for local (non-system) installs

## Usage

```
usage: differannotate [-h] -C GFF3 [-R FASTA] [--cname STR] -T GFF3 [GFF3 ...]
                      -N STR [STR ...] [-p INT] [--plot] [-e EXT] [-v]
                      [--temd]

A tool for comparing GFF3 annotations

optional arguments:
  -h, --help            show this help message and exit
  -C GFF3, --control GFF3
                        Control GFF3. All comparisons are relative to this
                        annotation.
  -R FASTA, --reference FASTA
                        Control reference to enable accurate base metrics
  --cname STR           Name of control GFF3
  -T GFF3 [GFF3 ...], --treat GFF3 [GFF3 ...]
                        Space separated list of GFF3 files for comparison
                        against the control
  -N STR [STR ...], --names STR [STR ...]
                        Space separated list of names for treatment GFF3 files
                        (name order must match file order)
  -p INT, --percent INT
                        Reciprocal percent overlap threshold [90]
  --plot                Plot venn diagrams of results
  -e EXT, --ext EXT     Figure extension [png]
  -v, --verbose         Enable verbose logging
  --temd                Analyze TE metadata
```

### Output

Tabular performance metrics will be printed to the CLI, and the following figures will be generated when possible.

- Base pair metrics for each category
  - `base*.[ext]` - Venn diagram of nucleotide logical relations
- Region metrics for each category
  - `region*.[ext]` - Venn diagram of region logical relations
  - `length_bp*.[ext]` - Boxplots showing length distribution across samples
  - `proportion*.[ext]` - Boxplots showing nucleotide proportion distribution across samples

### Example

```
differannoate -C tests/test_1.gff3 -R tests/test.fa -T tests/test_2.gff3 -N "Name with space" --plot
```

> Spaces can be used in sample names by encapsulating them in quotation marks.

## Citing
Zynda, G. J. (2020). Differannotate. GitHub repository. GitHub. Retrieved from https://github.com/zyndagj/differannotate

```
@misc{zynda2019differannotate,
    title     = {Differannotate},
    author    = {Zynda, Gregory J},
    url       = {https://github.com/zyndagj/differannotate},
    publisher = {GitHub},
    journal   = {GitHub repository},
    version   = {0.0.1},
    year      = {2020}
}
```
