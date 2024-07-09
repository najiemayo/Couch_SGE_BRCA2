# Authors: Steven N. Hart, Rohan David Gnanaolivu, Yohannes A. Tecleab
# Date: 7/9/2024
# Site Saturation Mutagenesis protocols

This repository contains code for high throughput SSM experiments. There are 4 primary directories

1. `htm_scripts`
    * All python scripts
        * AllPossible.py
            * Make a VCF for a particular sequence with every possible alternate allele
        * CountReads.py
            * Read in a SAM file and convert counts into a pkl object for use in EventPlots.md
2. `R`
    * All R code
        * EventPlots.md
            * Convert pickeled summary info into pretty graphs
3. `test`
    * Unit tests for all python code
4. `wdl`
    * WDL code for executing pipelines
5. `docs`
    * Instructions on how to use various tools in this package

The rest of this document should be reconfigured to describe the high level workflows and direct users to more specific
descriptions in the `docs` folder.

> NOTE: If you just want to write additional code for the project, make sure you load:
> `use_virtualenv("r-reticulate")` in R or
> `conda activate ~/SGE `
> for any python code.

# Setup environment

```shell
cd $PATH # Path to script directory
source ~/bin/activate
```

To run the analysis for another sample, see [GettingStarted.md](docs/GettingStarted.md).

### To run the pipeline

1. Convert Paired-End reads to SpliceAI

```shell
bash ./bash/run_ssm.sh -i $DATA_DIR -o $PWD/output -c $PWD/config/Gene.cfg
```

2. Run R Markdown scripts

```shell
bash ./bash/run_analysis_cfg.sh -i $DATA_DIR -o $PWD/output -c $PWD/config/Gene.cfg
```
