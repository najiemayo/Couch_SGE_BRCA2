# Getting Started

There are 3 processes that need to be accounted for when running this analysis:

1. Alignment
2. Parsing SAM files
3. Annotating VCF file with CAVA
4. HTML Reports

## Alignment

The first step in the process is aligning to the reference genome. To make this simpler, we have developed a WDL-based
pipeline to execute.

First, we need to make sure our environmental variables are set, so we will source `PKG_PROFILE`.

```bash
source PKG_PROFILE
```

Now that we have our environment variables set, we need to define our inputs. Edit the `RAD51C_HTM/wdl/input.json` file.
Specifically, you need to:

1. Add the path to the BAM file
2. Define the output directory
3. Set the `cutadapt` options for trimming

To determine what `cutadapt` parameters to use, see [here](WhatToClip.md).

### Run WDL pipeline to generate SAM files

```shell
cromwell run # TODO: fill in actual executable
```

## Parsing SAM files

Once you have the SAM files, you can then parse out the relevant information and counts
using [CountReads.py](../python/CountReads.py). While there are instructions for how to run inside the `__main__`
portion of that script, I'll l briefly summarize the necessary steps here.

1. Define your input dictionaries
2. Define your acceptance criteria
3. Loop through each line of each SAM file for reads of interest
4. Count if they meet the acceptance criteria

There are two ways this process can be run:

* In other python scripts through its class
* Directly on the command-line

To use it in your own python scripts:

```python
# Import the class
from python.main import SamParser

# Define your argument variables
args = {
    'samfile': '../data/example.sam',
    'chrom': 'chr17',
    'seq': 'CTTCTGTTCAGCACTAGATGATATTCTTGGGGGTGGAGTGCCCTTAATGAAAACAACAGAAATTTGTGGTGCACCAGGTGTTGGAAAAACACAATTATGGTAAAATAAA',
    'positions': [56772452],
    'bases': ['A'],
    'target_bounds': [56772452, 56772560]
}
# For the definitions of these variables, see the cli explanation below

# Instantiate the class
SP = SamParser(args)
# Get the results
res = SP.process_sam()
```

The object returned will contain a dictionary containing:
* 'sample_name':  sample_name,
* 'total_count':  total_count, Total number of reads in the library
* 'off_target': Reads that didn't start in the right posiiton
* 'cigar_count':  Fails CIGAR filter (150 contiguous bases)
* 'mm_count':     Has too many amino acid changes
* 'invalid_count': Required DNA changes not present
* 'frameshift_count':  How many frameshifts
* 'nochange_count': How many with no changes
* 'synonymous_count': How many synonymous?
* 'missense_count': How many missense
* 'stop_count': How many stop-gained mutations
* 'mixed_count': How many Synonmymous-Missense hybrid reads
* 'multi_syn': How many reads with multiple synonymous variants were found
* 'usable_count':  synonymous_count + missense_count + stop_count

If you wanted to run it on the command line, here's how you would do it.

``bash
usage: main.py [-h] -s SAMFILE -c CHROM -p POSITIONS [POSITIONS ...] -b
               {A,C,G,T} [{A,C,G,T} ...] -r SEQ -P TARGET_BOUNDS TARGET_BOUNDS
               -z {F,R} -d CODON_DISTANCE
               [-V {DEBUG,INFO,WARNING,ERROR,CRITICAL}]

optional arguments:
  -h, --help            show this help message and exit
  -s SAMFILE, --sam SAMFILE
                        SAM file containing reads
  -c CHROM, --chrom CHROM
                        Chromosome of target
  -p POSITIONS [POSITIONS ...], --positions POSITIONS [POSITIONS ...]
                        locations of the required changes
  -b {A,C,G,T} [{A,C,G,T} ...], --bases {A,C,G,T} [{A,C,G,T} ...]
                        Non reference bases that must be located at positions
  -r SEQ, --ref SEQ     Target sequence (NOT containing mutations!
  -P TARGET_BOUNDS TARGET_BOUNDS, --target_bounds TARGET_BOUNDS TARGET_BOUNDS
                        Target sequence position [low, high] corresponding to
                        the `seq` variable
  -z {F,R}, --strand {F,R}
                        Forward (F) or Reverse (R) strand
  -d CODON_DISTANCE, --distance_to_first_full_codon CODON_DISTANCE
                        How many bases from the start of the fragment to the
                        first full codon triplet (relative to the + strand)
  -V {DEBUG,INFO,WARNING,ERROR,CRITICAL}, --verbose {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Set the logging level

```

```bash
SAM="data/3N_D12.sam"
SEQ='TCATCTTTCTGTTGACAGTATGCAGTTGGCAGTAGATGTGCAGATACCAGAATGTTTTGGAGGAGTcGCAGGTGAAGCAGTTTTTATTGATACAGAaGG'
BOUNDS="56774036 56774134"
POSITIONS="56774102 56774132"
BASES="C A"
CHROM="chr17"

python python/main.py -s $SAM -c $CHROM -p $POSITIONS -b $BASES -P $BOUNDS -r $SEQ
```

## Annotating VCF file with CAVA

One of the additional files you will need for making reports is an annotated TSV file that has CAVA annotations. These
are useful when grouping variants as synonymous or nonsense. There is a [helper script](../python/AllPossible.py) to
generate a VCF file when you only have a sequence and a coordinate.

```bash
cromwell run wdl/workflows/SeqToAnnoVCF.wdl -i wdl/input2.json
```




