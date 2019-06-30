Binlorry is a tool for binning sequencing reads into folders or files by
attributes encoded in their headers or by size.

```
usage: binlorry -i INPUT [-t CSV_FILE] -o OUTPUT [-v VERBOSITY]
                         [--bin-by FIELD [FIELD ...]]
                         [--filter-by FILTER [FILTER ...]] [-n MIN] [-x MAX]
                         [-h] [--version]

Main options:
  -i INPUT, --input INPUT
                          FASTA/FASTQ of input reads or a directory which will
                          be recursively searched for FASTQ files (required)
  -o OUTPUT, --output OUTPUT
                          Output filename (or filename prefix)
  -v VERBOSITY, --verbosity VERBOSITY
                          Level of progress information: 0 = none, 1 = some, 2
                          = lots, 3 = full - output will go to stdout if reads
                          are saved to a file and stderr if reads are printed
                          to stdout (default: 1)

Binning/Filtering options:
  --bin-by FIELD [FIELD ...]
                          Specify header field(s) to bin the reads by. For
                          multiple fields these will be nested in order
                          specified.
  --filter-by FILTER [FILTER ...]
                          Specify header field and accepted values to filter
                          the reads by. Multiple--filter-by options can be
                          specified.
  -n MIN, --min-length MIN
                          Filter the reads by their length, specifying the
                          minimum length.
  -x MAX, --max-length MAX
                          Filter the reads by their length, specifying the
                          maximum length.

Help:
  -h, --help              Show this help message and exit
  --version               Show program's version number and exit
```

Examples:

`binlorry -i reads/ -o barcode --bin-by barcode --filter-by barcode BC01 BC02 -n 550 -x 750`

This would read all FASTQ or FASTA files in the directory `reads`, bin by the header field `barcode`, but only if this is `BC01` or `BC02` and if the length is between 550 and 750 nucleotides.
It would use the file name prefix `barcode` resulting in the files: `barcode_BC01.fastq` and `barcode_BC02.fastq`
