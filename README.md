# BinLorry

BinLorry is a flexible tool for binning and filtering sequencing reads into distinct files. Reads can be binned and filtered by any attributes encoded in their headers, documented in a CSV file or by length.

## Installing

Simply install with pip:

```
pip3 install binlorry
```

Run:
```
binlorry --help
```

### Install from repository

Clone the repository:
```
git clone https://github.com/rambaut/binlorry.git
```

Install:
```
pip3 install ./binlorry
```

### Run without installation

BinLorry can also be run directly from the repository clone, without installation:
```
git clone https://github.com/rambaut/binlorry.git
python binlorry/binlorry-runner.py -h
```
However, ensure that the ``pandas`` package is installed before use.


## Quick Usage Examples

```
binlorry -i reads/ -o barcode --bin-by barcode --filter-by barcode BC01 BC02 -n 550 -x 750
```

This would read all FASTQ or FASTA files in the directory `reads`, bin by the header field `barcode`, but only if this is `BC01` or `BC02` and if the length is between 550 and 750 nucleotides.
It would use the file name prefix `barcode` resulting in the files: `barcode_BC01.fastq` and `barcode_BC02.fastq`

```
binlorry -i my_file.fastq -t my_file.csv --out-report -o filtered --filter-by reference Type_1 -n 550 -x 750
```

The above example will take in reads from ``my_file.fastq`` and a csv report ``my_file.csv``. Assuming that ``my_file.csv`` has at least the structure shown below, and that the read names in the csv match those in the input read file, BinLorry will filter reads and output only those with Type_1 reference between 550 and 750 bases in length.

| read_name                             | reference | 
|:--------------------------------------|-----------:| 
| f66db89e-de96-4fa7-813a-6c5a89586100 | Type_1    | 
| a39069c5-c493-45f8-9fa8-49eccb5c1807 | Type_1    | 
| 868efa99-f4c1-4a68-87a9-196a44b997e0 | Type_2    | 


```
binlorry -i path/to/my_fastq_dir -t path/to/my_csv_dir \
--out-report -o path/to/binned/barcode \
--filter-by barcode BC01 --bin-by barcode -n 1000 -x 2000
```

Assuming you have reports in the csv dir corresponding to the read files in the fastq dir, binlorry will recursively search both directories, matching the csv and fastq files based on filename stem. This command will then filter reads only containing BC01 and output a csv report corresponding to the reads presented in the output fastq file.

## Command line interface
```
usage: binlorry -i INPUT [-t CSV_FILE] -o OUTPUT [-v VERBOSITY]
                         [--bin-by FIELD [FIELD ...]]
                         [--filter-by FILTER [FILTER ...]] [-n MIN] [-x MAX]
                         [-h] [--version]

Main options:
  -i INPUT, --input INPUT
                          FASTA/FASTQ of input reads or a directory which will
                          be recursively searched for FASTQ files (required)
  -t INPUT_CSV, --index-table INPUT_CSV
                           A CSV file with metadata fields for reads (otherwise these are assumed
                           to be in the read headers). This can also include a file and line number to improve performance. Assumes read name is first column of the csv.'
  -o OUTPUT, --output OUTPUT
                          Output filename (or filename prefix)
  -r REPORT, --out-report REPORT
                          Output a subsetted csv report along with the fastq. (Default: False)
                          Only implemented for use in conjunction with -t option.
  -f FORCE_OUTFILES, --force-output FORCE_OUTFILES
                          Output binned/ filtered files even if empty. (default: False)
                          Usage: only a single binning factor with a corresponding filter factor.
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
                          the reads by. Multiple filter-by options can be
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
