"""
Copyright 2019 Andrew Rambaut (a.rambaut@ed.ac.uk)
https://github.com/rambaut/Binlorry

This module contains the main script for Binlorry. It is executed when a user runs `binlorry`
(after installation) or `binlorry-runner.py` (directly from the source directory).

This file is part of Binlorry. Binlorry is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Binlorry is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Binlorry. If
not, see <http://www.gnu.org/licenses/>.
"""

import argparse
import os
import sys
import subprocess
import shutil
import re
from collections import defaultdict
from .misc import load_fasta_or_fastq, print_table, red, bold_underline, MyHelpFormatter, int_to_str
from .nanopore_read import NanoporeRead
from .version import __version__


def main():
    args = get_arguments()
    reads, check_reads, read_type = load_reads(args.input, args.verbosity, args.print_dest,
                                               args.check_reads)

    if args.native_barcodes:
        # construct a smaller set of search adapters with only the 12 barcodes to speed up the initial step
        # search_adapters = [a for a in ADAPTERS if '(full sequence)' not in a.name and '(forward)' not in a.name]
        search_adapters = NATIVE_BARCODES


    output_reads(reads, args.format, args.output, read_type, args.verbosity,
                 args.discard_middle, args.min_split_read_size, args.print_dest,
                 args.barcode_dir, args.barcode_labels, args.input, args.untrimmed, args.threads,
                 args.discard_unassigned)


def get_arguments():
    """
    Parse the command line arguments.
    """
    parser = argparse.ArgumentParser(description='Binlorry: a tool for binning sequencing reads into '
                                                 'folders based on header information or read properties.,
                                     formatter_class=MyHelpFormatter, add_help=False)
    main_group = parser.add_argument_group('Main options')
    main_group.add_argument('-i', '--input', required=True,
                            help='FASTA/FASTQ of input reads or a directory which will be '
                                 'recursively searched for FASTQ files (required)')
    main_group.add_argument('-o', '--output',
                            help='Filename for FASTA or FASTQ of trimmed reads (if not set, '
                                 'trimmed reads will be printed to stdout)')
    main_group.add_argument('--format', choices=['auto', 'fasta', 'fastq', 'fasta.gz', 'fastq.gz'],
                            default='auto',
                            help='Output format for the reads - if auto, the '
                                 'format will be chosen based on the output filename or the input '
                                 'read format')
    main_group.add_argument('-v', '--verbosity', type=int, default=1,
                            help='Level of progress information: 0 = none, 1 = some, 2 = lots, '
                                 '3 = full - output will go to stdout if reads are saved to '
                                 'a file and stderr if reads are printed to stdout')

    binning_group = parser.add_argument_group('Binning settings',
                                              'Control the binning of reads based on header fields')
    barcode_group.add_argument('-b', '--barcode_dir',
                               help='Reads will be binned based on their barcode and saved to '
                                    'separate files in this directory (incompatible with '
                                    '--output)')
    barcode_group.add_argument('--discard_unassigned', action='store_true',
                               help='Discard unassigned reads (instead of creating a "none" bin)')

    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                           help='Show this help message and exit')
    help_args.add_argument('--version', action='version', version=__version__,
                           help="Show program's version number and exit")

    args = parser.parse_args()

    if args.barcode_dir is not None and args.output is not None:
        sys.exit(
            'Error: only one of the following options may be used: --output, --barcode_dir')

    if args.untrimmed and args.barcode_dir is None:
        sys.exit('Error: --untrimmed can only be used with --barcode_dir')

    if args.barcode_dir is not None:
        args.discard_middle = True

    if args.output is None and args.barcode_dir is None:
        args.print_dest = sys.stderr
    else:
        args.print_dest = sys.stdout

    if args.threads < 1:
        sys.exit('Error: at least one thread required')

    return args


def load_reads(input_file_or_directory, verbosity, print_dest, check_read_count):

    # If the input is a file, just load reads from that file. The check reads will just be the
    # first reads from that file.
    if os.path.isfile(input_file_or_directory):
        if verbosity > 0:
            print('\n' + bold_underline('Loading reads'),
                  flush=True, file=print_dest)
            print(input_file_or_directory, flush=True, file=print_dest)
        reads, read_type = load_fasta_or_fastq(input_file_or_directory)
        if read_type == 'FASTA':
            reads = [NanoporeRead(x[2], x[1], '') for x in reads]
        else:  # FASTQ
            reads = [NanoporeRead(x[4], x[1], x[3]) for x in reads]
        check_reads = reads[:check_read_count]

    # If the input is a directory, assume it's an Albacore directory and search it recursively for
    # fastq files. The check reads will be spread over all of the input files.
    elif os.path.isdir(input_file_or_directory):
        if verbosity > 0:
            print('\n' + bold_underline('Searching for FASTQ files'),
                  flush=True, file=print_dest)
        fastqs = sorted([os.path.join(dir_path, f)
                         for dir_path, _, filenames in os.walk(input_file_or_directory)
                         for f in filenames
                         if f.lower().endswith('.fastq') or f.lower().endswith('.fastq.gz')])
        if not fastqs:
            sys.exit('Error: could not find fastq files in ' +
                     input_file_or_directory)
        reads = []
        read_type = 'FASTQ'
        check_reads = []
        check_reads_per_file = int(round(check_read_count / len(fastqs)))
        for fastq_file in fastqs:
            if verbosity > 0:
                print(fastq_file, flush=True, file=print_dest)
            file_reads, _ = load_fasta_or_fastq(fastq_file)
            file_reads = [NanoporeRead(x[4], x[1], x[3]) for x in file_reads]

            albacore_barcode = get_albacore_barcode_from_path(fastq_file)
            for read in file_reads:
                read.albacore_barcode_call = albacore_barcode
            reads += file_reads
            check_reads += file_reads[:check_reads_per_file]
        if verbosity > 0:
            print('', flush=True, file=print_dest)

    else:
        sys.exit('Error: could not find ' + input_file_or_directory)

    if verbosity > 0:
        print(int_to_str(len(reads)) + ' reads loaded\n\n',
              flush=True, file=print_dest)
    return reads, check_reads, read_type


def output_reads(reads, out_format, output, read_type, verbosity, discard_middle,
                 min_split_size, print_dest, barcode_dir, barcode_labels, input_filename,
                 untrimmed, threads, discard_unassigned):
    if verbosity > 0:
        trimmed_or_untrimmed = 'untrimmed' if untrimmed else 'trimmed'
        if barcode_dir is not None:
            verb = 'Saving '
            destination = 'barcode-specific files'
        elif output is None:
            verb = 'Outputting '
            destination = 'stdout'
        else:
            verb = 'Saving '
            destination = 'file'
        print(bold_underline(verb + trimmed_or_untrimmed + ' reads to ' + destination),
              flush=True, file=print_dest)

    if out_format == 'auto':
        if output is None:
            out_format = read_type.lower()
            if barcode_dir is not None and input_filename.lower().endswith('.gz'):
                out_format += '.gz'
        elif '.fasta.gz' in output.lower():
            out_format = 'fasta.gz'
        elif '.fastq.gz' in output.lower():
            out_format = 'fastq.gz'
        elif '.fasta' in output.lower():
            out_format = 'fasta'
        elif '.fastq' in output.lower():
            out_format = 'fastq'
        else:
            out_format = read_type.lower()

    gzipped_out = False
    gzip_command = 'gzip'
    if out_format.endswith('.gz') and (barcode_dir is not None or output is not None):
        gzipped_out = True
        out_format = out_format[:-3]
        if shutil.which('pigz'):
            if verbosity > 0:
                print('pigz found - using it to compress instead of gzip')
            gzip_command = 'pigz -p ' + str(threads)
        else:
            if verbosity > 0:
                print('pigz not found - using gzip to compress')

    # Output reads to barcode bins.
    if barcode_dir is not None:
        if not os.path.isdir(barcode_dir):
            os.makedirs(barcode_dir)
        barcode_files = {}
        barcode_read_counts, barcode_base_counts = defaultdict(
            int), defaultdict(int)
        for read in reads:
            barcode_name = read.barcode_call
            if discard_unassigned and barcode_name == 'none':
                continue
            if out_format == 'fasta':
                read_str = read.get_fasta(
                    min_split_size, discard_middle, untrimmed, barcode_labels)
            else:
                read_str = read.get_fastq(
                    min_split_size, discard_middle, untrimmed, barcode_labels)
            if not read_str:
                continue
            if barcode_name not in barcode_files:
                barcode_files[barcode_name] = \
                    open(os.path.join(barcode_dir,
                                      barcode_name + '.' + out_format), 'wt')
            barcode_files[barcode_name].write(read_str)
            barcode_read_counts[barcode_name] += 1
            if untrimmed:
                seq_length = len(read.seq)
            else:
                seq_length = read.seq_length_with_start_end_adapters_trimmed()
            barcode_base_counts[barcode_name] += seq_length
        table = [['Barcode', 'Reads', 'Bases', 'File']]

        for barcode_name in sorted(barcode_files.keys()):
            barcode_files[barcode_name].close()
            bin_filename = os.path.join(
                barcode_dir, barcode_name + '.' + out_format)

            if gzipped_out:
                if not os.path.isfile(bin_filename):
                    continue
                bin_filename_gz = bin_filename + '.gz'
                if os.path.isfile(bin_filename_gz):
                    os.remove(bin_filename_gz)
                try:
                    subprocess.check_output(gzip_command + ' ' + bin_filename,
                                            stderr=subprocess.STDOUT, shell=True)
                except subprocess.CalledProcessError:
                    pass
                bin_filename = bin_filename_gz

            table_row = [barcode_name, int_to_str(barcode_read_counts[barcode_name]),
                         int_to_str(barcode_base_counts[barcode_name]), bin_filename]
            table.append(table_row)

        if verbosity > 0:
            print('')
            print_table(table, print_dest, alignments='LRRL',
                        max_col_width=60, col_separation=2)

    # Output to all reads to stdout.
    elif output is None:
        for read in reads:
            read_str = read.get_fasta(min_split_size, discard_middle, barcode_labels=barcode_labels) if out_format == 'fasta' \
                else read.get_fastq(min_split_size, discard_middle, barcode_labels=barcode_labels)
            print(read_str, end='')
        if verbosity > 0:
            print('Done', flush=True, file=print_dest)

    # Output to all reads to file.
    else:
        if gzipped_out:
            out_filename = 'TEMP_' + str(os.getpid()) + '.fastq'
        else:
            out_filename = output
        with open(out_filename, 'wt') as out:
            for read in reads:
                read_str = read.get_fasta(min_split_size, discard_middle, barcode_labels=barcode_labels) if out_format == 'fasta' \
                    else read.get_fastq(min_split_size, discard_middle, barcode_labels=barcode_labels)
                out.write(read_str)
        if gzipped_out:
            subprocess.check_output(gzip_command + ' -c ' + out_filename + ' > ' + output,
                                    stderr=subprocess.STDOUT, shell=True)
            os.remove(out_filename)
        if verbosity > 0:
            print('\nSaved result to ' + os.path.abspath(output), file=print_dest)

    if verbosity > 0:
        print('', flush=True, file=print_dest)


def output_progress_line(completed, total, print_dest, end_newline=False, step=10):
    if step > 1 and completed % step != 0 and completed != total:
        return
    progress_str = int_to_str(completed) + ' / ' + int_to_str(total)
    if total > 0:
        percent = 100.0 * completed / total
    else:
        percent = 0.0
    progress_str += ' (' + '%.1f' % percent + '%)'

    end_char = '\n' if end_newline else ''
    print('\r' + progress_str, end=end_char, flush=True, file=print_dest)
