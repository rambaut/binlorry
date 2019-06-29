"""
Copyright 2019 Andrew Rambaut (a.rambaut@ed.ac.uk)
https://github.com/rambaut/Binlorry

This module contains the main script for BinLorry. It is executed when a user runs `binlorry`
(after installation) or `binlorry-runner.py` (directly from the source directory).

This file is part of BinLorry. BinLorry is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. BinLorry is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with BinLorry. If
not, see <http://www.gnu.org/licenses/>.
"""

import argparse
import gzip
import os
import sys

from .misc import bold_underline, MyHelpFormatter, int_to_str, \
    get_sequence_file_type, get_compression_type
from .read import Read
from .version import __version__


def main():
    args = get_arguments()

    process_files(args.input, args.format, args.output, args.verbosity, args.print_dest)

def get_arguments():
    """
    Parse the command line arguments.
    """
    parser = argparse.ArgumentParser(description='BinLorry: a tool for binning sequencing reads into '
                                                 'folders based on header information or read properties.',
                                     formatter_class=MyHelpFormatter, add_help=False)

    main_group = parser.add_argument_group('Main options')
    main_group.add_argument('-i', '--input', required=True,
                            help='FASTA/FASTQ of input reads or a directory which will be '
                                 'recursively searched for FASTQ files (required)')
    main_group.add_argument('-o', '--output',
                            help='Filename for FASTA or FASTQ of filtered reads')
    main_group.add_argument('--format', choices=['auto', 'fasta', 'fastq', 'fasta.gz', 'fastq.gz'],
                            default='auto',
                            help='Output format for the reads - if auto, the '
                                 'format will be chosen based on the output filename or the input '
                                 'read format')
    main_group.add_argument('-v', '--verbosity', type=int, default=1,
                            help='Level of progress information: 0 = none, 1 = some, 2 = lots, '
                                 '3 = full - output will go to stdout if reads are saved to '
                                 'a file and stderr if reads are printed to stdout')

    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                           help='Show this help message and exit')
    help_args.add_argument('--version', action='version', version=__version__,
                           help="Show program's version number and exit")

    args = parser.parse_args()

    if args.output is None:
        # output is to stdout so print messages to stderr
        args.print_dest = sys.stderr
    else:
        args.print_dest = sys.stdout

    return args


def process_files(input_file_or_directory, out_format, out_filename, verbosity, print_dest):

    read_files = []

    if os.path.isfile(input_file_or_directory):
        read_files.append(input_file_or_directory)

        if verbosity > 0:
            print('\n' + bold_underline('Processing reads'), flush=True, file=print_dest)

    # If the input is a directory, search it recursively for fastq files.
    elif os.path.isdir(input_file_or_directory):
        if verbosity > 0:
            print('\n' + bold_underline('Searching for FASTQ files'), flush=True, file=print_dest)

        read_files = sorted([os.path.join(dir_path, f)
                             for dir_path, _, filenames in os.walk(input_file_or_directory)
                             for f in filenames
                             if f.lower().endswith('.fastq') or f.lower().endswith('.fastq.gz') or
                             f.lower().endswith('.fasta') or f.lower().endswith('.fasta.gz')])
        if not read_files:
            sys.exit('Error: could not find FASTQ/FASTA files in ' + input_file_or_directory)

    else:
        sys.exit('Error: could not find ' + input_file_or_directory)

    with open(out_filename, 'wt') as out_file:
        for read_file in read_files:

            file_type = get_sequence_file_type(read_file)

            if get_compression_type(read_file) == 'gz':
                open_func = gzip.open
            else:  # plain text
                open_func = open

            if verbosity > 0:
                print(read_file, flush=True, file=print_dest)

            if file_type == 'FASTA':
                # For FASTA files we need to deal with line wrapped sequences...
                with open_func(read_file, 'wt') as in_file:

                    name = ''
                    sequence = ''

                    for line in in_file:
                        line = line.strip()

                        if not line:
                            continue

                        if line[0] == '>':  # Header line = start of new read
                            if name:
                                write_read(out_file, out_format, Read(name.split()[0], sequence, name))
                                sequence = ''
                            name = line[1:]
                        else:
                            sequence += line

                    if name:
                        write_read(out_file, out_format, Read(name.split()[0], sequence, name))

            else: # FASTQ
                with open_func(read_file, 'rt') as in_file:
                    for line in in_file:
                        full_name = line.strip()[1:]
                        sequence = next(in_file).strip()
                        next(in_file) # spacer line
                        qualities = next(in_file).strip()
                        write_read(out_file, out_format, Read(full_name, sequence, qualities))

        if verbosity > 0:
            print('', flush=True, file=print_dest)

        if verbosity > 0:
            print('\nSaved result to ' + os.path.abspath(out_filename), file=print_dest)


def write_read(out_file, out_format, read):
    read_str = read.get_fasta() if out_format == 'fasta' else read.get_fastq()
    out_file.write(read_str)


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
