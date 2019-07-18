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

import pandas as pd

from .misc import bold_underline, MyHelpFormatter, int_to_str, \
    get_sequence_file_type, get_compression_type
from .version import __version__


def main():
    '''
    Entry point for BinLorry. Gets arguments, processes them and then calls process_files function
    to do the actual work.
    :return:
    '''
    args = get_arguments()

    filters = []
    if args.filter_by:
        for filter_def in args.filter_by:
            values = filter_def[1:]
            filters.append({ 'field': filter_def[0], 'values': values})

    if args.verbosity > 0:
        for filter in filters:
            print("Filter reads unless " + filter['field'] + (" is one of: " if len(filter['values']) > 1 else " is "), end = '')
            for index, value in enumerate(filter['values']):
                print((", " if index > 0 else "") + value, end = '')
            print()

        if args.min_length and args.max_length:
            print("Filter reads unless length between " + str(args.min_length) + " and " + str(args.max_length))
        elif args.min_length:
            print("Filter reads unless length at least " + str(args.min_length))
        elif args.max_length:
            print("Filter reads unless length not more than " + str(args.max_length))

        if args.bin_by:
            print("Bin reads by ", end = '')
            for index, bin in enumerate(args.bin_by):
                print((", " if index > 0 else "") + bin, end = '')
            print()

    index_table = None

    if args.index_table_file:
        index_table = read_index_table(args.index_table_file)

    process_files(args.input, index_table, args.bin_by, filters,
                  getattr(args, 'min_length', 0), getattr(args, 'max_length', 1E10),
                  args.output,
                  args.verbosity, args.print_dest)

def read_index_table(index_table_file):
    unfiltered_index_table = pd.read_csv(index_table_file)
    return unfiltered_index_table

def filter_index_table(unfiltered_index_table, filters):
    # For each filter specified in args, filter the pd dataframe to only include the values given. 
    # In case of key error (i.e. filter in cmd line arg not present as a header in the csv) print informative error.
    filtered_index_table = unfiltered_index_table
    for filter in filters:
        try:
            filtered_index_table = unfiltered_index_table.loc[(unfiltered_index_table[filter["field"]].isin(filter["values"]))] 
        except:
            print("Check if csv has '{}' column.".format(filter["field"]))
    if len(filtered_index_table)==0:
        print("Warning: your filters have produced a dataframe of length 0.")
    return filtered_index_table

def process_files(input_file_or_directory, index_table, bins, filters,
                  min_length, max_length, output_prefix, verbosity, print_dest):
    """
    Core function to process one or more input files and create the required output files.

    Iterates through the reads in one or more input files and bins or filters them into the
    output files as required.
    """

    read_files = get_input_files(input_file_or_directory, verbosity, print_dest)

    if verbosity > 0:
        print('\n' + bold_underline('Processing reads in:'), flush=True, file=print_dest)
        for read_file in read_files:
            print(read_file, flush=True, file=print_dest)

    # A dictionary to store output files. These are created as needed when the first
    # read in each bin is being written. The prefix and suffix of the filenames are stored
    # here for convenience.
    out_files = {
        'prefix': output_prefix,
        'suffix': '.fastq'
    }

    counts = { 'read': 0, 'passed': 0, 'bins': {} }

    for read_file in read_files:

        file_type = get_sequence_file_type(read_file)

        if get_compression_type(read_file) == 'gz':
            open_func = gzip.open
        else:  # plain text
            open_func = open

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
                            write_read(out_files, filters, bins, min_length, max_length, name, sequence, None, counts,index_table)
                            sequence = ''
                        name = line[1:]
                    else:
                        sequence += line

                if name:
                    write_read(out_files, filters, bins, min_length, max_length, name, sequence, None, counts,index_table)

        else: # FASTQ
            with open_func(read_file, 'rt') as in_file:
                for line in in_file:
                    header = line.strip()[1:]
                    sequence = next(in_file).strip()
                    next(in_file) # spacer line
                    qualities = next(in_file).strip()

                    write_read(out_files, filters, bins, min_length, max_length, header, sequence, qualities, counts,index_table)

    for file in out_files:
        if hasattr(file, 'close'):
            file.close()

    if verbosity > 0:
        print('\nTotal reads read: ' + str(counts['read']), file=print_dest)
        print('\nReads written: ' + str(counts['passed']), file=print_dest)
        for file in counts['bins']:
            print(file + ': ' + str(counts['bins'][file]), file=print_dest)


def write_read(out_files, filters, bins, min_length, max_length, header, sequence, qualities, counts, index_table):
    '''
    Writes a read in either FASTQ or FASTA format depending if qualities are given
    :param out_file: The file to write to
    :param header:  The header line
    :param sequence: The sequence string
    :param qualities: The qualities string (None if FASTA)
    '''

    counts['read'] += 1

    length = len(sequence)

    if length > min_length and length < max_length:
        if index_table is not None:
            fields = get_csv_fields(index_table,header)
        
        else:
            fields = get_header_fields(header)
        
        if read_passes_filters(fields, index_table, filters):
            counts['passed'] += 1

            out_file = get_bin_output_file(fields, bins, out_files)
            if out_file:
                if qualities: # FASTQ
                    read_str = ''.join(['@', header, '\n', sequence, '\n+\n', qualities, '\n'])
                else:         # FASTA
                    read_str = ''.join(['>', header, '\n', sequence, '\n'])

                out_file.write(read_str)

                if not out_file.name in counts['bins']:
                    counts['bins'][out_file.name] = 0
                counts['bins'][out_file.name] += 1

    

def read_passes_filters(header_fields,index_table,filters):
    '''
    Returns true if the read passes all the filters
    :param fields:
    :param header_filters:
    :return:
    '''
    if index_table is not None:
        if filters:
            filtered_index_table = filter_index_table(index_table, filters)
            return header_fields["name"] in list(filtered_index_table["read"])
    else:
        for filter in filters:
    
            if filter['field'] in header_fields:
                if not header_fields[filter['field']] in filter['values']:
                    return False
            else:
                return False

    return True


def get_bin_output_file(fields, bins, out_files):
    '''
    This function decides which file to send this read to based on the current bins and filters. If
    the read's fields do not pass the filters then None is returned. Otherwise the file that is associated
    with the appropriate bin is returned. This will be created if it hasn't been already and the handle
    stored in `out_files`.

    :return: the appropriate output file if passes the filter, None if not
    # '''

    if bins:
        bin_name = ""

        for bin in bins:

            if bin in fields:
                bin_name += "_" + fields[bin]

        if len(bin_name) > 0:
            if not bin_name in out_files:
                out_files[bin_name] = open(out_files['prefix'] + bin_name + out_files['suffix'], "wt")

            return out_files[bin_name]

    if not 'unbinned' in out_files:
        out_files['unbinned'] = open(out_files['prefix'] + out_files['suffix'], "wt")

    return out_files['unbinned']


def get_csv_fields(index_table,header):
    name=header.split(' ')[0]
    row=index_table.query('read=="{}"'.format(name))
    fields = {'name': name}
    for i in row:
        fields[i]=list(row[i])[0]
    return fields


def get_header_fields(header):
    '''
    Splits a FASTA/FASTQ header line into key=value fields and returns them as
    an object.
    :param header:
    :return: An object of key value pairs
    '''

    parts = header.split(' ')

    fields = { 'name': parts[0] }

    for part in parts[1:]:
        (key, value) = part.split('=')
        fields[key] = value

    return fields


def get_input_files(input_file_or_directory, verbosity, print_dest):
    '''
    Takes a path to a single file or a directory and returns a list of file paths to be processed.
    :param input_file_or_directory: The input path
    :param verbosity: Verbosity level to report
    :param print_dest: Where to report (stdout or stderr)
    :return: An array of file paths to process
    '''

    input_files = []

    if os.path.isfile(input_file_or_directory):
        input_files.append(input_file_or_directory)

    # If the input is a directory, search it recursively for fastq files.
    elif os.path.isdir(input_file_or_directory):
        input_files = sorted([os.path.join(dir_path, f)
                              for dir_path, _, filenames in os.walk(input_file_or_directory)
                              for f in filenames
                              if f.lower().endswith('.fastq') or f.lower().endswith('.fastq.gz') or
                              f.lower().endswith('.fasta') or f.lower().endswith('.fasta.gz')])
        if not input_files:
            sys.exit('Error: could not find FASTQ/FASTA files in ' + input_file_or_directory)

    else:
        sys.exit('Error: could not find ' + input_file_or_directory)

    return input_files


def get_arguments():
    '''
    Parse the command line arguments.
    '''

    parser = argparse.ArgumentParser(description='BinLorry: a tool for binning sequencing reads into '
                                                 'files based on header information or read properties.',
                                     formatter_class=MyHelpFormatter, add_help=False)

    main_group = parser.add_argument_group('Main options')
    main_group.add_argument('-i', '--input', required=True,
                            help='FASTA/FASTQ of input reads or a directory which will be '
                                 'recursively searched for FASTQ files (required).')
    main_group.add_argument('-t', '--index-table', metavar='CSV_FILE', dest='index_table_file',
                           help='A CSV file with metadata fields for reads (otherwise these are assumed '
                                'to be in the read headers). This can also include a file and line number to '
                                'improve performance. Assumes read name is first column of the csv.')
    main_group.add_argument('-o', '--output', required=True,
                            help='Output filename (or filename prefix)')
    main_group.add_argument('-v', '--verbosity', type=int, default=1,
                            help='Level of output information: 0 = none, 1 = some, 2 = lots')

    bin_group = parser.add_argument_group('Binning/Filtering options')
    bin_group.add_argument('--bin-by', metavar='FIELD', nargs='+', dest='bin_by',
                            help='Specify header field(s) to bin the reads by. For multiple fields these '
                                 'will be nested in order specified. e.g. `--bin-by barcode reference`')
    bin_group.add_argument('--filter-by', metavar='FILTER', action='append', nargs='+', dest='filter_by',
                            help='Specify header field and accepted values to filter the reads by. Multiple '
                                 'instances of this option can be specified. e.g. `--filter-by reference BC01 '
                                 '--filter-by genotype Type1`')
    bin_group.add_argument('-n', '--min-length', metavar='MIN', type=int, dest='min_length',
                           help='Filter the reads by their length, specifying the minimum length.')
    bin_group.add_argument('-x', '--max-length', metavar='MAX', type=int, dest='max_length',
                           help='Filter the reads by their length, specifying the maximum length.')

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


