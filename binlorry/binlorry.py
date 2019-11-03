"""
Copyright 2019 Andrew Rambaut (a.rambaut@ed.ac.uk) & Áine O'Toole
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

from .misc import bold_underline, MyHelpFormatter, get_sequence_file_type, get_compression_type
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
        print(bold_underline("\nBins and filters:"),flush=True, file=args.print_dest)
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
            print("")

        if args.header_delimiters:
            print("Header delimiters ", args.header_delimiters)

    process_files(args.input_path, args.data_path, args.unordered_data,
                  args.bin_by, filters,
                  getattr(args, 'min_length', 0), getattr(args, 'max_length', 1E10),
                  getattr(args, 'header_delimiters', "="),
                  args.output, args.force_output,
                  args.verbosity, args.print_dest, args.out_report)

# def read_index_table(index_table_file):
#     unfiltered_index_table = pd.read_csv(index_table_file)
#     headers = unfiltered_index_table.columns.values
#     return unfiltered_index_table, headers

# def filter_index_table(unfiltered_index_table, filters, print_dest):
#     # For each filter specified in args, filter the pd dataframe to only include the values given.
#     # In case of key error (i.e. filter in cmd line arg not present as a header in the csv) print informative error.
#     for filter in filters:
#         try:
#             filtered_index_table = unfiltered_index_table.loc[(unfiltered_index_table[filter["field"]].isin(filter["values"]))]
#         except:
#             print("Check if csv has '{}' column.".format(filter["field"]))
#     return filtered_index_table

def process_files(input_path, input_report_path, unordered_data,
                  bins, filters, min_length, max_length, header_delimiters, output_prefix,
                  force_output, verbosity, print_dest, create_output_report):
    """
    Core function to process one or more input files and create the required output files.

    Iterates through the reads in one or more input files and bins or filters them into the
    output files as required.
    """

    read_files = get_input_files(input_path, verbosity, print_dest)

    if verbosity > 0:
        print(bold_underline('\nRead files found:'), flush=True, file=print_dest)
        for read_file in read_files:
            print(read_file, flush=True, file=print_dest)

    # A dictionary to store output files. These are created as needed when the first
    # read in each bin is being written. The prefix and suffix of the filenames are stored
    # here for convenience.
    output_files = {
        'prefix': output_prefix,
        'suffix': '.fastq' # todo - this should be fasta if appropriate
    }

    output_reports = None
    if create_output_report:
        output_reports = {
            'prefix': output_prefix,
            'suffix': '.csv'
        }

    if force_output:
        if len(bins) > 1:
            print(bold_underline("\nWarning: --force_output only works on a single bin/filter factor."))
            print(bold_underline("Usage: --force_output requires --bin_by to have matching --filter_by."))
        else:
            for filter in filters:
                for bin in bins:
                    if filter["field"] == bin:
                        print(bold_underline("\nForcing file creation of: "))
                        for value in filter["values"]:
                            filename = output_prefix + "_" + value + ".fastq"
                            with open(filename, "w"):
                                print(filename)
                                output_files[bin] = filename

                            if create_output_report:
                                report_name = filename.rstrip("fastq")+("csv")
                                with open(report_name, "w"):
                                    print(report_name)
                                    output_reports[bin] = report_name

    counts = { 'read': 0, 'passed': 0, 'bins': {} }

    report_dict = {}
    if input_report_path:
        report_dict = get_input_reports(read_files, input_report_path, verbosity, print_dest)

    # if '--unordered_data' is set then all the data tables should be loaded to be accessed as needed
    data_table = None
    if unordered_data:
        data_table = read_data_table(input_report_path, verbosity, print_dest)

    if verbosity > 0:
        print('\n' + bold_underline('Processing read file:'), flush=True, file=print_dest)

    for read_file in read_files:
        print(read_file, flush=True, file=print_dest)

        # is there a matching report file?
        report_file = report_dict[read_file]

        process_read_file(read_file, report_file, data_table, bins, filters, min_length, max_length,
                          header_delimiters, output_files, output_reports, counts, verbosity, print_dest)

    if output_files:
        for file in output_files:
            if hasattr(file, 'close'):
                file.close()

    if output_reports:
        for file in output_reports:
            if hasattr(file, 'close'):
                file.close()

    if verbosity > 0:
        print(bold_underline('\nSummary:')+'\nTotal reads read: ' + str(counts['read']), file=print_dest)
        print('Total reads written: ' + str(counts['passed']), file=print_dest)
        print(bold_underline('\nRead files written:'), file=print_dest)
        for file in counts['bins']:
            if output_reports:
                print(file + " with corresponding csv report file.", file=print_dest)
            else:
                print(file, file=print_dest)
        print("\n", file=print_dest)


def get_input_files(input_path, verbosity, print_dest):
    '''
    Takes a path to a single file or a directory and returns a list of file paths to be processed.
    :param input_file_or_directory: The input path
    :param verbosity: Verbosity level to report
    :param print_dest: Where to report (stdout or stderr)
    :return: An array of file paths to process
    '''
    input_files = []

    if os.path.isfile(input_path):
        input_files.append(input_path)

    # If the input is a directory, search it recursively for fastq files.
    elif os.path.isdir(input_path):
        input_files = sorted([os.path.join(dir_path, f)
                              for dir_path, _, filenames in os.walk(input_path)
                              for f in filenames
                              if f.lower().endswith('.fastq') or f.lower().endswith('.fastq.gz') or
                              f.lower().endswith('.fasta') or f.lower().endswith('.fasta.gz')])
        if not input_files:
            sys.exit('Error: could not find FASTQ/FASTA files in ' + input_path)

    else:
        sys.exit('Error: could not find ' + input_path)

    return input_files


def get_input_reports(input_files, input_report_path, verbosity, print_dest):
    '''
    returns a report dict with read file as key and report file as value (incl. paths).
    :param input_files: The array of input files
    :param input_report_path: The path to the report files
    :param verbosity: Verbosity level to report
    :param print_dest: Where to report (stdout or stderr)
    :return: An dictionary of report files
    '''
    report_dict = {}
    file_names = {}
    for i in input_files:
        file_names[i.split('/')[-1]] = i

    if len(input_files) == 1:
        if os.path.isfile(input_report_path):
            report_dict[file_names[0]] = input_report_path

    elif len(input_files) > 1:
        if os.path.isdir(input_report_path):

            for r,d,f in os.walk(input_report_path):
                for report_file in f:

                    for i in file_names:
                        file_stem = i.split('/')[-1].rstrip(".gz") \
                            .rstrip(".fastq").rstrip(".fasta") \
                            .rstrip(".FASTQ").rstrip(".FASTA")
                        if report_file.rstrip(".csv").rstrip(".CSV") == file_stem:
                            report_dict[file_names[i]]= r + '/' + report_file

    if verbosity > 0:
        print('\n' + bold_underline('Report files found:'), flush=True, file=print_dest)
        for i in sorted(report_dict):
            print(report_dict[i], flush=True, file=print_dest)

    return report_dict


def read_data_table(data_path, data_table, report_fields, print_dest):
    sys.exit('Unordered data files not implemented yet.')

    ### TODO implement reading of data tables into a dictionary keyed by read name
    data_table = {}

    line_number = 0

    with open_func(read_file, 'rt') as in_file:
        #report_fields = get_data_fields(in_data)

        for line in in_file:
            line_number += 1
            data, name = get_report_data(in_file, line, line_number, None, report_fields)

            data_table[name] = data

    # for read_file in read_files:
    #     try:
    #         report_file = report_dict[read_file]
    #
    #         full_index_table,headers = read_index_table(report_file)
    #
    #         index_table = filter_index_table(full_index_table, filters, print_dest)
    #         reads_with_reports[read_file] = index_table
    #     except:
    #         print("No corresponding report csv for {}!".format(read_file))
    #         continue
    #
    # print(bold_underline('\nHeaders:'), flush=True, file=print_dest)
    # print("\n".join(headers),flush=True, file=print_dest)

    return data_table


def process_read_file(read_file, report_file, data_table,
                      bins, filters, min_length, max_length, header_delimiters,
                      output_files, output_reports, counts, verbosity, print_dest):
    """
    Iterates through the reads in an input files and bins or filters them into the
    output files as required.
    """

    file_type = get_sequence_file_type(read_file)

    if get_compression_type(read_file) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open

    report_fields = None

    line_number = 0

    with open_func(read_file, 'rt') as in_file:
        in_report = None
        if report_file:
            in_report = open(report_file, 'rt')
            report_fields = get_report_fields(in_report)

            if 'fields' in output_reports:
                # todo check that the fields are the same
                pass
            else:
                output_reports['fields'] = report_fields

        if file_type == 'FASTA':
            # For FASTA files we need to deal with line wrapped sequences...

            header = ''
            sequence = ''

            for line in in_file:
                line = line.strip()

                if not line:
                    continue

                if line[0] == '>':  # Header line = start of new read
                    if header:
                        line_number += 1
                        data = get_read_data(header, header_delimiters, report_fields, in_report, data_table, line_number)
                        filter_bin_read(output_files, output_reports, filters, bins, min_length, max_length, data, header, sequence, None, counts, print_dest)
                        sequence = ''
                    header = line[1:]
                else:
                    sequence += line

            if header:
                line_number += 1
                data = get_read_data(header, header_delimiters, report_fields, in_report, data_table, line_number)
                filter_bin_read(output_files, output_reports, filters, bins, min_length, max_length, data, header, sequence, None, counts, print_dest)

        else: # FASTQ
            for line in in_file:
                header = line.strip()[1:]
                sequence = next(in_file).strip()
                next(in_file) # spacer line
                qualities = next(in_file).strip()

                line_number += 1
                data, name = get_read_data(header, header_delimiters, report_fields, in_report, data_table, line_number)
                filter_bin_read(output_files, output_reports, filters, bins, min_length, max_length, data, header, sequence, qualities, counts, print_dest)

        in_file.close()

    if output_files:
        for file in output_files:
            if hasattr(file, 'close'):
                file.close()


def get_report_fields(in_data):

    for line in in_data:
        line = line.strip()

        if not line:
            continue

        return line.split(',')


def get_read_data(header, header_delimiters, report_fields, in_data, data_table, line_number):

    name = header.split(' ')[0]

    if data_table:
        return data_table[name], name

    elif in_data:
        for line in in_data:
            line = line.strip()
            if not line:
                continue
            else:
                return get_report_data(in_data, line, line_number, name, report_fields)

    else:
        return get_header_data(header, header_delimiters), name


def get_report_data(file, line, line_number, name, report_fields):
    data = {}
    values = line.split(',')

    if len(values) != len(report_fields):
        sys.exit('input report file, ' + file.name + ', line ' + str(line_number + 1) + ', has the wrong number of fields')

    if name and values[0] != name:
        sys.exit("input report file, " + file.name + ", line " + str(line_number + 1) +
                 ", the first column doesn't match the name of the read")

    for i in range(len(values)):
        data[report_fields[i]] = values[i]

    return data, values[0]


def filter_bin_read(output_files, output_reports, filters, bins, min_length, max_length, data,
                    header, sequence, qualities, counts, print_dest):
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

        passed = read_passes_filters(data, filters)

        if passed:
            counts['passed'] += 1
            out_file = get_bin_output_file(data, bins, output_files, False)

            out_report = None
            if output_reports:
                out_report = get_bin_output_file(data, bins, output_reports, True)
                pass
            # todo - create output reports
            # col_names = ",".join(list(filtered_index_table.columns.values))
            # out_report = get_bin_output_file(data, bins, out_reports, col_names, True)
            #
            # row = filtered_index_table.loc[filtered_index_table["read_name"] == fields["read_name"]].values
            # row_str = ",".join([str(i) for i in list(row[0])])
            # out_report.write(row_str + "\n")

            if out_file:
                write_read(out_file, header, sequence, qualities)

                if out_report:
                    out_report.write(",".join(data.values()) + "\n")

                if not out_file.name in counts['bins']:
                    counts['bins'][out_file.name] = 0
                counts['bins'][out_file.name] += 1


def write_read(out_file, header, sequence, qualities):
    '''
    Writes a read in either FASTQ or FASTA format depending if qualities are given
    :param out_file: The file to write to
    :param header:  The header line
    :param sequence: The sequence string
    :param qualities: The qualities string (None if FASTA)
    '''

    if qualities: # FASTQ
        read_str = ''.join(['@', header, '\n', sequence, '\n+\n', qualities, '\n'])
    else:         # FASTA
        read_str = ''.join(['>', header, '\n', sequence, '\n'])

    out_file.write(read_str)


def read_passes_filters(header_fields, filters):
    '''
    Returns true if the read passes all the filters
    :param fields:
    :param header_filters:
    :return:
    '''
    for filter in filters:

        if filter['field'] in header_fields:
            if not header_fields[filter['field']] in filter['values']:
                return False
        else:
            return False

    return True


def  get_bin_output_file(fields, bins, output_files, report):
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
            if not bin_name in output_files:
                output_files[bin_name] = open(output_files['prefix'] + bin_name + output_files['suffix'], "wt")
                if report:
                    output_files[bin_name].write(",".join(output_files['fields']) +'\n')
            return output_files[bin_name]

    if not 'unbinned' in output_files:
        output_files['unbinned'] = open(output_files['prefix'] + output_files['suffix'], "wt")
        if report:
            output_files["unbinned"].write(",".join(output_files['fields']) +'\n')
    return output_files['unbinned']


def get_csv_fields(index_table,header):
    fields = {}
    name=header.split(' ')[0]
    try:
        row=index_table.query('read_name=="{}"'.format(name))
        fields = {'name': name}
        for i in row:
            fields[i]=list(row[i])[0]
    except:
        pass


    return fields

def get_header_data(header, header_delimiters):
    '''
    Splits a FASTA/FASTQ header line into key=value fields and returns them as
    an object.
    :param header:
    :return: An object of key value pairs
    '''

    parts = header.split(' ')

    fields = { 'name': parts[0] }

    for part in parts[1:]:
        for delimiter in header_delimiters:
            try:
                (key, value) = part.split(delimiter)
                fields[key] = value
            except:
                pass

    return fields

def get_arguments():
    '''
    Parse the command line arguments.
    '''
    parser = argparse.ArgumentParser(description='BinLorry: a tool for binning sequencing reads into '
                                                 'files based on header information or read properties.',
                                     formatter_class=MyHelpFormatter, add_help=False)

    main_group = parser.add_argument_group('Main options')
    main_group.add_argument('-i', '--input', dest='input_path', required=True,
                            help='FASTA/FASTQ of input reads or a directory which will be '
                                 'recursively searched for FASTQ files (required).')
    main_group.add_argument('-t', '--data', metavar='CSV_FILE', dest='data_path', default=None,
                            help='A CSV file with metadata fields for reads or a directory of csv files that will '
                                 'be recursively searched for names corresponding to a matching input FASTA/FASTQ '
                                 'files.')
    main_group.add_argument('-u', '--unordered_data', action='store_true',
                            help='The metadata tables are not in the same order as the reads - they will all be'
                                 'loaded and then looked up as needed (slower).')
    main_group.add_argument('-o', '--output', required=True,
                            help='Output filename (or filename prefix)')
    main_group.add_argument('-r', '--out-report',action='store_true',dest="out_report",
                            help='Output a report along with FASTA/FASTQ.')
    main_group.add_argument('-f', '--force-output',action='store_true',dest="force_output",
                            help='Output binned/filtered files even if empty.')
    main_group.add_argument('-v', '--verbosity', type=int, default=1,
                            help='Level of output information: 0 = none, 1 = some, 2 = lots')

    bin_group = parser.add_argument_group('Binning/Filtering options')
    bin_group.add_argument('--bin-by', metavar='FIELD', nargs='+', dest='bin_by',
                           help='Specify header field(s) to bin the reads by. For multiple fields these '
                                'will be nested in order specified. e.g. `--bin-by barcode reference`')
    bin_group.add_argument('--filter-by', metavar='FILTER', action='append', nargs='+', dest='filter_by',
                           help='Specify header field and accepted values to filter the reads by. Multiple '
                                'instances of this option can be specified. e.g. `--filter-by barcode BC01 BC02'
                                '--filter-by genotype Type1`')
    bin_group.add_argument('-n', '--min-length', metavar='MIN', type=int, dest='min_length',
                           help='Filter the reads by their length, specifying the minimum length.')
    bin_group.add_argument('-x', '--max-length', metavar='MAX', type=int, dest='max_length',
                           help='Filter the reads by their length, specifying the maximum length.')
    bin_group.add_argument('-d', '--header-delimiters', metavar='DELIM', dest='header_delimiters', default="=",
                           help='Delimiters to use when searching for key:value pairs in FASTA/FASTQ header.')


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


