#!/usr/bin/env python3
"""
This code is adapted from Porechop by Ryan Wick. The original licensing information is as follows:

Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Porechop

This file is part of Porechop. Porechop is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Porechop is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Porechop. If
not, see <http://www.gnu.org/licenses/>.
"""

import os
import shutil
import subprocess
import sys

from .misc import load_fasta_or_fastq, bold_underline, int_to_str


class Read(object):

    def __init__(self, name, seq, quals):
        self.name = name
        self.seq = seq
        self.quals = quals
        if len(quals) < len(seq):
            self.quals += '+' * (len(seq) - len(quals))

        self.reference_scores = {}

        self.best_reference = ('none', 0.0)

    def get_name(self, reference_labels=False):

        if reference_labels:
            return self.name + " reference=" + self.best_reference[0] + " ref_dist=" + self.best_reference[1]

        return self.name

    def get_seq(self):
        return self.seq

    def seq_length(self):
        return len(self.get_seq())

    def get_fasta(self, reference_labels=False):
        if not self.seq:  # Don't return empty sequences
            return ''
            
        return ''.join(['>', self.get_name(reference_labels), '\n', self.seq])

    def get_fastq(self, reference_labels=False):
        if not self.seq:  # Don't return empty sequences
            return ''

        return ''.join(['@', self.get_name(reference_labels), '\n', self.seq, '\n+\n', self.quals, '\n'])

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
            reads = [Read(x[2], x[1], '') for x in reads]
        else:  # FASTQ
            reads = [Read(x[4], x[1], x[3]) for x in reads]
        check_reads = reads[:check_read_count]

    # If the input is a directory, search it recursively for fastq files.
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
            file_reads = [Read(x[4], x[1], x[3]) for x in file_reads]

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


def output_reads(reads, out_format, output, read_type, verbosity,
                 reference_labels, print_dest, input_filename,
                threads, discard_unmapped):
    if verbosity > 0:

        if output is None:
            verb = 'Outputting '
            destination = 'stdout'
        else:
            verb = 'Saving '
            destination = 'file'
        print(bold_underline(verb + 'mapped reads to ' + destination),
              flush=True, file=print_dest)

    if out_format == 'auto':
        if output is None:
            out_format = read_type.lower()
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
    if out_format.endswith('.gz') and (output is not None):
        gzipped_out = True
        out_format = out_format[:-3]
        if shutil.which('pigz'):
            if verbosity > 0:
                print('pigz found - using it to compress instead of gzip')
            gzip_command = 'pigz -p ' + str(threads)
        else:
            if verbosity > 0:
                print('pigz not found - using gzip to compress')

    # Output to all reads to stdout.
    if output is None:
        for read in reads:
            read_str = read.get_fasta(reference_labels=reference_labels) if out_format == 'fasta' \
                else read.get_fastq(reference_labels=reference_labels)
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
                read_str = read.get_fasta(reference_labels=reference_labels) if out_format == 'fasta' \
                    else read.get_fastq(reference_labels=reference_labels)
                out.write(read_str)
        if gzipped_out:
            subprocess.check_output(gzip_command + ' -c ' + out_filename + ' > ' + output,
                                    stderr=subprocess.STDOUT, shell=True)
            os.remove(out_filename)
        if verbosity > 0:
            print('\nSaved result to ' + os.path.abspath(output), file=print_dest)

    if verbosity > 0:
        print('', flush=True, file=print_dest)

