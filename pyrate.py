#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   PyRate - Python tools for computing chemical reaction rates
#
#   Copyright (c) 2012 by Joshua W. Allen (jwallen@mit.edu)
#                         Yury V. Suleimanov (ysuleyma@mit.edu)
#                         William H. Green (whgreen@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a 
#   copy of this software and associated documentation files (the "Software"), 
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
#   and/or sell copies of the Software, and to permit persons to whom the 
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
#   DEALINGS IN THE SOFTWARE. 
#
################################################################################

"""
This is the main executable script for PyRate, a tool for computing chemical
reaction rates and other properties used in detailed kinetics models using
various methodologies and theories. To run PyRate, use the command ::

    $ python pyrate.py FILE

where ``FILE`` is the path to a PyRate input file describing the job to
execute. PyRate will run the specified job, writing the output to ``output.py``
and a log to both the console and to ``pyrate.log``, with both files appearing
in the same directory as the input file. Some additional command-line
arguments are available; run the command ::

    $ python pyrate.py -h

for more information.
"""

import os
import os.path
import argparse

################################################################################

def parseCommandLineArguments():
    """
    Parse the command-line arguments being passed to PyRate. This uses the
    :mod:`argparse` module, which ensures that the command-line arguments are
    sensible, parses them, and returns them.
    """

    parser = argparse.ArgumentParser(description=
    """
    PyRate is a tool for computing chemical reaction rates and other
    properties used in detailed kinetics models using various methodologies
    and theories.
    """)
    parser.add_argument('file', metavar='FILE', type=str, nargs=1,
        help='a file describing the job to execute')

    # Options for controlling the amount of information printed to the console
    # By default a moderate level of information is printed; you can either
    # ask for less (quiet), more (verbose), or much more (debug)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-q', '--quiet', action='store_const', const=logging.WARNING, default=logging.INFO, dest='verbose', help='only print warnings and errors')
    group.add_argument('-v', '--verbose', action='store_const', const=logging.DEBUG, default=logging.INFO, dest='verbose', help='print more verbose output')
    group.add_argument('-d', '--debug', action='store_const', const=0, default=logging.INFO, dest='verbose', help='print debug information')

    # Add options for controlling what directories files are written to
    parser.add_argument('-o', '--output-directory', type=str, nargs=1, default='',
        metavar='DIR', help='use DIR as output directory')

    return parser.parse_args()

################################################################################

if __name__ == '__main__':
    
    # Parse and validate the command-line arguments
    args = parseCommandLineArguments()
    
    # Determine the output directory
    # By default the directory containing the input file is used, unless an
    # alternate directory is specified using the -o flag
    if args.output_directory and os.path.isdir(args.output_directory[0]):
        outputDirectory = os.path.abspath(args.output_directory[0])
    else:
        outputDirectory = os.path.dirname(os.path.abspath(args.file[0]))
