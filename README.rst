***********************************************************
PyRate - Python tools for computing chemical reaction rates
***********************************************************

About PyRate
============

PyRate is a free, open-source Python toolkit for computing chemical reaction
rates using various methods and theories.

License
=======

PyRate is distributed under the terms of the 
`MIT license <http://www.opensource.org/licenses/mit-license>`_::

    Copyright (c) 2012 by Joshua W. Allen (jwallen@mit.edu)
                          Yury V. Suleimanov (ysuleyma@mit.edu)
                          William H. Green (whgreen@mit.edu)
    
    Permission is hereby granted, free of charge, to any person obtaining a 
    copy of this software and associated documentation files (the "Software"), 
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense, 
    and/or sell copies of the Software, and to permit persons to whom the 
    Software is furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
    THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
    DEALINGS IN THE SOFTWARE. 

Installation
============

Dependencies
------------

PyRate depends on several other packages in order to provide its full
functional capabilities. The following dependencies are required to run PyRate:

* `Python <http://www.python.org/>`_ (version 2.5.x or later, including any version of Python 3, is recommended)

* `NumPy <http://numpy.scipy.org/>`_ (version 1.5.0 or later is recommended)

* `SciPy <http://www.scipy.org/>`_ (version 0.9.0 or later is recommended)

* `Quantities <http://packages.python.org/quantities/index.html>`_ (version 0.9.0 or later is recommended)

* `Cython <http://www.cython.org/>`_ (version 0.15 or later is recommended)

* A standard C compiler

Some features require the availability of an external quantum chemistry package.
At least one of the following quantum chemistry packages must be installed to
use these features:

* `Gaussian <http://www.gaussian.com/>`_ (version 03 or later is recommended)

Getting PyRate
--------------

The best way to obtain a copy of the repository is to clone it using `git
<http://git-scm.com/>`_::

    $ git clone git://github.com/GreenGroup/PyRate.git

This enables you to easy update your local clone with the latest changes. If
you intend to contribute, you should fork the project on GitHub first, and
clone from your fork.

Compiling from Source
---------------------

PyRate is installed by invoking the ``setup.py`` script::

    $ python setup.py install

This will compile the Cython modules as part of the installation process, which
may take some time. (The speed boost gained during execution is well worth it!)
Note that you may need administrator privileges to install PyRate.

If you wish to use PyRate without installing, simply add the folder containing
this file to your ``PYTHONPATH`` environment variable and compile the Cython
modules in-place using the command ::

    $ python setup.py build_ext --inplace

A Makefile that wraps these commands has been provided. The Makefile also
provides a clean target for deleting all files created during the in-place
compilation.

Running the Unit Tests
----------------------

PyRate comes with a large suite of unit tests that ensure functionality is
working as intended. To run these tests, first install PyRate or compile it
from source in-place using the directions above. Then, simply invoke the entire
suite of unit tests using ::

    $ python setup.py test

This will run all of the unit tests in sequence, which may take some time. If
one or more unit tests fail, please report them on the GitHub issue tracker
(if they aren't already reported).
