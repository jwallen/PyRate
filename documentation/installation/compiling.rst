*********************
Compiling from Source
*********************

If you prefer, you can also download a PyRate release package manually and
install that way. You might also wish to follow the bleeding-edge development
package.

Dependencies
============

PyRate depends on several other packages in order to provide its full
functional capabilities. The following dependencies are required to run PyRate:

* `Python <http://www.python.org/>`_ (version 2.5.x or later, including any version of Python 3, is recommended)

* `NumPy <http://numpy.scipy.org/>`_ (version 1.5.0 or later is recommended)

* `SciPy <http://www.scipy.org/>`_ (version 0.9.0 or later is recommended)

* `Quantities <http://packages.python.org/quantities/index.html>`_ (version 0.9.0 or later is recommended)

* `Cython <http://www.cython.org/>`_ (version 0.15 or later is recommended)

* A standard C compiler

Compiliation
============

Once you have obtained a copy of PyRate, either from a released package or
from cloning the git repository, PyRate can be installed by invoking the 
``setup.py`` script as usual::

$ python setup.py install

This will compile the Cython modules as part of the installation process, which
may take some time. (The speed boost gained during execution is well worth it!)
Note that you may need administrator privileges to install PyRate.

If you wish to use PyRate without installing, simply add the folder containing
this file to your ``PYTHONPATH`` environment variable and compile the Cython
modules in-place using the command ::

$ python setup.py build_ext --inplace

A Makefile that wraps these commands has also been provided.
