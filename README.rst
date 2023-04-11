Plot ensemble dust properties
=============================

Plot the ensemble properties of difference dust extinction samples.
Basic idea is to have one place to put together the A(V), R(V), E(B-V),
FM90 UV extinction parameters, and any other parameters for different
samples of dust extinction curves.  Examples of such samples are those
from Gordon et al. (2003) and Valencic et al. (2004).

In Development!
---------------

Active development.
Data still in changing.
Use at your own risk.

Contributors
------------
Karl Gordon

License
-------

This code is licensed under a 3-clause BSD style license (see the
``LICENSE`` file).

Sample Data
-----------

The sample data is given in an IPAC ASCII data table in the data directory.
These tables should have at least A(V), E(B-V), and R(V) values.

Where needed, code is provided to convert the original data tables to a common
format for this repository (e.g., `data/process_valencic04.py`).
