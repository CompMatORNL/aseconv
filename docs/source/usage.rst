Usage
=====

.. _installation:

Installation
------------

To use aseconv, first install it using pip:

.. code-block:: console

   (.venv) $ pip install compmat-aseconv


.. _examples:

Examples
--------
    
- Convert a `FHI-aims`_ `1_nacl.in` file to `VASP`_ poscar file in Cartesian format after sorting atom orders by element names. Also, add constraint to atoms whose z coordinates are less than 2.5.

.. code-block:: console

    aseconv geo -t vasp --elsort -C --con "z<2.5" 1_nacl.in

    >>> Globbing '1_nacl.in' ...
     > Processing '1_nacl.in'...
     - Writing [Cl4Na4] '1_nacl_Cz2.5.poscar'...
 
**Input**: :download:`1_nacl.in <../../examples/1_nacl.in>`

.. literalinclude:: ../../examples/1_nacl.in

**Output**: `1_nacl_Cz2.5.poscar`

.. code-block:: text

    Na Cl
     1.0000000000000000
         5.6916940000000000    0.0000000000000000    0.0000000000000000
         0.0000000000000000    5.6916940000000000    0.0000000000000000
         0.0000000000000000    0.0000000000000000    5.6916940000000000
     Na  Cl
       4   4
    Selective dynamics
    Cartesian
      0.0000000000000000  0.0000000000000000  0.0000000000000000   F   F   F
      0.0000000000000000  2.8458470000000000  2.8458470000000000   T   T   T
      2.8458470000000000  0.0000000000000000  2.8458470000000000   T   T   T
      2.8458470000000000  2.8458470000000000  0.0000000000000000   F   F   F
      2.8458470000000000  0.0000000000000000  0.0000000000000000   F   F   F
      2.8458470000000000  2.8458470000000000  2.8458470000000000   T   T   T
      0.0000000000000000  0.0000000000000000  2.8458470000000000   T   T   T
      0.0000000000000000  2.8458470000000000  0.0000000000000000   F   F   F

.. _help:

Help
----

.. argparse::
   :filename: ../src/aseconv/main.py
   :func: acmparser
   :prog: aseconv
   :markdownhelp:

.. _VASP: https://www.vasp.at/
.. _FHI-aims: https://fhi-aims.org/