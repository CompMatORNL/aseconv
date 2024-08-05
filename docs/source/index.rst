.. CompMat - aseconv documentation master file, created by
   sphinx-quickstart on Wed Jul 17 13:24:53 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   H3(=) H4(-) H5(^) H6(")

#################
CompMat - aseconv
#################

**aseconv** is a Python CLI utility for geometry manipulation using ASE_ library. Check out the :doc:`usage` section for further information.

Features
========

service
-------

 Launch/stop a background processing server for many file processing to avoid initial python loading time. 
 
geo
---
 - Geometry manipulation functions. 
 
    - Geometry file type conversion: Supports the transformation of geometry files into different formats.
    - Extended file type support: Our code introduces support for additional file types and features beyond those supported by the ASE_ library.
    - Coordinate system conversion: Cartesian <-> Fractional.
    - Atom position translation
    - Rotation of axes
    - Supercell construction
    - Axis alignment: Users can align any axis of the structure to a target direction.
    - Surface generation: The code can generate surfaces based on Miller indices (hkl).
    - Strain implementation: Introduce strain into structures.
    - Consecutive operations are applied in the order of arguments.
    
 - K-path generator. 
    
    Using the SeeK-path_, generate `VASP`_/`Wannier90`_/`WannierTools`_ or `FHI-aims`_ 3D k-path input strings.
    It also automatically detects 2D slab structures, and generate 2D K-path.


Plugins
=======

 **aseconv** suports plugins. There are two types base plugin classes, :py:class:`~.pluginbase.AsecPlug` and :py:class:`~.pluginbase.AsecIO`.
 The :py:class:`~.pluginbase.AsecPlug` is for atom manipluation, and :py:class:`~.pluginbase.AsecIO` is for 
 read/write function for different file formats. If the name of a class starts with 'NoPlug', then the plugin is not loaded.
 The plugin files which contain plugin classes must be named as 'plug*.py' or 'io*.py'
 and should be placed in the plugin path set in an environmetal variable ``ASEC_PLUGIN_PATH``.
 Multiple plugin paths can be set in the variable with a delimiter `:`.
 
 One plugin class needs to register at least one argument to be passed. And the :py:mod:`~.pluginbase.AsecPlug.process`
 function must be defined.
 
- Example of :py:class:`~.pluginbase.AsecPlug` plugin.

 
.. literalinclude:: ../../src/aseconv/plugins/plugAbasic.py
   :pyobject: APlugSelect


- Example of :py:class:`~.pluginbase.AsecIO` plugin.

   For this plugin, if the read function is not defined, :py:mod:`ase.io.read` is used, instead.
   
   
.. literalinclude:: ../../src/aseconv/plugins/iormg.py
   :pyobject: RMGIO

 
 
 
.. toctree::
   :maxdepth: 4
   :caption: Contents:

   usage
   modules
   

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. _ASE: https://wiki.fysik.dtu.dk/ase/index.html
.. _VASP: https://www.vasp.at/
.. _wannier90: https://wannier.org/
.. _WannierTools: http://www.wanniertools.com/
.. _SeeK-path: https://seekpath.readthedocs.io/en/latest/maindoc.html
.. _FHI-aims: https://fhi-aims.org/
