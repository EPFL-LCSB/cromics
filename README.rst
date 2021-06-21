CROMICS

Papers
1. Angeles-Martinez and Hatzimanikatis. "Spatio-temporal modeling of the crowding conditions and metabolic variability in microbial communities".
2. Angeles-Martinez and Hatzimanikatis. "The influence of the crowding assumptions in biofilm simulations".
=======

Requirements
------------

You will need to have `Git-LFS <https://git-lfs.github.com/>`_ in order to properly download some binary files:

.. code:: bash

    git clone https://github.com/EPFL-LCSB/cromics.git /path/to/cromics
    cd /path/to/cromics
    git lfs install
    git lfs pull

The scripts have been developed with Matlab 2018b, and CPLEX 12.7.1 (freely downloadable with the `IBM Academic initiative <https://developer.ibm.com/academic/>`_) and successfully ran on several other versions of both softwares. However, it is important to respect the IBM compatibility specs sheets between Matlab, CPLEX, and the computer OS - available `on IBM's website <https://www.ibm.com/software/reports/compatibility/clarity/index.html>`_

This module requires matTFA.


License
=======
The software in this repository is put under an APACHE licensing scheme - please see the `LICENSE <https://github.com/EPFL-LCSB/cromics/blob/master/LICENSE>`_ file for more details.

This software uses open source components. 
