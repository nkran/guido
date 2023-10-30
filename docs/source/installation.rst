Installation
============

Install guido and azimuth (dependency) via PyPi using pip:

.. code-block:: bash

    $ pip install guido
    $ pip install git+https://github.com/Biomatters/Azimuth

Please note that guido requires bowtie to be installed on your system. Please follow the instructions on the `bowtie website
<http://bowtie-bio.sourceforge.net/index.shtml>`_  to install it.

Additionaly, guido requires tabix to be installed in order to search for off-targets.
Please follow the instructions on the `htslib website <https://github.com/samtools/htslib>`_ to install it.

Alternatively, you can install both bowtie and tabix by using conda:

.. code-block:: bash

    $ conda install -c bioconda bowtie
    $ conda install -c bioconda tabix
