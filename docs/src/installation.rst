Installation
-------------------------------------------------------------------------------

Freyja is entirely written in Python 3, but requires preprocessing by tools like iVar and `samtools <https://github.com/samtools/samtools>`_ mpileup to generate the required input data. We recommend using python3.7, but Freyja has been tested on python versions up to 3.10.

Install via Conda::

    conda install -c bioconda freyja


Local build from source::

    git clone https://github.com/andersen-lab/Freyja.git
    cd Freyja
    pip install -e .

Please note that pip will not install all the requirements.
You need to install some requirements via Conda by following these steps::

    conda create -n freyja
    conda activate freyja
    conda env update --file environment.yml

Docker::

    docker pull staphb/freyja
    docker run --rm -it staphb/freyja [command]

        