Creating custom barcodes
-------------------------------------------------------------------------------

Follow these steps to generate lineage‑specific barcodes with **BarcodeForge**.

1. Install **Nextflow** >=24.04.2 (see https://www.nextflow.io/docs/latest/install.html for installation instructions) and **Mamba/Conda** (see https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html).

1. Prepare a lineage tree

   Build a phylogenetic tree that includes every lineage you want covered in the barcode set
   (e.g. ``tree.nwk``).

2. Populate ``params.json``

    Configure the input and output paths. Adjust the file locations and any settings as needed before running the pipeline:

    .. code-block:: json

        {
          "alignment":        "data/aligned.fasta",
          "reference_genome": "data/reference.fasta",
          "tree_file":        "data/tree.nwk",
          "tree_file_format": "newick",
          "lineages":         "data/lineages.tsv",
          "barcode_prefix":   "RSVa",
          "outdir":           "results"
        }

    .. note::

       The parameters used above are described as follows:

       - alignment: Path to the aligned FASTA file required for barcode generation.
       - reference_genome: Path to the reference genome file used for alignment.
       - tree_file: Path to the phylogenetic tree file (supports Newick and Nexus formats).
       - tree_file_format: Format of the tree file ("newick" or "nexus").
       - lineages: Path to the TSV file containing lineage definitions.
       - barcode_prefix: A prefix string added to all generated barcodes (e.g. <prefix>-<lineage>).
       - outdir: Directory where the pipeline output (barcodes and logs) will be stored.

3. Run the pipeline

    .. important::
        
        The pipeline requires a working installation of **Nextflow** (see https://www.nextflow.io/docs/latest/install.html for installation instructions) and **Mamba/Conda**.

   .. code-block:: bash

      nextflow run https://github.com/andersen-lab/BarcodeForge.git \
        -profile <conda/mamba> \
        -params-file <path to the params.json file> \
        -latest

    .. note::
        The pipeline can use the following profiles:
        - conda: Uses Conda for package management.
        - mamba: Uses Mamba for package management.

        If you wish to use a specific version of the pipeline, replace ``-latest`` with the desired version tag (e.g. ``-r "v1.0.0"``).
        More information about ``-r`` can be found here: https://www.nextflow.io/docs/latest/cli.html#using-a-specific-revision

4. Retrieve the output

   The pipeline writes results to ``results/barcode/``:

   * ``barcodes.csv`` – barcode definitions for each lineage  
   * ``barcodes.html`` – the same barcodes in an interactive HTML format


Example
=======
The following example shows how to generate barcodes for the RSV-A lineage tree:

1. Download the following files to a folder named ``data``:

    .. code-block:: bash
        mkdir data
        wget https://raw.githubusercontent.com/andersen-lab/BarcodeForge/refs/heads/main/assets/test/input/tree.nwk -O data/tree.nwk
        wget https://raw.githubusercontent.com/andersen-lab/BarcodeForge/refs/heads/main/assets/test/input/aligned.fasta -O data/aligned.fasta
        wget https://raw.githubusercontent.com/andersen-lab/BarcodeForge/refs/heads/main/assets/test/input/reference.fasta -O data/reference.fasta
        wget https://raw.githubusercontent.com/andersen-lab/BarcodeForge/refs/heads/main/assets/test/input/lineages.tsv -O data/lineages.tsv

    - `RSV-A tree <https://github.com/andersen-lab/BarcodeForge/blob/main/assets/test/input/tree.nwk>`_
    - `RSV-A alignment <https://github.com/andersen-lab/BarcodeForge/blob/main/assets/test/input/aligned.fasta>`_
    - `RSV-A reference genome <https://github.com/andersen-lab/BarcodeForge/blob/main/assets/test/input/reference.fasta>`_
    - `RSV-A lineages per sample <https://github.com/andersen-lab/BarcodeForge/blob/main/assets/test/input/lineages.tsv>`_

2. Create a file named ``params.json`` **outside** the ``data`` folder with the following content:

    .. code-block:: json

        {
          "alignment":        "data/aligned.fasta",
          "reference_genome": "data/reference.fasta",
          "tree_file":        "data/tree.nwk",
          "tree_file_format": "newick",
          "lineages":         "data/lineages.tsv",
          "barcode_prefix":   "RSVa",
          "outdir":           "results"
        }

3. Run the pipeline:
    .. code-block:: bash

        nextflow run https://github.com/andersen-lab/BarcodeForge.git \
            -profile mamba \
            -params-file params.json \
            -latest

    .. note::
        If you wish to use a specific version of the pipeline, replace ``-latest`` with the desired version tag (e.g. ``-r "v1.0.0"``).
        More information about ``-r`` can be found here: https://www.nextflow.io/docs/latest/cli.html#using-a-specific-revision

4. Retrieve the output:
    The pipeline writes results to ``results/barcode/``:
    
    * ``barcodes.csv`` – barcode definitions for each lineage  
    * ``barcodes.html`` – the same barcodes in an interactive HTML format
