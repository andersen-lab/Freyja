Changelog
==========

All notable changes to this project will be documented in this file.

Versions
--------

- `Versions 2.0.x (2.0.0 – 2.0.3)`_
- `Versions 1.5.x (1.5.0 – 1.5.3)`_
- `Versions 1.4.x (1.4.1 – 1.4.9)`_
- `Versions 1.3.x (1.3.0 – 1.3.12)`_
- `Versions 1.2.x (1.2.0 – 1.2.1)`_
- `Versions 1.1.x (1.1.0)`_
- `Versions 1.0.x (1.0.0)`_

.. _Versions 2.0.x (2.0.0 – 2.0.3):

Version 2.0.3
-------------

Bug fixes
~~~~~~~~~

- Fixed Docker compatibility issues for pathogens other than SARS-CoV-2.
- Adjusted version reporting for ``freyja demix --version``.

Version 2.0.2
-------------

Improvements
~~~~~~~~~~~~

- Added ``covar`` support for accelerated ``freyja covariants`` analyses.
- Added barcode versioning for non-SARS-CoV-2 pathogens in the ``freyja-barcodes`` repository.
- Standardized barcode and lineage hierarchy nomenclature for non-SARS-CoV-2 pathogens.

Version 2.0.1
-------------

Improvements
~~~~~~~~~~~~

- Added control over the number of solver threads used by ``freyja demix``.
- Added verbose solver output options for ``freyja demix``.
- Added explicit checks for potential issues in lineage hierarchy YAML files.
- Updated amplicon statistics functionality, including BED-file validation for amplicon sequencing data.

Bug fixes
~~~~~~~~~

- Added minor bug fixes and documentation updates.

Version 2.0.0
-------------

Major changes
~~~~~~~~~~~~~

- Released Freyja 2.
- Expanded Freyja analyses beyond SARS-CoV-2 to support additional pathogens, including seasonal influenza, avian influenza, dengue, MPXV, and measles.

Improvements
~~~~~~~~~~~~

- Improved version checking.
- Updated barcode collapsing to support correct grouping for non-SARS-CoV-2 pathogens.
- Added documentation, including BarcodeForge documentation.

Bug fixes
~~~~~~~~~

- Added small bug fixes.

.. _Versions 1.5.x (1.5.0 – 1.5.3):

Version 1.5.3
-------------

Improvements
~~~~~~~~~~~~

- Improved flexibility and pandas robustness for ``freyja dash``.
- Added support for VCF outputs from ``bcftools`` and other variant callers.
- Added ``--autoadapt`` for automatic selection of the ``--adapt`` parameter.
- Added support for multiple new pathogens, including H5N1 avian influenza.

Bug fixes
~~~~~~~~~

- Fixed issues in ``freyja covariants``.

Version 1.5.2
-------------

Improvements
~~~~~~~~~~~~

- Added non-SARS-CoV-2 pathogen updating and demixing through the ``--pathogen`` flag.

Bug fixes
~~~~~~~~~

- Fixed bugs, including ``np.float()`` compatibility errors observed in some configurations.

Version 1.5.1
-------------

Improvements
~~~~~~~~~~~~

- Changed the default solver to ``clarabel``, with options to use ``ecos`` and ``osqp``.
- Improved ``freyja covariants``.
- Added plotting improvements.
- Switched barcodes to Feather format for faster loading and reduced file size.

Bug fixes
~~~~~~~~~

- Added small bug fixes.

Version 1.5.0
-------------

Improvements
~~~~~~~~~~~~

- Added ``--depthcutoff`` features to account for available genome coverage.
- Added de-aliasing and MRCA assignment functionality.
- Enabled both true and relaxed MRCA assignment.
- Incorporated collapsed MRCAs into summarization and plotting.
- Added small speed improvements.
- Added an option to set the ``variants`` threshold above zero to reduce output size.

.. _Versions 1.4.x (1.4.1 – 1.4.9):

Version 1.4.9
-------------

...

.. _Versions 1.3.x (1.3.0 – 1.3.12):

Version 1.3.12
--------------

...

.. _Versions 1.2.x (1.2.0 – 1.2.1):

Version 1.2.1
-------------

...

.. _Versions 1.1.x (1.1.0):

Version 1.1
-----------

...

.. _Versions 1.0.x (1.0.0):

Version 1.0
-----------

...