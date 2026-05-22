Freyja Changelog
================

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

Improvements
~~~~~~~~~~~~

- Improved error reporting to make failure modes easier to interpret.
- Updated ``dash`` and ``plot`` to prevent labels and legends from being cut off.
- Added spike or arbitrary genome-region coverage estimates.
- Added options to retain intermediate outputs from ``dash`` and ``plot``.
- Added lineage hierarchy input support for ``demix``, ``plot``, ``dash``, and ``relgrowthrate``.

Version 1.4.8
-------------

Improvements
~~~~~~~~~~~~

- Added beta adaptive-regularization functionality.

Bug fixes
~~~~~~~~~

- Added bug fixes.

Version 1.4.7
-------------

Bug fixes
~~~~~~~~~

- Fixed an issue with ``freyja update``.

Version 1.4.6
-------------

Improvements
~~~~~~~~~~~~

- Updated barcode generation to pull from the GISAID tree by default.
- Added an option to include amino-acid annotations in ``variants``.
- Updated ``plot-covariants`` to create heatmap plots showing the frequency of each covariant pattern.

Bug fixes
~~~~~~~~~

- Added additional minor fixes.

Version 1.4.5
-------------

Improvements
~~~~~~~~~~~~

- Added lineage collapse functionality to account for sequencing coverage using ``--depthcutoff`` in ``freyja demix``.

Bug fixes
~~~~~~~~~

- Added small fixes to ``covariants``.

Version 1.4.4
-------------

Bug fixes
~~~~~~~~~

- Fixed ``plot`` issues.
- Fixed metadata parsing issues.

Improvements
~~~~~~~~~~~~

- Refined methods for ``covariants`` and associated plotting in ``plot-covariants``.

Version 1.4.3
-------------

Improvements
~~~~~~~~~~~~

- Added beta ``barcode_build`` functionality.

Bug fixes
~~~~~~~~~

- Fixed issues in ``freyja covariants``.
- Added small bug fixes for ``freyja plot``.

Version 1.4.2
-------------

Bug fixes
~~~~~~~~~

- Fixed version reporting for ``freyja --version``.

Version 1.4.1
-------------

Improvements
~~~~~~~~~~~~

- Added read-analysis tools for mutation cluster detection and visualization.
- Added mutation-specific read extraction and filtering.
- Added barcode version checking through ``freyja demix --version``.

.. _Versions 1.3.x (1.3.0 – 1.3.12):

Version 1.3.12
--------------

Improvements
~~~~~~~~~~~~

- Switched updates to download from the repository by default.
- Added weekly plot intervals.
- Reversed legend order to better match plot order.
- Added support for lineages not observed in public data, including GISAID-only lineages.

Version 1.3.11
--------------

Improvements
~~~~~~~~~~~~

- Added relative growth-rate functionality.
- Improved ``freyja dash``.

Version 1.3.10
--------------

Improvements
~~~~~~~~~~~~

- Moved lineage metadata downloading into ``freyja update``.
- Defaulted to the ``cov-lineages`` lineage list to handle asynchronous lineage lists between UShER and ``cov-lineages``.
- Pulled lineage data from a patched ``lineage.yml`` from outbreak.info.

Version 1.3.9
-------------

Improvements
~~~~~~~~~~~~

- Added custom lineage grouping, colors, and ordering in ``freyja dash``.

Bug fixes
~~~~~~~~~

- Fixed dashboard bugs.

Version 1.3.8
-------------

Improvements
~~~~~~~~~~~~

- Improved robustness of ``freyja dash``.

Version 1.3.7
-------------

Improvements
~~~~~~~~~~~~

- Added dashboard output functionality.

Version 1.3.6
-------------

Improvements
~~~~~~~~~~~~

- Added sequencing coverage estimates.
- Added the ``--confirmedonly`` option to ``freyja demix``.

Version 1.3.5
-------------

Improvements
~~~~~~~~~~~~

- Removed warnings for newer Python versions.
- Added barcoding functions to support cases with many distinct mutations at a site over time.

Version 1.3.4
-------------

Improvements
~~~~~~~~~~~~

- Added local-directory output for barcode and lineage metadata files.

Version 1.3.3
-------------

Bug fixes
~~~~~~~~~

- Fixed a rare barcoding bug caused by multiple mutation reversions at the same site.

Version 1.3.2
-------------

Improvements
~~~~~~~~~~~~

- Added extension specification to ``aggregate``.

Bug fixes
~~~~~~~~~

- Added bug fixes.

Version 1.3.1
-------------

Improvements
~~~~~~~~~~~~

- Added VCF input format support to ``demix``.
- Added custom lineage metadata option to ``variants``.

Version 1.3
-----------

Improvements
~~~~~~~~~~~~

- Added bootstrap confidence intervals.
- Added color control for plotting.

Bug fixes
~~~~~~~~~

- Added bug fixes.

.. _Versions 1.2.x (1.2.0 – 1.2.1):

Version 1.2.1
-------------

Improvements
~~~~~~~~~~~~

- Added custom barcode option to ``demix``.

Bug fixes
~~~~~~~~~

- Fixed aggregation.

Version 1.2
-----------

Improvements
~~~~~~~~~~~~

- Added helper functions and visualization tools to the CLI.

Bug fixes
~~~~~~~~~

- Added minor bug fixes.

.. _Versions 1.1.x (1.1.0):

Version 1.1
-----------

Improvements
~~~~~~~~~~~~

- Cleaned up user outputs.
- Cleaned up barcoding.

.. _Versions 1.0.x (1.0.0):

Version 1.0
-----------

Improvements
~~~~~~~~~~~~

- Added unit tests.