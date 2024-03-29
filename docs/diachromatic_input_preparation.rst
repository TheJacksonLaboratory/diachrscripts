.. _RST_Diachromatic_input_preparation:

##############################
Diachromatic input preparation
##############################

*************
bowtie2 index
*************

Diachromatic uses bowtie2 to map reads to a reference sequence, which requires a reference index.
Such an index can be created with bowtie2 or a pre-calculated index can e.g. be downloaded from the
`iGenomes website <https://support.illumina.com/sequencing/sequencing_software/igenome.html>`_:

.. code-block:: console

    $ wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/UCSC/hg38/Homo_sapiens_UCSC_hg38.tar.gz
    $ mkdir Homo_sapiens_UCSC_hg38
    $ tar -xf Homo_sapiens_UCSC_hg38.tar.gz -C Homo_sapiens_UCSC_hg38

In addition to the bowtie2 index, the downloaded directory contains other data and is about ``30 GB`` in size.
The index consists of several files:

.. code-block:: console

    $ ls Homo_sapiens_UCSC_hg38/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/ | cat
    genome.1.bt2
    genome.2.bt2
    genome.3.bt2
    genome.4.bt2
    genome.fa
    genome.rev.1.bt2
    genome.rev.2.bt2

Only the path together with the file prefix of these files are passed to Diachromatic ``-i <PATH>/genome``.


***********
Digest file
***********

Diachromatic reports interactions as pairs of restriction fragments (or digest pairs)
along with the number of supporting read pairs.
This requires the coordinates of all possible digests in the entire reference genome.
These are passed to Diachromatic in the form of a text file, which we refer to as digest file.
For a given restriction enzyme and a reference genome,
a corresponding `digest map can be created with the GOPHER software <https://diachromatic.readthedocs.io/en/latest/digest.html>`__.
For instance, here are the first few lines of a digest file
for the ``hg38`` genome digested with the restriction enzyme ``HindIII``:

.. code-block:: console

    Chromosome      Fragment_Start_Position Fragment_End_Position   Fragment_Number 5'_Restriction_Site     3'_Restriction_Site     Length  5'_GC_Content   3'_GC_Content   5'_Repeat_Content       3'_Repeat_Content       Selected        5'_Probes       3'_Probes
    chr1    1       16007   1       None    HindIII 16007   0.000   0.000   0.000   0.003   F       0       0
    chr1    16008   24571   2       HindIII HindIII 8564    0.018   0.018   0.000   0.015   F       0       0
    chr1    24572   27981   3       HindIII HindIII 3410    0.046   0.046   0.000   0.044   F       0       0
    chr1    27982   30429   4       HindIII HindIII 2448    0.035   0.035   0.047   0.043   F       0       0

Each line represents one digest (or restriction fragment).
The first three columns contain the digest coordinates.
In the ``Selected`` column, digests selected for enrichment are marked with a ``T`` and all others with an ``F``.
Diachromatic passes the information about enriched digests through to the reported interactions.
In Diachromatic interaction files, an ``E`` corresponds to a ``T`` and an ``N`` to an ``F``.

.. code-block:: console

    chr11    9641153    9642657   E   chr11   47259263   47272706   N   5:4:2:7

Selecting enriched digests
==========================

If you have created the capture bait design for your own experiment using GOPHER,
then all digests in the digest file for which baits were ordered will be marked with a ``T``.
For cases where the coordinates of baited digests are known, we provide a script to create a corresponding digest file.

.. code-block:: console

    $ python diachrscripts/additional_scripts/ed_selector.py
        --enriched-digests-file <ENRICHED_DIGEST_COORDINATES>.bed
        --diachromatic-digest-file <DIGEST_FILE_TEMPLATE>.txt
        --out-prefix <CUSTOM_DIGEST_FILE_PREFIX>

This script is passed a BED file containing the coordinates of the digests that have been selected for enrichment.
Such information can be found, for example,
in the supplementary material of the corresponding publication (see examples below).
In addition, a digest file for the appropriate reference genome and restriction enzyme must be passed.
You can export such a digest file from GOPHER or you can use an existing digest file.
The script will only use the file as a template and the enrichment information will be completely overwritten.
For convenience, we've added a template digest file to this repository that just needs to be unzipped:

.. code-block:: console

    $ gunzip -k diachrscripts/additional_files/template_digest_file_hg38_HindIII.txt.gz

Mifsud et al. 2015
------------------

Supplementary Table 4 of the work published by
`Mifsud et al. 2015 <https://pubmed.ncbi.nlm.nih.gov/25938943/>`__
contains the coordinates and sequences of the baits used.
Save this table in text format and extract the coordinates.

.. code-block:: console

    $ cat mifsud_supplementary_table_4.txt | \
        awk '{if($1 ~ /^>/){split($0,a," ");split(a[1],b,":");gsub(/>C/,"c",b[1]);split(b[2],c,"-");print b[1]"\t"c[1]-1"\t"c[2]}}' \
        > mifsud_baits_hg19.bed

Use
`UCSC's LiftOver tool <https://genome.ucsc.edu/cgi-bin/hgLiftOver>`_
to convert the coordinates from ``hg19`` to ``hg38``.
Deselect ``Allow multiple output regions``.
37,604 bait coordinates were successfully converted to ``hg38``.
The conversion failed for the coordinates of 4 baits because the corresponding regions in ``hg38`` are either deleted or
partially deleted.
Furthermore, the coordinates of two baits are mapped to chromosome ``chr22_KI270879v1_alt`` of ``hg38``.
These can be removed as follows:

.. code-block:: console

    $ grep -v 'chr22_KI270879v1_alt' lift_over_results.bed  > mifsud_baits_hg38.bed

Next, extract the coordinates of all digests in the genome from the digest file template and write them to a BED file:

.. code-block:: console

    $ tail -n+2 diachrscripts/additional_files/template_digest_file_hg38_HindIII.txt \
    | awk '{print $1"\t"$2"\t"$3}' > all_hg38_HindIII_digests.bed

Then use
`bedtools <https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html>`_
to extract all digests that contain at least one bait completely:

.. code-block:: console

    $ intersectBed -wa -u -F 1.00 -a all_hg38_HindIII_digests.bed -b mifsud_baits_hg38.bed \
    > mifsud_baited_digests_hg38.bed

This results in 22,076 baited digests.

Finally, use our script to create a digest file in which digests that Mifsud et al. have selected for enrichment
are marked with a ``T`` and all others with an ``F``.

.. code-block:: console

    $ python diachrscripts/additional_scripts/ed_selector.py \
        --enriched-digests-file mifsud_baited_digests_hg38.bed \
        --diachromatic-digest-file \
            diachrscripts/additional_files/template_digest_file_hg38_HindIII.txt \
        --out-prefix mifsud_hg38_HindIII

This will produce the file ``mifsud_hg38_HindIII_diachromatic_digest_file.txt`` that can be used as input for
Diachromatic.
All 22,076 digests in the digest file were marked with a ``T``.
We have added the file ``mifsud_baited_digests_hg38.bed`` to this repository so that the digest file can be recreated
if needed.

Javierre et al. 2016
--------------------

For the work published by
`Javierre et al. 2016 <https://pubmed.ncbi.nlm.nih.gov/27863249/>`__,
the ``hg19`` coordinates of the baited digests can be downloaded from
`OFS <https://osf.io/e594p/>`__.
First, download an archive that expands into a *design folder* that is intended as input for the interaction caller
`CHiCAGO <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4908757/>`_:

.. code-block:: console

    $ wget -O human_PCHiC_hg19_HindIII_design.tar.gz https://osf.io/e594p/download
    $ tar -xf human_PCHiC_hg19_HindIII_design.tar.gz

Along with other files, this folder contains
`CHiCAGO's bait map file <https://bioconductor.org/packages/devel/bioc/vignettes/Chicago/inst/doc/Chicago.html>`_
that consists of the following columns:
``chr``, ``start``, ``end``, ``fragmentID``, ``geneName``.

.. code-block:: console

    $ head -n 4 Human_hg19/Digest_Human_HindIII_baits_e75_ID.baitmap
        1	831895	848168	218	RP11-54O7.16;RP11-54O7.1
        1	848169	850618	219	RP11-54O7.2
        1	850619	874081	220	AL645608.1;RP11-54O7.3;SAMD11
        1	889424	903640	223	KLHL17;NOC2L;PLEKHN1

Coordinates are available for a total of 22,076 baited digests.
Next, convert the bait map file into BED format:

.. code-block:: console

    $ awk '{print "chr"$1"\t"$2"\t"$3}' Human_hg19/Digest_Human_HindIII_baits_e75_ID.baitmap \
    > javierre_baited_digests_hg19.bed

Then use
`UCSC's LiftOver tool <https://genome.ucsc.edu/cgi-bin/hgLiftOver>`_
to convert the coordinates from ``hg19`` to ``hg38`` and save the resulting file as
``javierre_baited_digests_hg38.bed``.
Deselect ``Allow multiple output regions``.
22,056 digest coordinates were successfully converted  to ``hg38``.
The conversion failed for 20 digests
because ``hg19`` coordinates in ``hg38``
are either split or partially deleted.
Finally, use our script to create a digest file in which digests that Javierre et al. have selected for enrichment are marked
with a ``T`` and all others with an ``F``.

.. code-block:: console

    $ python diachrscripts/additional_scripts/ed_selector.py \
        --enriched-digests-file javierre_baited_digests_hg38.bed \
        --diachromatic-digest-file \
            diachrscripts/additional_files/template_digest_file_hg38_HindIII.txt \
        --out-prefix javierre_hg38_HindIII

This will produce the file ``javierre_hg38_HindIII_diachromatic_digest_file.txt`` that can be used as input for
Diachromatic.
22,008 digests in the digest file were marked with a ``T``.
48 input digest could not be found in the digest file.
We examined these cases in more detail (``--verbose``) and concluded that these cases are due to the LiftOver step.
We have added the file ``javierre_baited_digests_hg38.bed`` to this repository so that the digest file can be recreated
if needed. For the files ``javierre_baited_digests_hg38.bed`` and ``mifsud_baited_digests_hg38.bed``,
22,008 baited digests overlap.

Schoenfelder et al. 2015
------------------------

Supplementary Table 1 of the work published by
`Schoenfelder et al. 2015 <https://pubmed.ncbi.nlm.nih.gov/25752748/>`__
contains the coordinates and sequences of the baits used.
Save this table in text format and extract the coordinates.

.. code-block:: console

    $ tail -n+2 SuppTable1.txt | awk '{print $1"\t"$2-1"\t"$3}' > schoenfelder_baits_mm9.bed

Use
`UCSC's LiftOver tool <https://genome.ucsc.edu/cgi-bin/hgLiftOver>`_
to convert the coordinates from ``mm9`` to ``mm10``.
Deselect ``Allow multiple output regions``.
39,019 bait coordinates were successfully converted to ``mm10``.
The conversion failed for the coordinates of 2 baits because the corresponding regions in ``mm10`` are  deleted.
We save the BED file with the converted coordinates as ``schoenfelder_baits_mm10.bed``.

Next, extract the coordinates of all digests in the genome from the digest file template and write them to a BED file:

.. code-block:: console

    $ tail -n+2 diachrscripts/additional_files/template_digest_file_mm10_HindIII.txt \
    | awk '{print $1"\t"$2"\t"$3}' > all_mm10_HindIII_digests.bed

Then use
`bedtools <https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html>`_
to extract all digests that contain at least one bait completely:

.. code-block:: console

    $ intersectBed -wa -u -F 1.00 -a all_mm10_HindIII_digests.bed -b schoenfelder_baits_mm10.bed \
    > schoenfelder_baited_digests_mm10.bed

This results in 22,224 baited digests.

Finally, use our script to create a digest file in which digests that Schoenfelder et al. have selected for enrichment
are marked with a ``T`` and all others with an ``F``.

.. code-block:: console

    $ python diachrscripts/additional_scripts/ed_selector.py \
        --enriched-digests-file schoenfelder_baited_digests_mm10.bed \
        --diachromatic-digest-file \
            diachrscripts/additional_files/template_digest_file_mm10_HindIII.txt \
        --out-prefix schoenfelder_mm10_HindIII

This will produce the file ``schoenfelder_mm10_HindIII_diachromatic_digest_file.txt`` that can be used as input for
Diachromatic.
