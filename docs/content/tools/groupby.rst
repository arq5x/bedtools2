.. _groupby:

###############
*groupby*
###############
``bedtools groupby`` is a useful tool that mimics the "group by" clause in database 
systems. Given a file or stream that is sorted by the appropriate "grouping 
columns" (``-g``), ``groupby`` will compute summary statistics on another 
column (``-c``) in the file or stream. This will work with output from all 
BEDTools as well as any other tab-delimited file or stream.  As such, this is a 
generally useful tool for all command-line analyses, not just genomics related 
research.

.. note::
    When using ``bedtools groupby``, the input data must be ordered by the same
    columns as specified with the ``-grp`` argument, which establish which 
    columns should be used to define a group of similar data. For example, 
    if ``-grp`` is 1,2,3, the data should be pre-grouped accordingly. When 
    ``bedtools groupby`` detects changes in the group columns it then summarizes 
    all lines with that group.  For example, ``sort -k1,1 -k2,2 -k3,3 data.txt
    | bedtools groupby -g 1,2,3 -c 4 -o mean``.


==========================================================================
Usage and option summary
==========================================================================
Usage
::

    bedtools groupby [OPTIONS] -i <input> -g <group columns> -c <op. column> -o <operation>
  
or:
::

    groupBy [OPTIONS] -i <input> -g <group columns> -c <op. column> -o <operation>
  
===========================      ===============================================================================================================================================================================================================
Option                           Description
===========================      ===============================================================================================================================================================================================================
**-i**                           The input file that should be grouped and summarized. *Use "stdin" when using piped input*. 
                                 **Note: if -i is omitted, input is assumed to come from standard input (stdin)**
**-g (-grp)**					 Specifies which column(s) (1-based) should be used to group the input. Columns may be comma-separated with each column must be explicitly listed. Or, ranges (e.g. 1-4) are also allowed. *Default: 1,2,3*
**-c (-opCol)**                  Specify the column (1-based) that should be summarized. *Required*.
**-o (-op)**                     Specify the operation that should be applied to **opCol**.

                                 | Valid operations: 
                                 
                                 | **sum** - *numeric only*
                                 | **count** - *numeric or text*
                                 | **count_distinct** - *numeric or text*
                                 | **min** - *numeric only*
                                 | **max** - *numeric only*
                                 | **mean** - *numeric only*
                                 | **median** - *numeric only*
                                 | **mode** - *numeric or text*
                                 | **antimode** - *numeric or text*
                                 | **stdev** - *numeric only*
                                 | **sstdev** - *numeric only*
                                 | **collapse** (i.e., print a comma separated list) - *numeric or text*
                                 | **distinct** (i.e., print a comma separated list) - *numeric or text*
                                 | **concat** (i.e., print a comma separated list) - *numeric or text*
                                 | **freqasc** - *print a comma separated list of values observed and the number of times they were observed.* 
                                 |     Reported in ascending order of frequency*
                                 | **freqdesc** - *print a comma separated list of values observed and the number of times they were observed.*
                                 |     Reported in descending order of frequency*
                                 | **first** - *numeric or text*
                                 | **last** - *numeric or text*   
                                 | *Default: sum*
===========================      ===============================================================================================================================================================================================================





==========================================================================
Default behavior.
==========================================================================
Let's imagine we have three incredibly interesting genetic variants that we are 
studying and we are interested in what annotated repeats these variants overlap.

.. code-block:: bash

  $ cat variants.bed
  chr21  9719758 9729320 variant1
  chr21  9729310 9757478 variant2
  chr21  9795588 9796685 variant3

  $ bedtools intersect -a variants.bed -b repeats.bed -wa -wb > variantsToRepeats.bed
  $ cat variantsToRepeats.bed
  chr21  9719758 9729320 variant1   chr21  9719768 9721892 ALR/Alpha   1004  +
  chr21  9719758 9729320 variant1   chr21  9721905 9725582 ALR/Alpha   1010  +
  chr21  9719758 9729320 variant1   chr21  9725582 9725977 L1PA3       3288  +
  chr21  9719758 9729320 variant1   chr21  9726021 9729309 ALR/Alpha   1051  +
  chr21  9729310 9757478 variant2   chr21  9729320 9729809 L1PA3       3897  -
  chr21  9729310 9757478 variant2   chr21  9729809 9730866 L1P1        8367  +
  chr21  9729310 9757478 variant2   chr21  9730866 9734026 ALR/Alpha   1036  -
  chr21  9729310 9757478 variant2   chr21  9734037 9757471 ALR/Alpha   1182  -
  chr21  9795588 9796685 variant3   chr21  9795589 9795713 (GAATG)n    308   +
  chr21  9795588 9796685 variant3   chr21  9795736 9795894 (GAATG)n    683   +
  chr21  9795588 9796685 variant3   chr21  9795911 9796007 (GAATG)n    345   +
  chr21  9795588 9796685 variant3   chr21  9796028 9796187 (GAATG)n    756   +
  chr21  9795588 9796685 variant3   chr21  9796202 9796615 (GAATG)n    891   +
  chr21  9795588 9796685 variant3   chr21  9796637 9796824 (GAATG)n    621   +

  
We can see that variant1 overlaps with 3 repeats, variant2 with 4 and variant3 
with 6. We can use bedtools groupby to summarize the hits for each variant in 
several useful ways. The default behavior is to compute the *sum* of the opCol.

.. code-block:: bash

  $ bedtools groupby -i variantsToRepeats.bed -g 1,2,3 -c 9
  chr21 9719758 9729320 6353
  chr21 9729310 9757478 14482
  chr21 9795588 9796685 3604



==========================================================================
Computing the min and max.
==========================================================================
Now let's find the *min* and *max* repeat score for each variant. We do this 
by "grouping" on the variant coordinate columns (i.e. cols. 1,2 and 3) and 
ask for the min and max of the repeat score column (i.e. col. 9).

.. code-block:: bash

  $ bedtools groupby -i variantsToRepeats.bed -g 1,2,3 -c 9 -o min
  chr21 9719758 9729320 1004
  chr21 9729310 9757478 1036
  chr21 9795588 9796685 308
  
We can also group on just the *name* column with similar effect.

.. code-block:: bash

  $ bedtools groupby -i variantsToRepeats.bed -g 4 -c 9 -o min
  variant1 1004
  variant2 1036
  variant3 308
  
What about the *max* score? Let's keep the coordinates and the name of the 
variants so that we stay in BED format.

.. code-block:: bash

  $ bedtools groupby -i variantsToRepeats.bed -grp 1-4 -c 9 -o max
  chr21 9719758 9729320 variant1 3288
  chr21 9729310 9757478 variant2 8367
  chr21 9795588 9796685 variant3 891



==========================================================================
Computing the mean and median.
==========================================================================
Now let's find the *mean* and *median* repeat score for each variant.

.. code-block:: bash

  $ cat variantsToRepeats.bed | bedtools groupby -g 1-4 -c 9 -o mean
  chr21 9719758 9729320 variant1 1588.25
  chr21 9729310 9757478 variant2 3620.5
  chr21 9795588 9796685 variant3 600.6667

  $ bedtools groupby -i variantsToRepeats.bed -g 1-4 -c 9 -o median
  chr21 9719758 9729320 variant1 1030.5
  chr21 9729310 9757478 variant2 2539.5
  chr21 9795588 9796685 variant3 652


==========================================================================
Computing the mode and "antimode".
==========================================================================
Now let's find the *mode* and *antimode* (i.e., the least frequent) repeat 
score for each variant (in this case they are identical).

.. code-block:: bash

  $ bedtools groupby -i variantsToRepeats.bed -g 1-4 -c 9 -o mode
  chr21 9719758 9729320 variant1 1004
  chr21 9729310 9757478 variant2 1036
  chr21 9795588 9796685 variant3 308

  $ bedtools groupby -i variantsToRepeats.bed -g 1-4 -c 9 -o antimode
  chr21 9719758 9729320 variant1 1004
  chr21 9729310 9757478 variant2 1036
  chr21 9795588 9796685 variant3 308

  
  
==========================================================================
Computing the count of lines for a given group.
==========================================================================
Figure:

.. code-block:: bash

  $ bedtools groupby -i variantsToRepeats.bed -g 1-4 -c 9 -o count
  chr21 9719758 9729320 variant1 4
  chr21 9729310 9757478 variant2 4
  chr21 9795588 9796685 variant3 6


  
  
============================================================================
Collapsing: listing all of the values in the opCol for a given group.
============================================================================
Now for something different. What if we wanted all of the names of the repeats 
listed on the same line as the variants? Use the collapse option. This 
"denormalizes" things. Now you have a list of all the repeats on a single line.

.. code-block:: bash

  $ bedtools groupby -i variantsToRepeats.bed -grp 1-4 -c 9 -o collapse
  chr21 9719758 9729320 variant1 ALR/Alpha,ALR/Alpha,L1PA3,ALR/Alpha,
  chr21 9729310 9757478 variant2 L1PA3,L1P1,ALR/Alpha,ALR/Alpha,
  chr21 9795588 9796685 variant3 (GAATG)n,(GAATG)n,(GAATG)n,(GAATG)n,(GAATG)n,(GAATG)n,



==========================================================================
Computing frequencies: freqasc and freqdesc.
==========================================================================
What if we want to report each distinct value along with its number of 
occurrence (much like uniq -c)?  The ``freqasc`` and freqdesc`` operations
handle this.

.. code-block:: bash

  $ cat variantsToRepeats.bed | bedtools groupby -g 1 -c 8 -o freqdesc
  chr21 (GAATG)n:6,ALR/Alpha:5,L1PA3:2,L1P1:1,
  
  $ cat variantsToRepeats.bed | bedtools groupby -g 1 -c 8 -o freqasc
  chr21 L1P1:1,L1PA3:2,ALR/Alpha:5,(GAATG)n:6,
  
