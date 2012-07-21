###############
5.6 closestBed
###############
Similar to **intersectBed, closestBed** searches for overlapping features in A and B. In the event that
no feature in B overlaps the current feature in A, **closestBed** will report the *closest* (that is, least
genomic distance from the start or end of A) feature in B. For example, one might want to find which
is the closest gene to a significant GWAS polymorphism. Note that **closestBed** will report an
overlapping feature as the closest---that is, it does not restrict to closest *non-overlapping* feature.

==========================================================================
5.6.1 Usage and option summary
==========================================================================
**Usage:**
::
  closestBed [OPTIONS] -a <BED/GFF/VCF> -b <BED/GFF/VCF>
  
  
===========================      ===============================================================================================================================================================================================================
Option                           Description
===========================      ===============================================================================================================================================================================================================
**-s**				             Force strandedness. That is, find the closest feature in B overlaps A on the same strand. *By default, this is disabled*.
**-d**					         In addition to the closest feature in B, report its distance to A as an extra column. The reported distance for overlapping features will be 0.
**-t**					         How ties for closest feature should be handled. This occurs when two features in B have exactly the same overlap with a feature in A. *By default, all such features in B are reported*.
                                 
								 Here are the other choices controlling how ties are handled:
								      				 
								 *all-*   Report all ties (default).
								 
								 *first-*   Report the first tie that occurred in the B file.
								 
								 *last-*   Report the last tie that occurred in the B file.
===========================      ===============================================================================================================================================================================================================




==========================================================================
5.6.2 Default behavior
==========================================================================
**closestBed** first searches for features in B that overlap a feature in A. If overlaps are found, the feature
in B that overlaps the highest fraction of A is reported. If no overlaps are found, **closestBed** looks for
the feature in B that is *closest* (that is, least genomic distance to the start or end of A) to A. For
example, in the figure below, feature B1 would be reported as the closest feature to A1.
::
  Chromosome  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  BED FILE A                             *************
  
  BED File B         ^^^^^^^^                            ^^^^^^
  
  Result                                                 ======
  

For example:
::
  cat A.bed
  chr1  100  200

  cat B.bed
  chr1  500  1000
  chr1  1300 2000

  closestBed -a A.bed -b B.bed
  chr1  100  200  chr1  500  1000



==========================================================================
5.6.3 (-s)Enforcing "strandedness" 
==========================================================================
This option behaves the same as the -s option for intersectBed while scanning for the closest
(overlapping or not) feature in B. See the discussion in the intersectBed section for details.



==========================================================================
5.6.4 (-t)Controlling how ties for "closest" are broken 
==========================================================================
When there are two or more features in B that overlap the *same fraction* of A, **closestBed** will, by
default, report both features in B. Imagine feature A is a SNP and file B contains genes. It can often
occur that two gene annotations (e.g. opposite strands) in B will overlap the SNP. As mentioned, the
default behavior is to report both such genes in B. However, the -t option allows one to optionally
choose the just first or last feature (in terms of where it occurred in the input file, not chromosome
position) that occurred in B.

For example (note the difference between -l 200 and -l 300):
::
  cat A.bed
  chr1  100  101  rs1234

  cat B.bed
  chr1  0  1000  geneA  100  +
  chr1  0  1000  geneB  100  -

  closestBed -a A.bed -b B.bed
  chr1  100  101  rs1234  chr1  0  1000  geneA  100  +
  chr1  100  101  rs1234  chr1  0  1000  geneB  100  -

  closestBed -a A.bed -b B.bed -t all
  chr1  100  101  rs1234  chr1  0  1000  geneA  100  +
  chr1  100  101  rs1234  chr1  0  1000  geneB  100  -

  closestBed -a A.bed -b B.bed -t first
  chr1  100  101  rs1234  chr1  0  1000  geneA  100  +

  closestBed -a A.bed -b B.bed -t last
  chr1  100  101  rs1234  chr1  0  1000  geneB  100  -






==========================================================================
5.6.5 (-d)Reporting the distance to the closest feature in base pairs 
==========================================================================
ClosestBed will optionally report the distance to the closest feature in the B file using the **-d** option.
When a feature in B overlaps a feature in A, a distance of 0 is reported.
::
  cat A.bed
  chr1  100  200
  chr1  500  600

  cat B.bed
  chr1  500  1000
  chr1  1300 2000

  closestBed -a A.bed -b B.bed -d
  chr1  100  200  chr1  500  1000  300
  chr1  500  600  chr1  500  1000  0
