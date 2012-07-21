###############
5.20 bedToIgv
###############
**bedToIgv** creates an IGV (http://www.broadinstitute.org/igv/) batch script (see: http://
www.broadinstitute.org/igv/batch for details) such that a ¡°snapshot¡± will be taken at each features in a
feature file. This is useful as an efficient means for quickly collecting images of primary data at several
loci for subsequent screening, etc.

**NOTE: One must use IGV version 1.5 or higher.**



==========================================================================
5.20.1 Usage and option summary
==========================================================================
Usage:
::
  bedToIgv [OPTIONS] -i <BED/GFF/VCF> > <igv.batch>
  
  
===========================      ===============================================================================================================================================================================================================
 Option                           Description
===========================      ===============================================================================================================================================================================================================
**-path**				         The full path to which the IGV snapshots should be written. *Default: ./*		 
**-sess**					     The full path to an existing IGV session file to be loaded prior to taking snapshots. *Default is for no session to be loaded and the assumption is that you already have IGV open and loaded with your relevant data prior to running the batch script*.
**-sort**                        The type of BAM sorting you would like to apply to each image. **Valid sorting options**: *base, position, strand, quality, sample, and readGroup Default is to apply no sorting at all*.
**-clps**                        Collapse the aligned reads prior to taking a snapshot. *Default is to not collapse*.
**-name**                        Use the "name" field (column 4) for each image's filename. *Default is to use the "chr:start-pos.ext"*.
**-slop**                        Number of flanking base pairs on the left & right of the image.
**-img**                         The type of image to be created. **Valid options**: *png, eps, svg Default is png*.
===========================      ===============================================================================================================================================================================================================





==========================================================================
5.20.2 Default behavior
==========================================================================
Figure:
::
  bedToIgv -i data/rmsk.hg18.chr21.bed | head -9
  snapshotDirectory ./
  goto chr21:9719768-9721892
  snapshot chr21:9719768-9721892.png
  goto chr21:9721905-9725582
  snapshot chr21:9721905-9725582.png
  goto chr21:9725582-9725977
  snapshot chr21:9725582-9725977.png
  goto chr21:9726021-9729309
  snapshot chr21:9726021-9729309.png

  
  

==========================================================================
5.20.3 Using a bedToIgv batch script within IGV.
==========================================================================
Once an IGV batch script has been created with **bedToIgv**, it is simply a matter of running it from
within IGV.

For example, first create the batch script:
::
  bedToIgv -i data/rmsk.hg18.chr21.bed > rmsk.igv.batch
  
Then, open and launch the batch script from within IGV. This will immediately cause IGV to begin
taking snapshots of your requested regions.

