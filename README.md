# scCNV_heatmap v1.1.1 - plotting a ploidy sorted heatmap from scCNV data (10x Genomics)
<pre>
A bar with ploidies will be plotted on the right side of the heatmap.

Note: Current version uses the GRCh37 (hg19) genome. If a different species or genome build is desired, 
the "chromosome_sizes" list can be adjusted in the skript.

Dependencies: csv, argparse, cairosvg

Usage:
-o  Output filename
-c  node_unmerged_cnv_calls.bed file from the 10x pipeline output
-p  per_cell_summary_metrics.csv file from the 10x pipeline output

Optional:
-e  Set to "no" to include noisy cells marked by the 10x pipeline. Default is "yes"

<img src="https://raw.githubusercontent.com/StefanKurtenbach/scCNV_heatmap/master/heatmap.png" width = "300">

</pre>
