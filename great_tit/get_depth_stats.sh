#!/bin/bash

set -e
set -u
set -o pipefail

if [ -f depth_summary.txt ] ; then
    rm depth_summary.txt
fi

if [ -f depth_hist.txt ] ; then
    rm depth_hist.txt
fi

if [ ! -d plots ] ; then
    mkdir plots
fi


METRICS_DIR=$1

ls $1/*.wgsmetrics_file.txt | while read file; do basename $file *.wgsmetrics_file.txt | cut -f1 -d '.' >> samples.txt; done

echo "Sample" > header.txt
cat header.txt samples.txt > sample_names.txt
rm samples.txt
rm header.txt

metrics_file=($( ls $METRICS_DIR/*.wgsmetrics_file.txt | cut -f1 ))
grep ^GENOME_TERRITORY ${metrics_file[1]} > depth_metrics.txt # get the header information
grep -m 1 ^coverage ${metrics_file[1]} > histogram_header.txt

for file in ${metrics_file[@]}
do
    grep -A 1 ^GENOME_TERRITORY $file | tail -n 1 >> depth_metrics.txt 
done

echo "Sample" > sample_rows.txt
grep -A 250 ^coverage $METRICS_DIR/*.wgsmetrics_file.txt | cut -f1 -d'.' | cut -f7 -d'/' >> sample_rows.txt 
cat $METRICS_DIR/*.wgsmetrics_file.txt | grep ^coverage -A 251 | grep -v ^coverage > depth_hist.txt

cat histogram_header.txt depth_hist.txt > depth_histograms.txt


paste sample_names.txt depth_metrics.txt > depth_summary.txt
paste sample_rows.txt depth_histograms.txt | grep -v -e "--" > depth_hist.txt

rm histogram_header.txt
rm depth_metrics.txt
rm depth_histograms.txt
rm sample_names.txt
rm sample_rows.txt


echo '---' > depth_summary_report.md
echo "title: 'Plots and summaries of coverage for each sample'"  >> depth_summary_report.md
echo 'fontsize: 8pt' >>  depth_summary_report.md
echo 'geometry: margin=0.4in' >>  depth_summary_report.md
echo '---' >> depth_summary_report.md


echo 'This report plots the coverage metrics gathered by the CollectWgsMetrics tool from Picard.
The definition for the reported metrics can be found
[here](https://broadinstitute.github.io/picard/picard-metric-definitions.html#CollectWgsMetrics.WgsMetrics). ' >> depth_summary_report.md
echo '  '  >> depth_summary_report.md
echo '  '  >> depth_summary_report.md
# Plot the depth means and histograms
echo "Table 1: Coverage metrics produced by Picard's CollectWgsMetrics  " >> depth_summary_report.md
echo '  ' >> depth_summary_report.md

Rscript plot_depth_summaries.R depth_summary.txt depth_hist.txt 

Rscript plot_z_to_aut.R z_to_aut_cov.txt Z_to_aut_coverage_bgi.pdf


mv *.pdf plots/
mv depth_hist.txt plots/
mv depth_summary.txt plots/
mv depth_summary_report.md plots/


echo '  ' >> depth_summary_report.md
echo '  ' >> depth_summary_report.md
echo '![Mean depth of coverage (+/- standard deviation) across whole 
genome calcualated with the CollectWgsMetrics. Only
reads with a minimum mapping quality of 20 and sites with a minimum base quality of 
20 were considered.](mean_depth.pdf)  ' >>  plots/depth_summary_report.md
echo ' ' >>  plots/depth_summary_report.md


echo '![Histrogram of genome coverage calculated using CollectWgsMetrics .
 Only reads with a minimum mapping quality of 20 and sites with a minimum base 
 quality of 20 were considered. Histograms show a maximum values of 140 for ease of
 visualisation.](depth_histograms.pdf)  ' >>  plots/depth_summary_report.md 
echo ' ' >>  plots/depth_summary_report.md


echo '![The fraction of bases that attained at least 5x and up to 100x 
sequence coverage in post-filtering bases. Only reads with a minimum mapping 
quality of 20 and sites with a minimum base quality of 20 were considered.](
fraction_covered.pdf)  ' >>  plots/depth_summary_report.md 
echo ' ' >>  plots/depth_summary_report.md
echo ' ' >>  plots/depth_summary_report.md 

echo '![The ratio of coverages on the Z chromosome to autosomes. The horizontal dashed
line shows the expected 0.5 for female samples](Z_to_aut_coverage_bgi.pdf)'  >>  plots/depth_summary_report.md

tar czf plots.tar.gz plots
rm -r plots/


