peak_fn=$1
outdir=$2

mkdir -p $outdir/topological_dist
python2 scripts/topological_dist/Peak_distribution_on_utr_cds.py $peak_fn mm10 100 >$outdir/topological_dist/dist.data
Rscript scripts/topological_dist/plot.R $outdir/topological_dist/dist.data $outdir/topological_dist/dist.pdf

mkdir -p $outdir/homer
findMotifsGenome.pl $peak_fn mm10 $outdir/homer \
	-rna -len 5,6,7 \
	-p 4 -size 100 -S 10
