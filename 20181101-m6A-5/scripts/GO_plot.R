argv = commandArgs(TRUE)
sig_fn=argv[1]
bg_fn=argv[2]
outdir=argv[3]
out_name=argv[4]

source('scripts/GO_enrichr_plot.R')
#1.1for customized search
background_gene_table = NULL
for(fn in strsplit(bg_fn, ",")[[1]]){
	background_gene_table = rbind.data.frame(background_gene_table, read.table(fn), stringsAsFactors=F)
}
observed_gene_table = NULL
for(fn in strsplit(sig_fn, ",")[[1]]){
	observed_gene_table <- rbind.data.frame(observed_gene_table, read.table(fn), stringsAsFactors=F)
}

# Remove duplicates and put it in correct type for enrichR
background_gene_list <- unique(background_gene_table[,1])
observed_gene_list <- unique(observed_gene_table[,1])

GO_table <- get_enrichr_GO(background_gene_list, observed_gene_list)

save_GO_result(GO_table,output_dir=outdir, output_prefix=out_name)
# You can adjust the aspect ratio and bar width according to the size of the result 
pdf(file=paste0(outdir,'/',out_name,'.pdf',sep=''), width=10, height=7)
generate_go_plot_dev(GO_table, cutoff=0.05,num_top = 5,term_size=10, bar_width = 0.6, aspect_ratio = 1.2)
dev.off()

