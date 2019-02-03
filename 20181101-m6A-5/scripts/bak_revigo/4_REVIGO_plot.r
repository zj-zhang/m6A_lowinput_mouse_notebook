

# A plotting R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800


# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );
library( RColorBrewer );
library( ggridges )
# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );


# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.



revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0071336","regulation of hair follicle cell proliferation", 0.012,-5.635, 0.847, 0.477,-1.2909,0.504,0.000),
                     c("GO:1901799","negative regulation of proteasomal protein catabolic process", 0.190, 5.428, 0.693, 1.531,-1.2909,0.870,0.026),
                     c("GO:0006656","phosphatidylcholine biosynthetic process", 0.208, 0.856,-5.112, 1.568,-1.2909,0.917,0.034),
                     c("GO:0060784","regulation of cell proliferation involved in tissue homeostasis", 0.006,-3.970, 0.322, 0.301,-1.2909,0.517,0.234),
                     c("GO:2000101","regulation of mammary stem cell proliferation", 0.006,-1.391, 2.532, 0.301,-1.2909,0.510,0.234),
                     c("GO:0060723","regulation of cell proliferation involved in embryonic placenta development", 0.012,-1.571, 4.408, 0.477,-1.2909,0.481,0.326),
                     c("GO:0070663","regulation of leukocyte proliferation", 1.171,-4.592, 3.200, 2.310,-1.2909,0.363,0.326),
                     c("GO:1901382","regulation of chorionic trophoblast cell proliferation", 0.023,-3.036, 1.311, 0.699,-1.2909,0.487,0.344),
                     c("GO:0071863","regulation of cell proliferation in bone marrow", 0.035,-6.529, 2.398, 0.845,-1.2909,0.476,0.355),
                     c("GO:2000291","regulation of myoblast proliferation", 0.058,-5.543, 1.807, 1.041,-1.2909,0.462,0.370),
                     c("GO:0070344","regulation of fat cell proliferation", 0.063,-4.213, 1.532, 1.079,-1.2909,0.460,0.373),
                     c("GO:1901722","regulation of cell proliferation involved in kidney development", 0.081,-2.510, 3.011, 1.176,-1.2909,0.443,0.381),
                     c("GO:0033688","regulation of osteoblast proliferation", 0.133,-6.306, 3.427, 1.380,-1.2909,0.438,0.398),
                     c("GO:0060251","regulation of glial cell proliferation", 0.138,-2.668, 3.899, 1.398,-1.2909,0.420,0.400),
                     c("GO:0010464","regulation of mesenchymal cell proliferation", 0.202,-3.871, 2.455, 1.556,-1.2909,0.425,0.414),
                     c("GO:0060043","regulation of cardiac muscle cell proliferation", 0.202,-3.066, 4.964, 1.556,-1.2909,0.335,0.414),
                     c("GO:2000495","regulation of cell proliferation involved in compound eye morphogenesis", 0.331,-3.872, 4.921, 1.763,-1.2909,0.336,0.434),
                     c("GO:0072091","regulation of stem cell proliferation", 0.340,-4.780, 2.562, 1.778,-1.2909,0.408,0.435),
                     c("GO:2000177","regulation of neural precursor cell proliferation", 0.427,-5.647, 4.478, 1.875,-1.2909,0.400,0.445),
                     c("GO:0048145","regulation of fibroblast proliferation", 0.490,-5.426, 3.178, 1.934,-1.2909,0.395,0.451),
                     c("GO:1901645","regulation of synoviocyte proliferation", 0.544,-5.150, 4.626, 1.978,-1.2909,0.391,0.456),
                     c("GO:1904002","regulation of sebum secreting cell proliferation", 0.544,-5.588, 3.902, 1.978,-1.2909,0.357,0.456),
                     c("GO:0035206","regulation of hemocyte proliferation", 1.111,-3.624, 3.762, 2.286,-1.2909,0.356,0.493),
                     c("GO:0050678","regulation of epithelial cell proliferation", 1.714,-4.790, 3.762, 2.474,-1.2909,0.348,0.518),
                     c("GO:1904073","regulation of trophectodermal cell proliferation", 0.006,-1.884, 5.420, 0.301,-1.2909,0.469,0.523),
                     c("GO:1903860","negative regulation of dendrite extension", 0.012, 2.005, 7.146, 0.477,-1.0073,0.674,0.528),
                     c("GO:0042127","regulation of cell proliferation", 8.898,-4.470, 3.708, 3.188,-1.2909,0.333,0.596),
                     c("GO:2000136","regulation of cell proliferation involved in heart morphogenesis", 0.098,-3.868, 5.464, 1.255,-1.2909,0.390,0.660),
                     c("GO:0003420","regulation of growth plate cartilage chondrocyte proliferation", 0.012,-3.052, 5.899, 0.477,-1.2909,0.425,0.696));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- -as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
one.data$log10_p_value <- as.numeric( as.character(one.data$log10_p_value) );
one.data$log10_p_value <- -one.data$log10_p_value
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );

#k-means clustering of data
num.of.cluster=2         ###change###
#cluster=kmeans(one.data[,c("plot_X","plot_Y")],2)$cluster
cluster=as.character(kmeans(one.data[,"plot_X"],num.of.cluster)$cluster)
one.data$cluster <- cluster;
#head(one.data);

# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, log10_p_value, colour = cluster, size = plot_size), alpha = I(0.6) ) + scale_size_area();
#p1 <- p1 + scale_colour_gradientn( colours = c("purple","orange"), limits = c( 1, 2) );
p1 <- p1 + scale_color_manual(values = brewer.pal(num.of.cluster, "Dark2"));      #color the dots by cluster
p1 <- p1 + geom_point( aes(plot_X, log10_p_value, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.1) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(3, 12)); #scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.02, ];   #choose the ones that are really independent (smaller dispensability = more independent)
ex$log10_p_value<-ex$log10_p_value-0.3        #plot_X and log10_p_value will be used as coordinates for labels
p1 <- p1 + geom_text( data = ex, aes(plot_X, log10_p_value, label = description), colour = I(alpha("black", 0.85)), size =  5);     #add the labels onto the plot
p1 <- p1 + labs (y = "Corrected -log10(p-value)", x = "Semantic distance",size="Scaled\nterm size", 
                 col="Term semantic distance" );
p1 <- p1 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$log10_p_value) - min(one.data$log10_p_value);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);  #we add some margin (one.x_range/10) on each side 
p1 <- p1 + ylim(0,max(one.data$log10_p_value)+one.y_range/5);      #we add some margin (one.y_range/5) on the top
p1<-p1 + theme_classic() +guides(col=FALSE)+#theme(legend.position="bottom") 
  theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.title = element_text(size = 15),legend.title=element_text(size=15))
# --------------------------------------------------------------------------
# Output the plot to screen

#p1;

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

#setwd("~/Desktop")
ggsave("Revigo_input_all_region_logit_JC_A3SS_pvalue.txt_BP.pdf",width=5, height=4, units='in')
