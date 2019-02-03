

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


args = commandArgs(TRUE)

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data = read.csv(args[1], header=T)

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
num.of.cluster=3         ###change###
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
ex <- one.data [ one.data$dispensability < 0.01, ];   #choose the ones that are really independent (smaller dispensability = more independent)
ex$log10_p_value<-ex$log10_p_value-0.3        #plot_X and log10_p_value will be used as coordinates for labels
p1 <- p1 + geom_text( data = ex, aes(plot_X, log10_p_value, label = description), colour = I(alpha("black", 0.85)), size =  5);     #add the labels onto the plot
p1 <- p1 + labs (y = "Corrected -log10(p-value)", x = "Semantic distance",size="Scaled\nterm size", 
                 col="Term semantic distance" );
p1 <- p1 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$log10_p_value) - min(one.data$log10_p_value);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);  #we add some margin (one.x_range/10) on each side 
p1 <- p1 + ylim(min(one.data$log10_p_value)-one.y_range/5,max(one.data$log10_p_value)+one.y_range/5);      #we add some margin (one.y_range/5) on the top
p1<-p1 + theme_classic() +guides(col=FALSE)+#theme(legend.position="bottom") 
  theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.title = element_text(size = 15),legend.title=element_text(size=15))
# --------------------------------------------------------------------------
# Output the plot to screen

#p1;

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

ggsave(args[2],width=7, height=7, units='in')
