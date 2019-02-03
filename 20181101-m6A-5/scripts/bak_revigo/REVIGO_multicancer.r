

# A plotting R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800


# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );
# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );


# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.



revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability","cluster");
revigo.data <- rbind(c("GO:0001128","RNA polymerase II transcription coactivator activity involved in preinitiation complex assembly", 0.006,-3.584, 4.241, 0.301,-1.3981,0.986,0.000,0),
                     c("GO:0003705","transcription factor activity, RNA polymerase II distal enhancer sequence-specific binding", 0.549,-1.439, 3.842, 1.982,-1.3981,0.986,0.000,0),
                     c("GO:0004869","cysteine-type endopeptidase inhibitor activity", 0.335,-4.185, 0.802, 1.771,-2.2963,0.963,0.000,-2),
                     c("GO:0035198","miRNA binding", 0.087, 8.327, 0.217, 1.204,-15.00,0.777,0.000,6),
                     c("GO:1902945","metalloendopeptidase activity involved in amyloid precursor protein catabolic process", 0.670, 0.295,-7.351, 2.068,-2.4628,0.812,0.000,-2),
                     c("GO:0043024","ribosomal small subunit binding", 0.075,-0.468, 5.837, 1.146,-1.8632,0.966,0.004,-2),
                     c("GO:0017124","SH3 domain binding", 0.699,-2.346, 4.631, 2.086,-2.2830,0.959,0.005,-2),
                     c("GO:0001540","beta-amyloid binding", 0.191,-1.709, 2.416, 1.531,-1.3981,0.981,0.005,-4),
                     c("GO:0005506","iron ion binding", 0.959,-0.551, 7.187, 2.223,-2.0376,0.969,0.006,-2),
                     c("GO:0051010","microtubule plus-end binding", 0.064,-4.138,-2.141, 1.079,-2.0769,0.947,0.027,-4),
                     c("GO:0070577","lysine-acetylated histone binding", 0.104,-1.951,-0.525, 1.279,-1.3981,0.975,0.028,0),
                     c("GO:0070530","K63-linked polyubiquitin binding", 0.116,-1.071,-0.967, 1.322,-1.4899,0.975,0.028,0),
                     c("GO:0030331","estrogen receptor binding", 0.237,-4.750,-0.236, 1.623,-1.4899,0.960,0.030,0),
                     c("GO:0045296","cadherin binding", 1.710,-2.671, 0.916, 2.473,-2.2720,0.946,0.037,-4),
                     c("GO:0031267","small GTPase binding", 1.774,-3.352,-3.783, 2.489,-1.7377,0.952,0.041,-2),
                     c("GO:1990259","histone-glutamine methyltransferase activity", 0.012,-0.682,-5.879, 0.477,-2.2963,0.891,0.098,-2),
                     c("GO:0098505","G-rich strand telomeric DNA binding", 0.046, 4.123, 4.333, 0.954,-3.0562,0.868,0.151,0),
                     c("GO:0005525","GTP binding", 2.225, 3.636, 6.651, 2.587,-1.7279,0.875,0.157,-2),
                     c("GO:0004577","N-acetylglucosaminyldiphosphodolichol N-acetylglucosaminyltransferase activity", 0.012,-0.006,-5.851, 0.477,-1.8632,0.885,0.185,-2),
                     c("GO:0004385","guanylate kinase activity", 0.075,-1.484,-5.980, 1.146,-1.6445,0.882,0.207,-2),
                     c("GO:0003723","RNA binding", 9.262, 6.346, 3.489, 3.205,-15.00,0.812,0.236,6),
                     c("GO:0001069","regulatory region RNA binding", 0.006, 6.896,-3.032, 0.301,-15.00,0.809,0.257,6),
                     c("GO:0033204","ribonuclease P RNA binding", 0.006, 6.050,-2.731, 0.301,-15.00,0.809,0.257,6),
                     c("GO:0030627","pre-mRNA 5'-splice site binding", 0.012, 5.572,-1.569, 0.477,-3.9087,0.777,0.268,6),
                     c("GO:0042835","BRE binding", 0.012, 7.168,-2.370, 0.477,-15.00,0.802,0.268,6),
                     c("GO:0061752","telomeric repeat-containing RNA binding", 0.012, 6.374,-2.209, 0.477,-15.00,0.802,0.268,6),
                     c("GO:0098808","mRNA cap binding", 0.012, 7.674,-2.460, 0.477,-1.7518,0.802,0.268,6),
                     c("GO:0034584","piRNA binding", 0.017, 6.209,-1.510, 0.602,-15.00,0.797,0.275,6),
                     c("GO:0070883","pre-miRNA binding", 0.023, 7.426,-1.523, 0.699,-15.00,0.794,0.280,6),
                     c("GO:0071208","histone pre-mRNA DCP binding", 0.023, 7.926,-1.534, 0.699,-15.00,0.794,0.280,6),
                     c("GO:0033592","RNA strand annealing activity", 0.029, 8.310,-1.139, 0.778,-2.9153,0.781,0.284,6),
                     c("GO:0034513","box H/ACA snoRNA binding", 0.029, 6.834,-1.172, 0.778,-1.8632,0.791,0.284,6),
                     c("GO:0002151","G-quadruplex RNA binding", 0.040, 5.709,-0.458, 0.903,-15.00,0.787,0.290,6),
                     c("GO:1990247","N6-methyladenosine-containing RNA binding", 0.040, 6.841,-0.626, 0.903,-15.00,0.787,0.290,6),
                     c("GO:0008312","7S RNA binding", 0.046, 6.276,-0.293, 0.954,-15.00,0.785,0.293,6),
                     c("GO:0035197","siRNA binding", 0.046, 8.377,-0.445, 0.954,-15.00,0.785,0.293,6),
                     c("GO:0070878","primary miRNA binding", 0.046, 7.635,-0.411, 0.954,-15.00,0.785,0.293,6),
                     c("GO:0035613","RNA stem-loop binding", 0.064, 6.016, 0.181, 1.079,-15.00,0.781,0.300,6),
                     c("GO:0004722","protein serine/threonine phosphatase activity", 0.393, 0.818,-7.361, 1.839,-2.0869,0.770,0.302,-2),
                     c("GO:0008094","DNA-dependent ATPase activity", 0.462, 1.486,-7.503, 1.908,-2.1238,0.764,0.307,-2),
                     c("GO:0000498","base pairing with RNA", 0.167, 6.356, 0.805, 1.462,-15.00,0.768,0.321,6),
                     c("GO:0000339","RNA cap binding", 0.098, 7.902, 0.413, 1.255,-15.00,0.775,0.324,6),
                     c("GO:0070034","telomerase RNA binding", 0.098, 7.367, 0.329, 1.255,-15.00,0.775,0.324,6),
                     c("GO:0005544","calcium-dependent phospholipid binding", 0.324, 1.160, 7.849, 1.756,-1.4899,0.969,0.326,-2),
                     c("GO:0017091","AU-rich element binding", 0.133, 7.768, 0.884, 1.380,-15.00,0.771,0.332,6),
                     c("GO:1904678","alpha-aminoacyl-tRNA binding", 0.167, 6.603, 1.210, 1.462,-15.00,0.768,0.338,6),
                     c("GO:0017020","myosin phosphatase regulator activity", 0.006,-4.336, 1.803, 0.301,-2.0376,0.953,0.346,-2),
                     c("GO:0031492","nucleosomal DNA binding", 0.168, 5.121, 4.778, 1.477,-1.6018,0.843,0.350,0),
                     c("GO:0016922","ligand-dependent nuclear receptor binding", 0.133,-3.495,-1.419, 1.380,-1.3981,0.961,0.351,-2),
                     c("GO:0019843","rRNA binding", 0.318, 6.225, 1.372, 1.748,-15.00,0.758,0.356,6),
                     c("GO:0034336","misfolded RNA binding", 0.167, 8.021, 1.110, 1.462,-15.00,0.768,0.356,6),
                     c("GO:0034583","21U-RNA binding", 0.167, 6.965, 0.814, 1.462,-15.00,0.768,0.356,6),
                     c("GO:1990605","GU repeat RNA binding", 0.167, 5.965, 0.770, 1.462,-15.00,0.768,0.356,6),
                     c("GO:0036002","pre-mRNA binding", 0.168, 7.413, 1.466, 1.477,-15.00,0.768,0.356,6),
                     c("GO:0030515","snoRNA binding", 0.173, 7.654, 1.841, 1.491,-15.00,0.767,0.357,6),
                     c("GO:0017069","snRNA binding", 0.208, 7.876, 1.607, 1.568,-15.00,0.765,0.362,6),
                     c("GO:0000049","tRNA binding", 0.283, 7.417, 1.970, 1.699,-15.00,0.771,0.372,6),
                     c("GO:0008199","ferric iron binding", 0.064,-1.502, 6.676, 1.079,-1.8244,0.971,0.385,-2),
                     c("GO:0017166","vinculin binding", 0.058,-1.783,-3.192, 1.041,-1.8632,0.949,0.386,-2),
                     c("GO:0003727","single-stranded RNA binding", 0.445, 6.751, 1.899, 1.892,-15.00,0.752,0.387,6),
                     c("GO:0031624","ubiquitin conjugating enzyme binding", 0.196,-2.664,-4.046, 1.544,-1.3553,0.960,0.391,-2),
                     c("GO:0003725","double-stranded RNA binding", 0.364, 6.850, 1.731, 1.806,-15.00,0.756,0.391,6),
                     c("GO:0043531","ADP binding", 0.196, 2.382, 7.651, 1.544,-1.6018,0.894,0.397,-2));
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
one.data$cluster <- as.numeric( as.character(one.data$cluster) );
#head(one.data);


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, log10_p_value, colour = cluster, size = plot_size), alpha = I(0.6) ) + scale_size_area();
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( -5, 8) );
p1 <- p1 + geom_point( aes(plot_X, log10_p_value, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.1) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(3, 12)); #scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.001, ]; 
ex$log10_p_value<-ex$log10_p_value-0.3
p1 <- p1 + geom_text( data = ex, aes(plot_X, log10_p_value, label = description), colour = I(alpha("black", 0.85)), size =  5);
p1 <- p1 + labs (y = "Corrected -log10(p-value)", x = "Semantic distance",size="Scaled\nterm size", 
                 col="Term semantic distance" );
p1 <- p1 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$log10_p_value) - min(one.data$log10_p_value);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(0,max(one.data$log10_p_value)+one.y_range/10);
p1<-p1 + theme_classic() +guides(col=FALSE)+#theme(legend.position="bottom") 
  theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.title = element_text(size = 15),legend.title=element_text(size=15))
# --------------------------------------------------------------------------
# Output the plot to screen

#p1;

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

ggsave("revigo-plot_BRCA_shared.pdf",width=5, height=4, units='in')