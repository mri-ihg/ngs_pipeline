#
#
# calcOnOffTargetCoverageProfilePlot.R
#
# Plots two data series with the same X/Y scale provided from outside
# and saves it to PNG 
#
# WARNING: X POINTS MUST BE THE SAME
#
# usage
# 	Rscript calcOnOffTargetCoverageProfilePlot.R "PLOTTITLE" XMIN XMAX "XAXISLABEL" YMIN YMAX "YAXISLABEL" DATASERIESXY1.dat DATASERIESXY2.dat OUTFILE.png
# 
# 2016 Riccardo Berutti
#
#

# Get arguments
args<-commandArgs(TRUE);
    # Title
    title<-args[1];
    # X axis
    xmin<-as.numeric(args[2]);
    xmax<-as.numeric(args[3]);
    xtitle<-args[4];
    # Y axis
    ymin<-as.numeric(args[5]);
    ymax<-as.numeric(args[6]);
    ytitle<-args[7];
    # Data
    f1<-args[8];
    f2<-args[9];
    # Outfile
    outpng<-args[10];

# Set plot ticks and labels
myticks<-seq(xmin,xmax,  (xmax-xmin+1));
mylabels<-seq(xmin,xmax, (xmax-xmin+1)/10);

#Read series 
series1<-read.table(file=f1, header=FALSE, nrows=(xmax-xmin+1));
series2<-read.table(file=f2, header=FALSE, nrows=(xmax-xmin+1));

#Output PNG
png( outpng );

# Plot Series 1
plot(series1$V1, series1$V2 , main=title, xlab=xtitle, ylab=ytitle, xlim=c(xmin, xmax), col="red", pch=19);

# Plot Series 2 points
points(series2$V2, col="blue", pch=1);

# Set title
title(title);

# Draw legend
myleg <- c("On target", "Off target");
legend('topright', myleg , col=c('red', 'blue'), pch=c(19,1), bty='n', cex=.75);

# Close plotting device
dev.off();



