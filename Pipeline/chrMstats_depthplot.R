
margs<-commandArgs(trailingOnly=TRUE);

title <- margs[1];
xtitle <- margs[2];
ytitle <- margs[3];
outfile <- margs[4];

series<-read.table(file='stdin', header=FALSE, sep=" ");
#x11();

png(outfile, width=2000, height=1000, units="px");

plot(series$V1, series$V2,  main=title, xlab=xtitle, ylab=ytitle, xaxt="n", xlim=c(1,16569), ylim=c(-3250, 10200), yaxt="n");
#axis(1, c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,11000,12000,13000,14000, 15000, 16000, 16569));
axis(1, c(1,2000,4000,6000,8000,10000,12000,14000,16569), las=2);
axis(2, c(1,1000,2000,3000,4000,5000,6000,7000,8000, 9000, 10000), las=1);


#rect(xleft, ybottom, xright, ytop, density = NULL, angle = 45,
#     col = NA, border = NULL, lty = par("lty"), lwd = par("lwd"),
#     ...)

#text(x, y = NULL, labels = seq_along(x), adj = NULL,
#     pos = NULL, offset = 0.5, vfont = NULL,
#     cex = 1, col = NULL, font = NULL, ...

rect(1, -10 , 576, -300, col="yellow" );
text(((1+576)/2)-30, -3000, labels="Control_region", srt=90);
rect(57, -10 , 372, -300, col="red" );
text(((57+372)/2)-30, -1500, labels="HVS-II", srt=90);
rect(438, -10 , 574, -300, col="yellow" );
text(((438+574)/2)-30, -3000, labels="HVS-III", srt=90);
rect(577, -10 , 647, -300, col="red" );
text(((577+647)/2)-30, -1500, labels="tRNA-Phe", srt=90);
rect(648, -10 , 1601, -300, col="yellow" );
text(((648+1601)/2)-30, -3000, labels="12S_rRNA", srt=90);
rect(1602, -10 , 1670, -300, col="red" );
text(((1602+1670)/2)-30, -1500, labels="tRNA-Val", srt=90);
rect(1671, -10 , 3106, -300, col="yellow" );
text(((1671+3106)/2)-30, -3000, labels="16S_rRNA", srt=90);
rect(3108, -10 , 3229, -300, col="red" );
text(((3108+3229)/2)-30, -1500, labels="16S_rRNA.2", srt=90);
rect(3230, -10 , 3304, -300, col="yellow" );
text(((3230+3304)/2)-30, -3000, labels="tRNA-Leu2", srt=90);
rect(3307, -10 , 4262, -300, col="red" );
text(((3307+4262)/2)-30, -1500, labels="ND1", srt=90);
rect(4263, -10 , 4331, -300, col="yellow" );
text(((4263+4331)/2)-30, -3000, labels="tRNA-Ile", srt=90);
rect(4329, -10 , 4400, -300, col="red" );
text(((4329+4400)/2)-30, -1500, labels="tRNA-Gln", srt=90);
rect(4402, -10 , 4469, -300, col="yellow" );
text(((4402+4469)/2)-30, -3000, labels="tRNA-Met", srt=90);
rect(4470, -10 , 5511, -300, col="red" );
text(((4470+5511)/2)-30, -1500, labels="ND2", srt=90);
rect(5512, -10 , 5579, -300, col="yellow" );
text(((5512+5579)/2)-30, -3000, labels="tRNA-Trp", srt=90);
rect(5587, -10 , 5655, -300, col="red" );
text(((5587+5655)/2)-30, -1500, labels="tRNA-Ala", srt=90);
rect(5657, -10 , 5729, -300, col="yellow" );
text(((5657+5729)/2)-30, -3000, labels="tRNA-Asn", srt=90);
rect(5761, -10 , 5826, -300, col="red" );
text(((5761+5826)/2)-30, -1500, labels="tRNA-Cys", srt=90);
rect(5826, -10 , 5891, -300, col="yellow" );
text(((5826+5891)/2)-30, -3000, labels="tRNA-Tyr", srt=90);
rect(5904, -10 , 7445, -300, col="red" );
text(((5904+7445)/2)-30, -1500, labels="CO1", srt=90);
rect(7446, -10 , 7514, -300, col="yellow" );
text(((7446+7514)/2)-30, -3000, labels="tRNA-Ser2", srt=90);
rect(7518, -10 , 7585, -300, col="red" );
text(((7518+7585)/2)-30, -1500, labels="tRNA-Asp", srt=90);
rect(7586, -10 , 8269, -300, col="yellow" );
text(((7586+8269)/2)-30, -3000, labels="CO2", srt=90);
rect(8295, -10 , 8364, -300, col="red" );
text(((8295+8364)/2)-30, -1500, labels="tRNA-Lys", srt=90);
rect(8366, -10 , 8572, -300, col="yellow" );
text(((8366+8572)/2)-30, -3000, labels="ATP8", srt=90);
rect(8527, -10 , 9207, -300, col="red" );
text(((8527+9207)/2)-30, -1500, labels="ATP6", srt=90);
rect(9207, -10 , 9990, -300, col="yellow" );
text(((9207+9990)/2)-30, -3000, labels="CO3", srt=90);
rect(9991, -10 , 10058, -300, col="red" );
text(((9991+10058)/2)-30, -1500, labels="tRNA-Gly", srt=90);
rect(10059, -10 , 10404, -300, col="yellow" );
text(((10059+10404)/2)-30, -3000, labels="ND3", srt=90);
rect(10405, -10 , 10469, -300, col="red" );
text(((10405+10469)/2)-30, -1500, labels="tRNA-Arg", srt=90);
rect(10470, -10 , 10766, -300, col="yellow" );
text(((10470+10766)/2)-30, -3000, labels="ND4L", srt=90);
rect(10760, -10 , 12137, -300, col="red" );
text(((10760+12137)/2)-30, -1500, labels="ND4", srt=90);
rect(12138, -10 , 12206, -300, col="yellow" );
text(((12138+12206)/2)-30, -3000, labels="tRNA-His", srt=90);
rect(12207, -10 , 12265, -300, col="red" );
text(((12207+12265)/2)-30, -1500, labels="tRNA-Ser1", srt=90);
rect(12266, -10 , 12336, -300, col="yellow" );
text(((12266+12336)/2)-30, -3000, labels="tRNA-Leu1", srt=90);
rect(12337, -10 , 14148, -300, col="red" );
text(((12337+14148)/2)-30, -1500, labels="ND5", srt=90);
rect(14149, -10 , 14673, -300, col="yellow" );
text(((14149+14673)/2)-30, -3000, labels="ND6", srt=90);
rect(14674, -10 , 14742, -300, col="red" );
text(((14674+14742)/2)-30, -1500, labels="tRNA-Glu", srt=90);
rect(14747, -10 , 15887, -300, col="yellow" );
text(((14747+15887)/2)-30, -3000, labels="CYB", srt=90);
rect(15888, -10 , 15953, -300, col="red" );
text(((15888+15953)/2)-30, -1500, labels="tRNA-Thr", srt=90);
rect(15956, -10 , 16023, -300, col="yellow" );
text(((15956+16023)/2)-30, -3000, labels="tRNA-Pro", srt=90);
rect(16024, -10 , 16383, -300, col="red" );
text(((16024+16383)/2)-30, -1500, labels="HVS-I", srt=90);
rect(16024, -10 , 16569, -300, col="yellow" );
text(((16024+16569)/2)-30, -3000, labels="Control_region.2", srt=90);


title(title);
#, col=sample(colours(), 1));

dev.off();






