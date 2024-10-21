#######################################################################################
#######################################################################################
#######################################################################################

#Author: WILLIAN T.A.F. SILVA
#E-mail: willian.silva@evobiolab.com
#Website: http://www.evobiolab.com/
#Copyright: CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)

#######################################################################################
#######################################################################################
#######################################################################################

#######################################################################################
#Plot style.

linethickness<-3
pointsize<-1.7
fontsize<-pointsize

#######################################################################################
#ggplot2 theme.

myggplottheme<-theme(title=element_text(size=10,face="bold"),
                     axis.title=element_text(size=10,face="bold"),
                     axis.text=element_text(size=10),
                     axis.text.x=element_text(angle=60,size=8,vjust=0.5),
                     legend.position="none",
                     legend.title=element_text(size=10,face="bold"),
                     legend.text=element_text(size=10),
                     legend.key=element_blank(),
                     panel.grid=element_line(colour="gray90"),
                     panel.grid.major.x=element_blank(),
                     panel.grid.minor.x=element_blank(),
                     panel.background=element_rect(fill="white",colour="black"),
                     panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                     strip.background=element_rect(colour="black",
                                                   fill="white"))

#######################################################################################
#R base plots.

blankplot<-function(){
  plot(1,type="n",axes=F,xlab="",ylab="")
}

plot_pointtypes<-function(){
	plot(1:25,rep(0,25),pch=1:25,col="black",
		main="Point types",xlab="pch",ylab=NA,yaxt="n",
		cex=pointsize,cex.main=pointsize,cex.lab=pointsize,cex.axis=pointsize,
		col.axis="black") 
}

plot_linetypes<-function(){
	plot(1:6,1:6,pch=NA,col=NULL,
		main="Line types",xlab=NA,ylab="lty",xaxt="n",
		cex=pointsize,cex.main=pointsize,cex.lab=pointsize,cex.axis=pointsize,
		col.axis="black")
	#axis(2,col="blue")
	abline(h=1:6,lty=1:6,col="black",lwd=linethickness)
}

plot_rainbowcolors<-function(){
  plot(1:10,rep(10,10),pch=15,col=rainbow(10),
       main="Rainbow colors",xlab="rainbow(n)[y]",ylab="rainbow(n)",ylim=c(1,10),#yaxt="n",
       cex=pointsize,cex.main=pointsize,cex.lab=pointsize,cex.axis=pointsize,
       col.axis="black")
  for(i in 1:9){
    points(1:i,rep(i,i),pch=15,col=rainbow(i),cex=pointsize)
  }
}

ggplot_colors<-function(){
  grDevices::colors()
}

#######################################################################################
#Convex hull around data points in a particular color (specified by linecolor).

plot_convexhull<-function(xcoord,ycoord,linecolor="black",linetype=1,linewidth=1,fillcolor="blue"){
  hpts<-chull(x=xcoord,y=ycoord)
  hpts<-c(hpts,hpts[1])
  #lines(xcoord[hpts],ycoord[hpts],col=linecolor,lty=linetype,lwd=linewidth)
  polygon(xcoord[hpts],ycoord[hpts],col=fillcolor,border=linecolor,lty=linetype,lwd=linewidth)
}  

#######################################################################################
#Find the element of a vector with the closest value.

closest<-function(myvalue,myvector){
  which(abs(myvalue-myvector)==min(abs(myvalue-myvector)))
}  

#######################################################################################
#Manhattan plot.

plot_manhattan<-function(
    chrnamevector=paste0("Chr",1:10),
    chrpositionvector=1:10,
    values=runif(10),
    abovethreshold=Inf,
    belowthreshold=-Inf,
    colors=c("blue","red"),
    xlabel="Chromosome",
    ylabel="Value",
    plottitle="Manhattan plot"){
  
  #Load libraries.
  library(ggplot2)
  library(tidyverse)
  
  #Plot theme.
  myggplottheme<-theme(title=element_text(size=10,face="bold"),
                       axis.title=element_text(size=10,face="bold"),
                       axis.text=element_text(size=10),
                       axis.text.x=element_text(angle=60,size=8,vjust=0.5),
                       legend.position="none",
                       legend.title=element_text(size=10,face="bold"),
                       legend.text=element_text(size=10),
                       legend.key=element_blank(),
                       panel.grid=element_line(colour="gray90"),
                       panel.grid.major.x=element_blank(),
                       panel.grid.minor.x=element_blank(),
                       panel.background=element_rect(fill="white",colour="black"),
                       panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                       strip.background=element_rect(colour="black",
                                                     fill="white"))
  
  #Format data.
  DATA<-data.frame(Chromosome=chrnamevector,
                   Position=chrpositionvector,
                   Value=values)
  
  #Define transparency (alpha value).
  if(abovethreshold<Inf | belowthreshold>-Inf){
    DATA$Alpha<-0.3
  }else{
    DATA$Alpha<-1
  }
  if(abovethreshold<Inf){
    DATA$Alpha[DATA$Value>=abovethreshold]<-1
  }
  if(belowthreshold>-Inf){
    DATA$Alpha[DATA$Value<=belowthreshold]<-1
  }
  
  #Define chromosome colors.
  #RECTCOLORS<-rep(c("gray","white"),length(unique(DATA$Chromosome)))
  CHRCOLORS<-rep(colors,length(unique(DATA$Chromosome)))
  COL<-0
  for(C in unique(DATA$Chromosome)){
    COL<-COL+1
    #DATA$RectColor[DATA$Chromosome==C]<-RECTCOLORS[COL]
    DATA$ChrColor[DATA$Chromosome==C]<-CHRCOLORS[COL]
  }
  
  #Calculate range of gray/white background.
  RECTCOLORS<-rep(c("gray","white"),length(unique(DATA$Chromosome)))
  RECTX<-data.frame(Rectxmin=rep(NA,length(unique(DATA$Chromosome))),
                    Rectxmax=NA,
                    Rectcolor=NA)
  for(C in 1:length(unique(DATA$Chromosome))){
    RECTX$Rectxmin[C]<-min(DATA$Position[DATA$Chromosome==unique(DATA$Chromosome)[C]])
    RECTX$Rectxmax[C]<-max(DATA$Position[DATA$Chromosome==unique(DATA$Chromosome)[C]])
    RECTX$Rectcolor[C]<-RECTCOLORS[C]
  }
  
  #Calculate center of each chromosome on the x axis.
  AXIS_SET<-DATA %>% 
    group_by(Chromosome) %>%
    summarize(center=mean(Position))
  
  #Plot.
  p<-ggplot(data=DATA)+
    geom_rect(data=RECTX,aes(xmin=Rectxmin,xmax=Rectxmax,ymin=-Inf,ymax=Inf),
              fill=RECTX$Rectcolor,
              alpha=0.3)+
    #geom_rect(aes(xmin=Position,xmax=dplyr::lead(Position),ymin=-Inf,ymax=Inf),
    #          fill=DATA$RectColor,
    #          alpha=0.3)+
    geom_point(aes(Position,Value),
               color=DATA$ChrColor,
               alpha=DATA$Alpha,
               shape=18,
               size=1)+
    scale_x_continuous(label=AXIS_SET$Chromosome,breaks=AXIS_SET$center)+
    scale_y_continuous(limits=c(min(DATA$Value,na.rm=TRUE),max(DATA$Value,na.rm=TRUE)))+
    geom_hline(yintercept=belowthreshold,color="black",linetype="dashed")+
    geom_hline(yintercept=abovethreshold,color="black",linetype="dashed")+
    labs(x=xlabel,y=ylabel)+
    ggtitle(plottitle)+
    myggplottheme
  
  return(p)
}

#######################################################################################
#Greek letters and mathematical symbols (plots).

char_alpha<-intToUtf8(945)
char_beta<-intToUtf8(946)
char_gamma<-intToUtf8(947)
char_delta<-intToUtf8(948)
char_epsilon<-intToUtf8(949)
char_zeta<-intToUtf8(950)
char_eta<-intToUtf8(951)
char_theta<-intToUtf8(952)
char_iota<-intToUtf8(953)
char_kappa<-intToUtf8(954)
char_lambda<-intToUtf8(955)
char_mu<-intToUtf8(956)
char_nu<-intToUtf8(957)
char_xi<-intToUtf8(958)
char_omicron<-intToUtf8(959)
char_pi<-intToUtf8(960)
char_rho<-intToUtf8(961)
char_sigma<-intToUtf8(963)
char_tau<-intToUtf8(964)
char_upsilon<-intToUtf8(965)
char_phi<-intToUtf8(966)
char_chi<-intToUtf8(967)
char_psi<-intToUtf8(968)
char_omega<-intToUtf8(969)
rightarrow<-intToUtf8(8594)
leftarrow<-intToUtf8(8592)
leftrightarrow<-intToUtf8(8596)
infinity<-intToUtf8(8734)
forall<-intToUtf8(8704)
partialdiff<-intToUtf8(8706)
thereexists<-intToUtf8(8707)
thereexistsnot<-intToUtf8(8708)
emptyset<-intToUtf8(8709)
elementof<-intToUtf8(8712)
elementofnot<-intToUtf8(8713)
contains<-intToUtf8(8715)
containsnot<-intToUtf8(8716)
plusorminus<-intToUtf8(177)
cdot<-intToUtf8(8729)
angle<-intToUtf8(8736)
logicaland<-intToUtf8(8743)
logicalor<-intToUtf8(8744)
integral<-intToUtf8(8747)
therefore<-intToUtf8(8756)
because<-intToUtf8(8757)
notsign<-intToUtf8(172)
lessthanorequal<-intToUtf8(8804)
greaterthanorequal<-intToUtf8(8805)
notequal<-intToUtf8(8800)

#######################################################################################
#



#######################################################################################
#



#######################################################################################
#



#######################################################################################
#



#######################################################################################
#


