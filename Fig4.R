#####################################################
#
#    Code for paper "Evaluation of the effects of vaccination regimes on the transmission dynamics of COVID-19 pandemic"
#    Author: Ichiro Nakamoto
#####################################################




library(reshape2)
library(ggplot2)
library(cowplot)
library(RConics)
library(latex2exp)
library(dplyr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(RConics)
library(latex2exp)
library(dplyr)

R0computation=function(nu,omega,omega1,mu,beta,gamma,rho1,rho2,rho3,epsilonV1,epsilonV2,epsilonV3,epsilon1,epsilon2,epsilon3){
  a=mu/(nu+mu)
  b=nu/(nu+mu)
  c=beta/(gamma+mu)
  d=1/(omega+rho1+mu)
  e=omega/(omega1+rho2+mu)
  f=omega1/(rho3+mu)
  
  R0=c*(a*(1+nu*epsilonV1*d+nu*epsilonV2*d*e+nu*epsilonV3*d*e*f)+b*(rho1*epsilon1*d+rho2*epsilon2*d*e+rho3*epsilon3*d*e*f))
  
  
  return(R0)
}



R0computation.vectorized=Vectorize(R0computation)

epsilon1.epsilon2.plots=function(omega1,mu,beta,gamma,rho2,rho3,epsilonV1,epsilonV2,epsilonV3,epsilon1,epsilon2,epsilon3){
  omega=seq(from=0,to=1,length.out=300)
  nu=seq(from=0,to=1,length.out=300)
  
  rho1.vector=c(1/(0.1*52),1/(0.25*52),1/(0.5*52))
  rho1.vector.labels=c("immunity of dose 1 wanes to suscep. in 5w","immunity of dose 1 wanes to suscep. in 13w","immunity of dose 1 wanes to suscep. in 26w")
  data.list=list()
  for (i in 1:length(rho1.vector)){
    data.list[[i]]=outer(X=nu,Y=omega,
                         FUN=R0computation.vectorized,
                         beta=3.4,gamma=7/5,mu=1/(50*52),nu=nu,
                         rho1=rho1.vector[i],epsilonV1=epsilonV1,
                         epsilonV2=epsilonV2,epsilonV3=epsilonV3,rho2=rho2,
                         rho3=rho3,epsilon1=epsilon1,epsilon2=epsilon2,epsilon3=epsilon3)
  }
  for.plot.list=list()
  for (i in 1:length(rho1.vector)){
    for.plot.list[[i]]=melt(data.list[[i]])
    for.plot.list[[i]]$Var1=omega[for.plot.list[[i]]$Var1]
    for.plot.list[[i]]$Var2=nu[for.plot.list[[i]]$Var2]
    names(for.plot.list)[[i]]=rho1.vector.labels[i]
  }
  
  for.plot=bind_rows(for.plot.list,.id=c("rho1"))
  for.plot$rho1=factor(for.plot$rho1,levels=c("immunity of dose 1 wanes to suscep. in 5w","immunity of dose 1 wanes to suscep. in 13w","immunity of dose 1 wanes to suscep. in 26w"))
  ggplot()+geom_tile(data=for.plot,aes(x=Var2,y=Var1,fill=value))+
    scale_fill_binned(name=expression(R[""]),low = "mediumpurple", high = "mediumblue")+
    facet_wrap(~rho1,ncol=3)+scale_x_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0))+scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0))+
    xlab(expression(" dose 1 rate(" * nu[] * ")"))+
    ylab(expression(atop("dose 2 rate(" * omega[] * ")")))+theme(legend.text=element_text(size=14))+
    theme(legend.margin=margin(c(0,0,0,0)))+theme(plot.margin=unit(c(0,0,1,0),"cm"))+theme(text=element_text(size=14))+
    labs(caption = "  ",plot.caption = element_text(hjust = 0.2)) 
  
}






output.plot.a0=epsilon1.epsilon2.plots(omega1=0.0000005,beta=3.4, rho2=1/(0.5*52),rho3=1/(0.5*52),gamma=1.4,mu=1/(50*52),
                                      epsilon1=0.6,epsilon2=0.5,epsilon3=0.5,
                                      epsilonV1=0.8,epsilonV2=0.7,epsilonV3=0.6)

output.plot.a=epsilon1.epsilon2.plots(omega1=0.0000005,beta=3.4, rho2=1/(0.5*52),rho3=1/(0.5*52),gamma=1.4,mu=1/(50*52),
                                         epsilon1=0.5,epsilon2=0.2,epsilon3=0.1,
                                         epsilonV1=0.4,epsilonV2=0.1,epsilonV3=0.08)


output.plot.b=epsilon1.epsilon2.plots(omega1=0.0000005,beta=3.4, rho2=1/(1*52),rho3=1/(1*52),gamma=1.4,mu=1/(50*52),
                                         epsilon1=0.4,epsilon2=0.2,epsilon3=0.1,
                                         epsilonV1=0.3,epsilonV2=0.05,epsilonV3=0.05)

output.plot.c=epsilon1.epsilon2.plots(omega1=0.0000005,beta=3.4, rho2=1/(1.5*52),rho3=1/(1.5*52),gamma=1.4,mu=1/(50*52),
                                         epsilon1=0.4,epsilon2=0.2,epsilon3=0.1,
                                         epsilonV1=0.2,epsilonV2=0.03,epsilonV3=0.02)



combined=plot_grid(output.plot.a0+theme(plot.margin=unit(c(0.5,1,0,1),"cm"))+
                     ggtitle(expression(paste("severe variants"))),
                   output.plot.a+theme(plot.margin=unit(c(0.5,1,0,1),"cm"))+
                            ggtitle(expression(paste(rho[2]==3*rho[0]))),
                   output.plot.b+theme(plot.margin=unit(c(0.5,1,0,1),"cm"))+
                     ggtitle(expression(paste(rho[2]==1.5*rho[0]))),
                  output.plot.c+theme(plot.margin=unit(c(0.5,1,1,1),"cm"))+
                    ggtitle(expression(paste(rho[2]==rho[0]))),
                   nrow=4,labels = "AUTO",label_fontface="plain")+ 
      draw_figure_label(label = "Fig.4  Impact of vaccination rate on the reproduction number of COVID-19. \n
                                (A) Case where severe variant occurs; (B) Immunity of dose 2&3 wanes in 26 weeks;  \n
                                (C) Immunity of dose 2&3 wanes in 52 weeks;  (D) Immunity of dose 2&3 wanes in 78 weeks.",
                         size = 12,position = "bottom.left",lineheight = 0.4)





pdf("Fig4.pdf",width=12,height=15)
#ggsave("Fig4.svg", width = 12, height = 10)
print(combined)
dev.off()

