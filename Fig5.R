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

R0computation=function(omega,omega1,nu,mu,beta,gamma,rho1,rho2,rho3,epsilonV1,epsilonV2,epsilonV3,epsilon1,epsilon2,epsilon3){
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

epsilon1.epsilon2.plots=function(nu,mu,beta,gamma,rho1,rho3,epsilonV1,epsilonV2,epsilonV3,epsilon1,epsilon2,epsilon3){
  omega=seq(from=0,to=0.2,length.out=300)
  omega1=seq(from=0,to=0.2,length.out=300)
  rho2.vector=c(1/(0.25*52),1/(0.5*52),1/(1*52))
  rho2.vector.labels=c("immunity of dose 2 wanes to suscep. in 13w","immunity of dose 2 wanes to suscep. in 26w","immunity of dose 2 wanes to suscep. in 52w")
  data.list=list()
  for (i in 1:length(rho2.vector)){
    data.list[[i]]=outer(X=omega,Y=omega1,
                         FUN=R0computation.vectorized,
                         beta=3.4,gamma=1.4,mu=1/(50*52),nu=nu,
                         rho2=rho2.vector[i],epsilonV1=epsilonV1,
                         epsilonV2=epsilonV2,epsilonV3=epsilonV3,rho1=rho1,
                         rho3=rho3,epsilon1=epsilon1,epsilon2=epsilon2,epsilon3=epsilon3)
  }
  for.plot.list=list()
  for (i in 1:length(rho2.vector)){
    for.plot.list[[i]]=melt(data.list[[i]])
    for.plot.list[[i]]$Var1=omega[for.plot.list[[i]]$Var1]
    for.plot.list[[i]]$Var2=omega1[for.plot.list[[i]]$Var2]
    names(for.plot.list)[[i]]=rho2.vector.labels[i]
  }
  
  for.plot=bind_rows(for.plot.list,.id=c("rho2"))
  for.plot$rho2=factor(for.plot$rho2,levels=c("immunity of dose 2 wanes to suscep. in 13w","immunity of dose 2 wanes to suscep. in 26w","immunity of dose 2 wanes to suscep. in 52w"))
  ggplot()+geom_tile(data=for.plot,aes(x=Var1,y=Var2,fill=value))+
    scale_fill_binned(name=expression(R[""]),n.breaks = 3, low = "plum", high = "blue")+
    facet_wrap(~rho2,ncol=3)+scale_x_continuous(breaks=c(0,0.05,0.1,0.15,0.2))+scale_y_continuous(breaks=c(0,0.05,0.1,0.15,0.2))+
    xlab(expression(" dose 2 vaccination rate (" * omega[] * ")"))+
    ylab(expression(atop("dose 3 vaccination rate (" * omega[1] * ")")))+theme(legend.text=element_text(size=14))+
    theme(legend.margin=margin(c(0,0,0,0)))+theme(plot.margin=unit(c(0,0,1,0),"cm"))+theme(text=element_text(size=14))+
    labs(caption = "  ",plot.caption = element_text(hjust = 0.2)) 
}





output.plot.a=epsilon1.epsilon2.plots(nu=0.01,beta=3.4, rho1=1/(0.5*52),rho3=1/(1*52),gamma=1.4,mu=1/(50*52),
                                         epsilon1=0.5,epsilon2=0.3,epsilon3=0.1,
                                         epsilonV1=0.1,epsilonV2=0.05,epsilonV3=0.03)


output.plot.b=epsilon1.epsilon2.plots(nu=0.015,beta=3.4, rho1=1/(0.5*52),rho3=1/(1*52),gamma=1.4,mu=1/(50*52),
                                         epsilon1=0.5,epsilon2=0.3,epsilon3=0.08,
                                         epsilonV1=0.1,epsilonV2=0.05,epsilonV3=0.02)

output.plot.c=epsilon1.epsilon2.plots(nu=0.02,beta=3.4, rho1=1/(0.5*52),rho3=1/(1*52),gamma=1.4,mu=1/(50*52),
                                         epsilon1=0.5,epsilon2=0.3,epsilon3=0.05,
                                         epsilonV1=0.1,epsilonV2=0.05,epsilonV3=0.01)



combined=plot_grid(output.plot.a+theme(plot.margin=unit(c(0,1,0,1),"cm"))+
                            ggtitle(expression(paste(nu[]==nu[0]))),
                   output.plot.b+theme(plot.margin=unit(c(0,1,0,1),"cm"))+
                     ggtitle(expression(paste(nu[]==1.5*nu[0]))),
                  output.plot.c+theme(plot.margin=unit(c(0,1,0.5,1),"cm"))+
                    ggtitle(expression(paste(nu[]==2*nu[0]))),
                   nrow=3,labels = "AUTO",label_fontface="plain",hjust=-0.9)+
  draw_figure_label(label = "Fig.5  Impact of vaccination rate interaction on the dynamic reproduction number of COVID-19 pandemic. \n
                    In all scenarios, immunity of dose2&3 wanes to susceptibility in 26 and 52 weeks respectively.",
                         size = 12,position = "bottom.left",lineheight = 0.4)

pdf("Fig5.pdf",width=12,height=12)
#ggsave("Fig5.svg", width = 12, height = 12)
print(combined)
dev.off()

