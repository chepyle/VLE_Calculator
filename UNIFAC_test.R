# collection of UNIFAC structural groups and parameters
# found here: http://www.aim.env.uea.ac.uk/aim/info/UNIFACgroups.html
# used data from these sources:
# H. K. Hansen, P. Rasmussen, A. Fredenslund, M. Schiller, and J Gmehling (1991) Ind. Eng. Chem. Res. 30, 2352-2355.
# R. Wittig, J. Lohmann, and J. Gmehling (2003) Ind. Eng. Chem. Res. 42, 183-188.
# K. Balslev and J. Abildskov (2002) Ind. Eng. Chem. Res. 41, 2047-205.
# C. Peng, M. N. Chan, and C. K. Chan (2001) Environ. Sci. Technol. 35, 4495-4501


#checked on example 8-14 in Properties of Gases and Liquids, Fifth Edition (Bruce E. Poling, John M. Prausnitz, John P. O’Connell)


source('UNIFAC_calc.R')
library(reshape2)
library(ggplot2)
#test case for UNIFAC calculation:
P<-1.022 # bar
molecules=c('water','methanol')

## txy plot:
txy<-melt(Txy(P,molecules),
          id.vars = c('Temperature','x2','y2'))
p<-ggplot(txy,aes(x=value,y=Temperature,color=variable))+geom_line()+geom_point()
print(p)

## xyplot
xy<-melt(Txy(P,molecules),id.vars=c('x1','y1'),value.name=('Temperature'))
p<-ggplot(xy,aes(x=x1,y=y1,color=Temperature))+geom_line()+geom_point()+geom_abline(slope=1,intercept=0)
print(p)

P<-1.022 # bar
molecules=c('chloroform','methanol')

## txy plot:
txy<-Txy(P,molecules)
txy<-melt(Txy(P,molecules),
          id.vars = c('Temperature','x2','y2'))
p<-ggplot(txy,aes(x=value,y=Temperature,color=variable))+geom_line()+geom_point()
print(p)

## xyplot
xy<-melt(Txy(P,molecules),id.vars=c('x1','y1'),value.name=('Temperature'))
p<-ggplot(xy,aes(x=x1,y=y1,color=Temperature))+geom_line()+geom_point()+geom_abline(slope=1,intercept=0)
print(p)
