library(latex2exp)
ind.index = sort(sample(1:nrow(X_test),nrow(X_test),replace = F))
X_pred=X_test[ind.index,]
n=nrow(X_pred)
Niter=nmcmc/thining
Median_mat<- array(0,dim = c(n,K,Niter))
#u<- runif(100)
#u<-c(0.15,0.30,0.45,0.6,0.75,.90)
u<-c(0.20,.80)
b<-rep(0,length(u))
x.grid<- seq(0,1,length.out = 100001)
Median_mat.loc = array(0,dim = c(n,K,length(u)))
Median_mat.loc.th = array(0,dim = c(n,K,length(u)))
xi=apply(XI_mat,MARGIN = c(1,2),mean)
tau=apply(TAU_mat,MARGIN = c(1,2),mean)
beta=apply(BETA_mat,MARGIN = c(1,2),mean)
for(j in 1: length(u)){
  
  Median_mat.loc[,,j] = get_quantile(xi = xi,tau = tau,beta = beta,
                                     test_knots = test_knots,
                                     prob = u[j],X_pred = X_pred)
}
#Median_mat.loc

# My_Theme2 =  theme(#axis.title.x=element_blank(),
#   axis.text.x=element_text(size = 20,colour = "black"),
#   #axis.ticks.x=element_blank(),
#   axis.text.y = element_text(size = 20,colour = "black"),
#   axis.title.y = element_text(size = 32), 
#   axis.title.x = element_text(size = 32),
#   plot.title = element_text(size=32,hjust = 0.5, face="bold"),
#   legend.title = element_blank(), legend.text = element_text(size = 32),legend.position = 'none')

A = (y_test>Median_mat.loc[,,1] & y_test<Median_mat.loc[,,2])
A[A]=1
b = rowMeans(A)
which(b>0.80)
l = min(which(b>0.80))
data.patient = data.frame(location = 1:83,step1 = Median_mat.loc[l,,1], 
                          step2= Median_mat.loc[l,,2],
                          obs =y_test[l,])
data.patient$col = as.factor(A[l,])

# ggplot(data.patient)+geom_step(aes(x=location-0.75,y=step1),size=3.5,col=4) +
#   geom_step(aes(x=location-0.75,y=step2),size=3.5,col=4) + ylim(c(0.3,0.9)) +
#   geom_point(aes(x=location,y=obs,col=col),size=4) + scale_color_manual(values = c(2,3)) +
#   theme_classic() + My_Theme2 + labs(title = "Detection of abnormalities for a given individual",y=TeX('FA'),x='location')

