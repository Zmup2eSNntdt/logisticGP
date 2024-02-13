#######################
# nmcmc=25000
# burnin=30000
# thining=5
# # run1<- get_updates(xi=xi,beta = beta, tau=tau,h_loc=h_loc, X=X2,y=y1,
# #                    test_knots = test_knots,
# #                    verbose = T, burnin = burnin,nmcmc = nmcmc,thining = thining)
# 
# XI_mat = run1$XI_mat.out[,,2001:5000]
# TAU_mat = run1$TAU_mat.out[,,2001:5000]
# BETA_mat = run1$BETA_mat.out[,,2001:5000]
# 
# n <- nrow(X) #no. of individuals
# K <- ncol(y) #no. of Location
# p <- ncol(X) #no. of covariates
#######################
#save("run1", file = "parameters_beta_grid_1001.RData")
#######################
library(reshape2)
library(ggplot2)
library(ggridges)

set.seed(738)
grid<- seq(0,1,length.out = 10001)
par(mfrow=c(3,3))
location.1<-sort(sample(1:K,4,replace = F)) # select some locations randomly
density_mat<- matrix(0,nr=length(grid),nc=length(location.1))
ind.index = sort(sample(1:nrow(X_test),2,replace = F))
X_new=matrix(X2[ind.index,],nc=p)
for(i in 1:length(location.1))
{

  density.array<- PostFun_loc(XI_mat, grid, X_new, BETA_mat, TAU_mat,
                              location = location.1[i],
                              test_knots = test_knots)
  
  density<- density.array[,,2]
  density_q2<- apply(density,1,function(x) {quantile(x,0.50)})
  density_mat[,i]<- density_q2/sum(density_q2*(grid[2]-grid[1]))
  #print(i)
}


density.df0 = data.frame(cbind(grid,density_mat))
for (i in 2:ncol(density.df0)) {
  names(density.df0)[i] = paste('',location.1[(i-1)])
}
density.df1 = melt(density.df0,id.vars=grid)

# My_Theme2 =  theme(#axis.title.x=element_blank(),
#   axis.text.x=element_text(size = 20,colour = "black"),
#   #axis.ticks.x=element_blank(),
#   axis.text.y = element_text(size = 20,colour = "black"),
#   axis.title.y = element_text(size = 32), 
#   axis.title.x = element_text(size = 32),
#   plot.title = element_text(size=32,hjust = 0.5, face="bold"),
#   legend.title = element_blank(), legend.text = element_text(size = 32))

# ggplot(data = density.df1, aes(x=grid, y=variable, height = value, group= variable)) +
#   geom_density_ridges(stat = 'identity',scale=1.5) + theme_classic() + My_Theme2 + 
#   labs(title = 'Conditional density estimates for a given individual',
#        y = 'location')












