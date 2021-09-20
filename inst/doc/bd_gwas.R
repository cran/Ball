## ---- message=FALSE, warning=FALSE--------------------------------------------
library(mvtnorm)

num <- 100
snp_num <- 200
k <- 100
rho <- 0.5
freq0 <- 0.75
d <- 3

set.seed(2021)

ar1 <- function (p, rho = 0.5) 
{
    Sigma <- matrix(0, p, p)
    for (i in 1:p) {
        for (j in 1:p) {
            Sigma[i, j] <- rho^(abs(i - j))
        }
    }
    return(Sigma)
}

mean0 <- rep(0, k)
mean1 <- rep(0.1 * d, k)
mean2 <- rep(-0.1 * d, k)

cov0 <- ar1(p = k, rho = rho)
cov1 <- ar1(p = k, rho = rho - 0.1 * d)
cov2 <- ar1(p = k, rho = rho + 0.1 * d)

p1 <- freq0 ^ 2
p2 <- 2 * freq0 * (1 - freq0)
n1 <- round(num * p1)
n2 <- round(num * p2)
n3 <- num - n1 - n2
x0 <- rmvnorm(n1, mean = mean0, sigma = cov0)
x1 <- rmvnorm(n2, mean = mean1, sigma = cov1)
x2 <- rmvnorm(n3, mean = mean2, sigma = cov2)
x <- rbind(x0, x1, x2)
head(x[, 1:6])

## -----------------------------------------------------------------------------
effect_snp <- c(rep(0, n1), rep(1, n2),
                rep(2, n3))
noise_snp <- sapply(2:snp_num, function(j) {
  sample(
    0:2,
    size = num,
    replace = TRUE,
    prob = c(p1, p2, 1 - p1 - p2)
  )
})
snp <- cbind(effect_snp, noise_snp)
head(snp[, 1:6])

## -----------------------------------------------------------------------------
library(Ball)
res <- bd.gwas.test(x = x, snp = snp)

## -----------------------------------------------------------------------------
str(res)

## ---- echo=FALSE, fig.align='center', eval=FALSE------------------------------
#  library(ggplot2)
#  library(ggpubr)
#  size_15 <- 9
#  size_20 <- 12
#  
#  df_1_75<-data.frame(d=c(1:5,1:5,1:5),
#                 power=c(0,0.98,1,1,1,0,1,1,1,1,0.04,1,1,1,1),
#                 group=c(rep("n=500",5),rep("n=750",5),rep("n=1000",5)))
#  df_1_75$group<-factor(df_1_75$group, levels = c("n=500", "n=750", "n=1000"))
#  
#  p1 <- ggplot(data=df_1_75, aes(x=d, y=power, colour=group,linetype=group)) +
#    geom_point()+
#    geom_line(size=1) +
#    scale_linetype_manual(values = c('dotdash', 'dotted', 'dashed'))+
#    labs(x="d",y="power",title = "power curve \n rho=0.5,p=0.75",
#         caption = "Figure 1")+
#    labs(colour="sample size",linetype="sample size")+
#    theme_classic()+
#    theme(plot.title = element_text(hjust = 0.5,size=size_20),
#          plot.caption = element_text(hjust = 0.5,size=size_15),
#          axis.title.x =element_text(size=size_15),
#          axis.title.y=element_text(size=size_15),
#          axis.text.x = element_text(size=size_15),
#          axis.text.y = element_text(size=size_15),
#          legend.title = element_text(size=size_15),
#          legend.text = element_text(size=size_15),
#          legend.position = "bottom")
#  
#  
#  df_1_95<-data.frame(d=c(1:5,1:5,1:5),
#                      power=c(0.01,0.03,0.5,1,1,0,0.15,0.81,1,1,0.01,0.21,0.99,1,1),
#                      group=c(rep("n=500",5),rep("n=750",5),rep("n=1000",5)))
#  df_1_95$group<-factor(df_1_95$group, levels = c("n=500", "n=750", "n=1000"))
#  
#  p2 <- ggplot(data=df_1_95, aes(x=d, y=power, colour=group,linetype=group)) +
#    geom_point()+
#    geom_line(size=1) +
#    scale_linetype_manual(values = c('dotdash', 'dotted', 'dashed'))+
#    labs(x="d",y="power", title = "power curve \n rho=0.5,p=0.95",
#         caption = "Figure 2")+
#    labs(colour="sample size",linetype="sample size")+
#    theme_classic()+
#    theme(plot.title = element_text(hjust = 0.5,size=size_20),
#          plot.caption = element_text(hjust = 0.5,size=size_15),
#          axis.title.x =element_text(size=size_15),
#          axis.title.y=element_text(size=size_15),
#          axis.text.x = element_text(size=size_15),
#          axis.text.y = element_text(size=size_15),
#          legend.title = element_text(size=size_15),
#          legend.text = element_text(size=size_15),
#          legend.position = "bottom")
#  
#  df_2_75<-data.frame(d=c(1:5,1:5,1:5),
#                      power=c(0,1,1,1,1,0,1,1,1,1,0,1,1,1,1),
#                      group=c(rep("n=500",5),rep("n=750",5),rep("n=1000",5)))
#  df_2_75$group<-factor(df_2_75$group, levels = c("n=500", "n=750", "n=1000"))
#  
#  p3 <- ggplot(data=df_2_75, aes(x=d, y=power, colour=group,linetype=group)) +
#    geom_point()+
#    geom_line(size=1) +
#    scale_linetype_manual(values = c('dotdash', 'dotted', 'dashed'))+
#    labs(x="d",y="power",title = "power curve \n rho=0,p=0.75",
#         caption = "Figure 3")+
#    labs(colour="sample size",linetype="sample size")+
#    theme_classic()+
#    theme(plot.title = element_text(hjust = 0.5,size=size_20),
#          plot.caption = element_text(hjust = 0.5,size=size_15),
#          axis.title.x =element_text(size=size_15),
#          axis.title.y=element_text(size=size_15),
#          axis.text.x = element_text(size=size_15),
#          axis.text.y = element_text(size=size_15),
#          legend.title = element_text(size=size_15),
#          legend.text = element_text(size=size_15),
#          legend.position = "bottom")
#  
#  
#  
#  df_2_95<-data.frame(d=c(1:5,1:5,1:5),
#                      power=c(0,0.02,0.56,1,1,0,0.04,0.93,1,1,0,0.11,1,1,1),
#                      group=c(rep("n=500",5),rep("n=750",5),rep("n=1000",5)))
#  df_2_95$group<-factor(df_2_95$group, levels = c("n=500", "n=750", "n=1000"))
#  
#  p4 <- ggplot(data=df_2_95, aes(x=d, y=power, colour=group,linetype=group)) +
#    geom_point()+
#    geom_line(size=1) +
#    scale_linetype_manual(values = c('dotdash', 'dotted', 'dashed'))+
#    labs(x="d",y="power",title = "power curve \n rho=0,p=0.95",
#         caption = "Figure 4")+
#    labs(colour="sample size",linetype="sample size")+
#    theme_classic()+
#    theme(plot.title = element_text(hjust = 0.5,size=size_20),
#          plot.caption = element_text(hjust = 0.5,size=size_15),
#          axis.title.x =element_text(size=size_15),
#          axis.title.y=element_text(size=size_15),
#          axis.text.x = element_text(size=size_15),
#          axis.text.y = element_text(size=size_15),
#          legend.title = element_text(size=size_15),
#          legend.text = element_text(size=size_15),
#          legend.position = "bottom")
#  
#  p <- ggarrange(p1, p2, p3, p4, nrow = 2, widths = 8, heights = 8,
#                 common.legend = TRUE, ncol = 2, legend = "bottom")
#  ggexport(p, filename = "kbd_gwas.png")

