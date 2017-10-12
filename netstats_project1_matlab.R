###visualize netstates###
netstats<-readMat('./wb_stats.mat')
btn_cent = mean(netstats$wb.bt)
ccoef = mean(netstats$ccoef)
eglob = mean(netstats$wb.eglob)
stats.df <- data.frame(btn_cent = apply(netstats$null.bt,2,mean),
                       ccoef = apply(netstats$null.ccoef,2,mean),
                       eglob = netstats$null.eglob[1,])
stats.df.mt <- melt(stats.df)
names(stats.df.mt) <- c('netstat','value')
ggplot(stats.df.mt[1:100,],aes(netstat,value))+ 
                                geom_boxplot() + 
                               geom_point(aes(y=btn_cent),,color = 'red')
ggplot(stats.df.mt[101:200,],aes(netstat,value))+ 
  geom_boxplot() + 
  geom_point(aes(y=ccoef),color = 'red')
ggplot(stats.df.mt[201:300,],aes(netstat,value))+ 
  geom_boxplot() + 
  geom_point(aes(y=eglob),color = 'red')
