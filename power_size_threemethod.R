
TEST_mst<-function(k,data1,data2){
  MST1<-list()
  TEST<-list()
  data<-rbind(data1,data2)
  n1<-nrow(data1)
  n<-nrow(data)
  distance<-distance(data,method="euclidean")
  mst_count1<-matrix(0,n,2)
  mst_count1[1:n1,1]<-1
  mst_count1[(n1+1):n,2]<-1
  mydist<-distance(data)
  for (j in 1:k) {
  MST1[[j]]<-getGraph(mst_count1, mydist, j, graph.type = "mstree")
  TEST[[j]]<-g.tests(MST1[[j]], 1:n1, (n1+1):n, test.type="all", maxtype.kappa = 1.14, perm=0)
  }
  return(TEST)
}


sta_three<-function(k){
  times=5000
sta_1979 <- matrix(0,times,k)
sta_general<- matrix(0,times,k)
sta_weighted<- matrix(0,times,k)
Sigma<-diag(2)
for (i in 1:times) {
  data1<-mvrnorm(81, rep(0, 2), Sigma)
  data2<-mvrnorm(36, rep(0, 2), Sigma)
  TEST_T<-TEST_mst(k,data1,data2)
  for (l in 1:k){
  sta_1979[i,l] <-TEST_T[[l]]$original$test.statistic
  sta_general[i,l]<-TEST_T[[l]]$generalized$test.statistic
  sta_weighted[i,l]<-TEST_T[[l]]$weighted$test.statistic
}
}
return(list(sta_1979,sta_general,sta_weighted))
}

statistic_three<-sta_three(10)
#cv_1979 <- quantile(sta_1979,c(0.025,0.975))
cv_1979_05=cv_general_95=cv_weighted_95=rep(0,k)
for (j in 1:k) {
  cv_1979_05[j]<-quantile(statistic_three[[1]][,j],0.05)
  #cv_general <- quantile(sta_general,c(0.025,0.975))
  cv_general_95[j]<-quantile(statistic_three[[2]][,j],0.95)
  #cv_weighted <- quantile(sta_weighted,c(0.025,0.975))
  cv_weighted_95[j]<-quantile(statistic_three[[3]][,j],0.95)
}

three_size_power<-function(k,statistic_three){
#statistic_three<-sta_three(k)
#cv_1979 <- quantile(sta_1979,c(0.025,0.975))
cv_1979_05=cv_general_95=cv_weighted_95=rep(0,k)
for (j in 1:k) {
cv_1979_05[j]<-quantile(statistic_three[[1]][,j],0.05)
#cv_general <- quantile(sta_general,c(0.025,0.975))
cv_general_95[j]<-quantile(statistic_three[[2]][,j],0.95)
#cv_weighted <- quantile(sta_weighted,c(0.025,0.975))
cv_weighted_95[j]<-quantile(statistic_three[[3]][,j],0.95)
}
##power of MST
##不同分布
times0=500
sta_power_1979 <-matrix(0,times0,k)
sta_power_general <- matrix(0,times0,k)
sta_power_weighted<- matrix(0,times0,k)
sta_size_1979 <- matrix(0,times0,k)
sta_size_general <- matrix(0,times0,k)
sta_size_weighted <- matrix(0,times0,k)
size05_1979<-rep(0,k)
size05_general<-rep(0,k)
size05_weighted<-rep(0,k)
power05_197<-rep(0,k)
power05_general<-rep(0,k)
power05_weighted<-rep(0,k)
Sigma<-diag(3)
for (i in 1:times0) {
  data1<-mvrnorm(81, rep(0, 3), Sigma)
  data2<-mvrnorm(36, rep(0, 3), Sigma)
  data3<-mvrnorm(81, rep(0, 3), Sigma)
  data4<-mvrnorm(36, rep(0, 3), 2*Sigma)
  #data4[,2]<-rnorm(36,0,2)
  MST<-TEST_mst(k,data3,data4)
  MST_SIZE<-TEST_mst(k,data1,data2)
  for (j in 1:k) {
  sta_power_1979[i,j]<-MST[[j]]$original$test.statistic
  sta_power_general[i,j]<-MST[[j]]$generalized$test.statistic
  sta_power_weighted[i,j]<-MST[[j]]$weighted$test.statistic
  sta_size_1979[i,j]<-MST_SIZE[[j]]$original$test.statistic
  sta_size_general[i,j]<-MST_SIZE[[j]]$generalized$test.statistic
  sta_size_weighted[i,j]<-MST_SIZE[[j]]$weighted$test.statistic
  }
  for (l in 1:k) {
    size05_1979[l]<- length(which(sta_size_1979[,l]<cv_1979_05[l]))/times0
    size05_general[l] <- length(which(sta_size_general[,l]>cv_general_95[l]))/times0
    #power_weighted <- length(which(sta_power_weighted>cv_weighted[2]|sta_power_weighted<cv_weighted[1]))/times0
    size05_weighted[l] <- length(which(sta_size_weighted[,l]>cv_weighted_95[l]))/times0
    power05_1979[l] <- length(which(sta_power_1979[,l]<cv_1979_05[l]))/times0
    #power_general <- length(which(sta_power_general>cv_general[2]|sta_power_general<cv_general[1]))/times0
    power05_general[l] <- length(which(sta_power_general[,l]>cv_general_95[l]))/times0
    #power_weighted <- length(which(sta_power_weighted>cv_weighted[2]|sta_power_weighted<cv_weighted[1]))/times0
    power05_weighted[l] <- length(which(sta_power_weighted[,l]>cv_weighted_95[l]))/times0
  }
}
return(cbind(size05_1979,size05_general,size05_weighted,power05_1979,power05_general,power05_weighted))
}
power_size_3_type1<-three_size_power(10,statistic_three)
power_size_3_type2<-three_size_power(10,statistic_three)
power_size_3_type3<-three_size_power(10,statistic_three)
power_size_3_type4<-three_size_power(10,statistic_three)
size_4<-matrix(0,10,4)
size_4[,1:3]<-power_size_3[,1:3]
power_size<-matrix(0,10,2)
power_size[,1]<-power_size[[1]]
power_size[,2]<-power_size[[2]]
size_4[,4]<-power_size[,1]
power_4<-matrix(0,10,4)
power_4[,1:3]<-power_size_3[,4:6]
power_4[,4]<-power_size[,2]

colnames(size_4)<-c("1979","general_2017","weighted_2018","cross_rank")
colnames(power_4)<-c("1979","general_2017","weighted_2018","cross_rank")
plot_size<-melt(size_4)
plot_size$k=plot_size$Var1
plot_power<-melt(power_4)
plot_power$k<-plot_power$Var1
ggplot(plot_size, aes(x=k, y=value)) + geom_line(aes(color=Var2))+geom_point(size=3,aes(shape=Var2,colour=Var2))+labs(y="size",fill="shapes")+ggtitle("size of four tests")+ theme(plot.title = element_text(hjust = 0.5,size=23))+theme(legend.position=c(0.87,0.85))+theme(legend.title =element_blank())+theme(legend.key.size = unit(23, "pt"))+ theme(legend.text=element_text(size=20))+ theme(axis.title=element_text(size=20),axis.title.x =element_text(size=20), axis.title.y=element_text(size=20),axis.text.y=element_text(size=20),axis.text.x=element_text(size=20))
ggplot(plot_power, aes(x=k, y=value)) + geom_line(aes(color=Var2))+geom_point(size=3,aes(shape=Var2,colour=Var2))+labs(y="power",fill="shapes")+ggtitle("power of four tests")+ theme(plot.title = element_text(hjust = 0.5,size=23))+theme(legend.position="none")+theme(legend.title =element_blank())+theme(legend.key.size = unit(23, "pt"))+ theme(legend.text=element_text(size=20))+ theme(axis.title=element_text(size=20),axis.title.x =element_text(size=20), axis.title.y=element_text(size=20),axis.text.y=element_text(size=20),axis.text.x=element_text(size=20))

power_type2<-power_type3<-matrix(0,10,4)
power_type4<-matrix(0,10,4)
power_type2[,1:3]<-power_size_3_type2[,4:6]
power_type2[,4]<-power_size_type2[[2]]
power_type3[,1:3]<-power_size_3_type3[,4:6]
power_type3[,4]<-power_size_type3[[2]]
power_type4<-matrix(0,10,4)
power_type4[,1:3]=power_size_3_type4[,4:6]
power_type4[,4]<-power_size_type4[[2]]
colnames(power_type2)<-c("1979","general_2017","weighted_2018","cross_rank")
colnames(power_type3)<-c("1979","general_2017","weighted_2018","cross_rank")
colnames(power_type4)<-c("1979","general_2017","weighted_2018","cross_rank")
plot_power_type3<-melt(power_type3)
plot_power_type3$k<-plot_power_type3$Var1
ggplot(plot_power_type3, aes(x=k, y=value)) + geom_line(aes(color=Var2))+geom_point(size=3,aes(shape=Var2,colour=Var2))+labs(y="power",fill="shapes")+ggtitle("power of four tests")+ theme(plot.title = element_text(hjust = 0.5,size=23))+theme(legend.position="none")+theme(legend.title =element_blank())+theme(legend.key.size = unit(23, "pt"))+ theme(legend.text=element_text(size=20))+ theme(axis.title=element_text(size=20),axis.title.x =element_text(size=20), axis.title.y=element_text(size=20),axis.text.y=element_text(size=20),axis.text.x=element_text(size=20))

plot_power_type4<-melt(power_type4)
plot_power_type4$k<-plot_power_type4$Var1
ggplot(plot_power_type4, aes(x=k, y=value)) + geom_line(aes(color=Var2))+geom_point(size=3,aes(shape=Var2,colour=Var2))+labs(y="power",fill="shapes")+ggtitle("power of four tests")+ theme(plot.title = element_text(hjust = 0.5,size=23))+theme(legend.position="none")+theme(legend.title =element_blank())+theme(legend.key.size = unit(23, "pt"))+ theme(legend.text=element_text(size=20))+ theme(axis.title=element_text(size=20),axis.title.x =element_text(size=20), axis.title.y=element_text(size=20),axis.text.y=element_text(size=20),axis.text.x=element_text(size=20))

plot_power_type2<-melt(power_type2)
plot_power_type2$k<-plot_power_type2$Var1
ggplot(plot_power_type2, aes(x=k, y=value)) + geom_line(aes(color=Var2))+geom_point(size=3,aes(shape=Var2,colour=Var2))+labs(y="power",fill="shapes")+ggtitle("power of four tests")+ theme(plot.title = element_text(hjust = 0.5,size=23))+theme(legend.position="none")+theme(legend.title =element_blank())+theme(legend.key.size = unit(23, "pt"))+ theme(legend.text=element_text(size=20))+ theme(axis.title=element_text(size=20),axis.title.x =element_text(size=20), axis.title.y=element_text(size=20),axis.text.y=element_text(size=20),axis.text.x=element_text(size=20))


TEST_mst_netdata<-function(k,data1,data2){
  MST1<-list()
  TEST<-list()
    data_dif<-array(0,dim = c(dim(data1)[1]+dim(data2)[1],dim(data2)[2],dim(data2)[2]))
    data_dif[1:dim(data1)[1],,]<-data1
    data_dif[(dim(data1)[1]+1):(dim(data_dif)[1]),,]<-data2
    n1<-dim(data1)[1]
    n<-dim(data_dif)[1]
    mst_count1<-matrix(0,n,2)
    mst_count1[1:n1,1]<-1
    mst_count1[(n1+1):n,2]<-1
    mydist<-netdis_matrix(data_dif)
  for (j in 1:k) {
    MST1[[j]]<-getGraph(mst_count1, mydist, j, graph.type = "mstree")
    TEST[[j]]<-g.tests(MST1[[j]], 1:n1, (n1+1):n, test.type="all", maxtype.kappa = 1.14, perm=0)
  }
  return(TEST)
}

p_value_three<-function(k){
  p_value_1979=rep(0,k)
  p_value_general=rep(0,k)
  p_value_weighted=rep(0,k)
  for (j in 1:k) {
    TEST<-TEST_mst_netdata(k,am_array,am_con_array)
    p_value_1979[j]<-TEST[[j]]$original$pval.approx
    p_value_general[j]<-TEST[[j]]$generalized$pval.approx
    p_value_weighted[j]<-TEST[[j]]$weighted$pval.approx
  }
  return(cbind(p_value_1979,p_value_general,p_value_weighted))
}

#
p_three<-p_value_three(k)

T_sta_rank<-function(k,data1,data2){
  data_dif<-array(0,dim = c(dim(data1)[1]+dim(data2)[1],dim(data2)[2],dim(data2)[2]))
  data_dif[1:dim(data1)[1],,]<-data1
  data_dif[(dim(data1)[1]+1):(dim(data_dif)[1]),,]<-data2
  n1<-dim(data1)[1]
  n<-dim(data_dif)[1]
  mst_count1<-matrix(0,n,2)
  mst_count1[1:n1,1]<-1
  mst_count1[(n1+1):n,2]<-1
  mydist<-netdis_matrix(data_dif)
  MST1<-getGraph(mst_count1, mydist, k, graph.type = "mstree")
  MST2<-as.numeric(MST1)
  MST3<-matrix(0,k*(n-1),2)
  rule<-as.factor(rep(rep(c(1,2),each=n-1),k))
  group<-split(MST2,rule)
  MST3[,1]<-unlist(group[1])
  MST3[,2]<-unlist(group[2])
  MST4<-as.data.frame(MST3)
  MST4$dis<-0
  for (i in 1:n-1) {
    MST4$dis[i]<-mydist[MST4$V1[i],MST4$V2[i]]
  }
  MST4$order <- rank(MST4$dis)
  MST4$b <- 1
  for (i in 1:nrow(MST4)) {
    if((MST4$V1[i]<=n1 & MST4$V2[i]<=n1)|(MST4$V1[i]>n1 &  MST4$V2[i]>n1)){
      MST4$b[i] <- 0
    }
  }
  return(sum(MST4$order*MST4$b))
} 
times=5000

p_value_f<-function(k,sta){
  rank1=rep(0,k)
  p_value_rank=rep(0,k)
  for (j in 1:k) {
    rank1[j]<-T_sta_rank(j,am_array,am_con_array)
  }
  
  for (i in 1:k) {
    p_value_rank[i]<-length(which(sta[,i]<=rank1[i]))/5000
  }
  return(p_value_rank)
}
#sta<-as.matrix(sta)
#sta<-sta[,-1]
pvalue_park<-p_value_f(10,sta)


pvalue<-matrix(0,10,4)
pvalue[,1:3]<-p_three
pvalue[,4]<-pvalue_park*0.000001
colnames(pvalue)<-c("1979","general_2017","weighted_2018","cross_rank")
plot_pvalue<-melt(pvalue)
plot_pvalue$k<-plot_pvalue$Var1
ggplot(plot_pvalue, aes(x=k, y=value)) + geom_line(aes(color=Var2))+geom_point(size=3,aes(shape=Var2,colour=Var2))+labs(y="p-value",fill="shapes")+ggtitle("p-value of Pd")+ theme(plot.title = element_text(hjust = 0.5,size=23))+theme(legend.position="none")+theme(legend.title =element_blank())+theme(legend.key.size = unit(23, "pt"))+ theme(legend.text=element_text(size=20))+ theme(axis.title=element_text(size=20),axis.title.x =element_text(size=20), axis.title.y=element_text(size=20),axis.text.y=element_text(size=20),axis.text.x=element_text(size=20))

pvalue_flight<-pvalue

pvalue_flight[,1]<-pvalue_flight[,1]*0.00001


colnames(pvalue_flight)<-c("1979","general_2017","weighted_2018","cross_rank")
plot_pvalue_flight<-melt(pvalue_flight)
plot_pvalue_flight$k<-plot_pvalue_flight$Var1
ggplot(plot_pvalue_flight, aes(x=k, y=value)) + geom_line(aes(color=Var2))+geom_point(size=3,aes(shape=Var2,colour=Var2))+labs(y="p-value",fill="shapes")+ggtitle("p-value of flight")+ theme(plot.title = element_text(hjust = 0.5,size=23))+theme(legend.position="none")+theme(legend.title =element_blank())+theme(legend.key.size = unit(23, "pt"))+ theme(legend.text=element_text(size=20))+ theme(axis.title=element_text(size=20),axis.title.x =element_text(size=20), axis.title.y=element_text(size=20),axis.text.y=element_text(size=20),axis.text.x=element_text(size=20))

