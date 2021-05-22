#function of our test statistic
T_sta_cross<-function(k,data1,data2){
  #data<-c(data1,data2)
  data<-rbind(data1,data2)
  n1<-nrow(data1)
  n<-nrow(data)
  #n1<-length(data1)
  #n<-length(data)
  mst_count1<-matrix(0,n,2)
  mst_count1[1:n1,1]<-1
  mst_count1[(n1+1):n,2]<-1
  data<-as.data.frame(data)
  mydist<-distance(data)
  MST1<-getGraph(mst_count1, mydist, k, graph.type = "mstree")
  MST2<-as.numeric(MST1)
  MST3<-matrix(0,k*(n-1),2)
  rule<-as.factor(rep(rep(c(1,2),each=n-1),k))
  group<-split(MST2,rule)
  MST3[,1]<-unlist(group[1])
  MST3[,2]<-unlist(group[2])
  #for (i in 1:k) {
  #MST3[(k-1)*(n-1)+1:k*(n-1),1]<-MST2[(k-1)*(n-1)+1:k*(n-1)]
  #MST3[(k-1)*(n-1)+1:k*(n-1),2]<-MST2[k*n:(k*(n-1))]
  #}
  MST4<-as.data.frame(MST3)
  MST4$dis<-0
  for (i in 1:nrow(MST4)) {
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

#find the distribution under H0
size_f<-function(k){
times <-5000
sta <- matrix(0,times,k)
Sigma<-diag(2)
for (i in 1:times) {
  data1<-mvrnorm(81,rep(0,2),Sigma)
  data2<-mvrnorm(36, rep(0,2), Sigma)
  for (j in 1:k) {
  sta[i,j] <-T_sta_cross(j,data1,data2)
}
}
return(sta)
}

sta<-size_f(3)
cv_rank<-matrix(0,2,k)
cv_rank_95<-matrix(0,2,k)
for (j in 1:k) {
cv_rank[,j]<-quantile(sta[,j],c(0.025,0.975))
}

##size and power
power_f<-function(k){
  sta<-size_f(k)
  cv_rank<-matrix(0,2,k)
  cv_rank_95<-matrix(0,2,k)
  for (j in 1:k) {
    cv_rank[,j]<-quantile(sta[,j],c(0.025,0.975))
  }
times0=500
size<-rep(0,k)
power<-rep(0,k)
rank_test_size<-matrix(0,times0,k)
rank_test_power<-matrix(0,times0,k)
for (i in 1:times0) {
  data1<-mvrnorm(81, rep(0, 2), Sigma)
  data2<-mvrnorm(36, rep(0, 2), Sigma)
  data3<-mvrnorm(81, rep(0, 2), Sigma)
  data4<-mvrnorm(36, rep(1, 2), Sigma)
  for (j in 1:k) {
  rank_test_size[i,j]<-T_sta_cross(j,data1,data2)
  rank_test_power[i,j]<-T_sta_cross(j,data3,data4)
  }
}
for (l in 1:k) {
size[l]<-length(which(rank_test_size[,l]<cv_rank[1,l]|rank_test_size[,l]>cv_rank[2,l]))/times0
power[l]<-length(which(rank_test_power[,l]<cv_rank[1,l]|rank_test_power[,l]>cv_rank[2,l]))/times0
}
return(list(size,power,sta))
}
power_size_type1<-power_f(k)
power_size_type2<-power_f(k)
power_size_type3<-power_f(k)
power_size_type4<-power_f(k)


# function of our test statistic for network data
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
  MST1<-getGraph(mst_count1, mydist,k, graph.type = "mstree")
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
#p-value of the real data
p_value_f<-function(k,sta){
  rank1=rep(0,k)
  p_value_rank=rep(0,k)
  for (j in 1:k) {
    rank1[j]<-T_sta_rank(j,am_array,am_con_array)
   # rank1[j]<-T_sta_rank(k,am_array,am_con_array)
    p_value_rank[j]<-length(which(sta[,j]<=rank1[j]))/times
  }
  return(p_value_rank)
}

pvalue_park<-p_value_f(10,sta)

