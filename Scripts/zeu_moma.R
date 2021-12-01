Zeu_moma<-function(chl_npq, dep_chl) {
  Zeu<-NA
  
  # interpolate Fchla per meter
  dep_int<-c(0:500)
  chl_int<-NA
  
  chl_int<-approx(dep_chl,chl_npq,dep_int)$y
  
  chl_int[which(is.na(chl_int)==T)]<-0
  
  # integrate Fchla
  chl_sum<-NA
  chl_sum<-cumsum(chl_int)
  
  
  # Calculate the Ze from Morel with iteration
  
  for (i in dep_int) {
    #print(i)
    ze<-NA
    ze<-912.5*(chl_sum[which(dep_int==i)]^(-0.839))
    
    if(length(ze)==0){
      ze<-NA}
    
    if (ze < i & (is.na(ze)==F)) {
      Zeu <- ze
      
      if(ze > 102) {
        for (j in dep_int) {
          #print(j)
          ze<-NA
          ze<-426.3*(chl_sum[which(dep_int==j)]^(-0.547))
          
          if(length(ze)==0){
            ze<-NA}
          
          if (ze < j & (is.na(ze)==F)) {
            Zeu<-ze
            break
          }
        }
      }
      
      break
    }
    
  }
  
  
  return(Zeu)
}
