plot_profiles <- function(data, group){
  
  total <- length(unique(pull(argo, lovbio)))
  pages <- ceiling(total/6)
  
  pdf("Output/Test_pdf.pdf", paper = "a4")
  for(i in c(1:pages)){
    print(ggplot(argo)+
            geom_point(aes(x = chla, y = -depth, colour = "Chla"))+
            geom_point(aes(x = fluo * 2, y = - depth, colour = "Fluo"))+
            scale_color_brewer(palette = "Set1")+
            ggforce::facet_wrap_paginate(.~code * lovbio, scales = "free", ncol = 3, nrow = 2, page = i))
  }
  dev.off()
  }


plot_profiles <- function(data, group){
    ggplot(argo)+
            geom_point(aes(x = chla, y = fluo * 2))+
            scale_color_brewer(palette = "Set1")+
            facet_wrap(.~code, scales = "free")
  }
  dev.off()
}
