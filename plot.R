rm(ls=list())
library(dplyr)
library(tidyr)
library(ggplot2)
library(xlsx)
library(grid)
setwd("/home/o0/Desktop/pom/")
# dat  = read.xlsx(file = "result_pom.xlsx",sheetIndex = 1)
# 
# colnames(dat) = c("beta","Tp","Tm","Tw")
# dat$beta = factor(round(dat$beta,2))
# dat$id = 1:nrow(dat)
# 
# df1= dat %>% gather("item",value,-c("id","beta")) %>% 
#   bind_cols(data.frame(item_id=rep(1:3,each=nrow(dat))))
# 
# df1$item = as.factor(df1$item)
m1=plothi(maf=0.05)
m2=plothi(maf=0.3)
m3=plothi(maf=0.35)
grid.newpage()
pushViewport(viewport(layout = grid.layout(3,3)))
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
print(m1, vp = vplayout(1,1))   ###将（1,1)和(1,2)的位置画图c
print(m2, vp = vplayout(1,2))###将(2,1)的位置画图b
print(m3,vp = vplayout(1,3))

# m=list()
# for(i in 1:9){
#   m[[i]] = plothi(maf=0.05*i)
#   print(m[[i]],vp=vplayout(ceiling(i/3),i%%3))
# }



plothi<-function(maf){
  dat  = read.xlsx(file = "result_pom.xlsx",sheetIndex =maf*20)
  
  colnames(dat) = c("beta","Tp","Tm","Tw")
  dat$beta = factor(round(dat$beta,2))
  dat$id = 1:nrow(dat)
  
  df1= dat %>% gather("item",value,-c("id","beta")) %>% 
    bind_cols(data.frame(item_id=rep(1:3,each=nrow(dat))))
  df1$item = as.factor(df1$item)
  if(maf==0.05){
    a<-ggplot(df1,aes(x=beta, y=value, fill=factor(item,levels = levels(df1$item)[c(2,1,3)])))+
      geom_bar(stat='identity',
               position=position_dodge())+ #本来不用reverse是从上到下，反过来
      # scale_fill_discrete(name='legend', #图例项（或者用scale_fill_discrete)
      #                     labels=levels(df1$item)[c(2,1,3)])+
      guides(fill=guide_legend(title=NULL))+
      scale_fill_manual(breaks = c("Tp", "Tm","Tw"),
                        values=alpha(c("seagreen4", "dodgerblue2", "red2"), .8))+
      theme(legend.position="right")+#图例位置
      theme(legend.text = element_text(size = 6))+
      theme(legend.key.size=unit(0.3,'cm'))+
      theme(legend.key.width=unit(0.35,'cm'))+      # ggtitle(paste('MAF = ',maf))+
      labs(title=paste('MAF = ',maf))+
      theme(plot.title = element_text(hjust = 0.5))+
      theme(axis.title.x =element_text(size=8), axis.title.y=element_text(size=8),
            ,title =element_text (size=10))+
      scale_x_discrete(name=expression(beta), #x轴坐标名称
                       labels= c("ln1.1","ln1.2","ln1.3","ln1.4","ln1.5"))+ #离散的标签
      scale_y_continuous(name=paste("R",'(%)')) #y轴坐标名称
  }else{
    a<-ggplot(df1,aes(x=beta, y=value, fill=factor(item,levels = levels(df1$item)[c(2,1,3)])))+
      geom_bar(stat='identity',
               position=position_dodge())+ #本来不用reverse是从上到下，反过来
      # scale_fill_discrete(name='legend', #图例项（或者用scale_fill_discrete)
      #                     labels=levels(df1$item)[c(2,1,3)])+
      guides(fill=guide_legend(title=NULL))+
      scale_fill_manual(breaks = c("Tp", "Tm","Tw"),
                        values=alpha(c("seagreen4", "dodgerblue2", "red2"), .8))+
      theme(legend.position="right")+#图例位置
      theme(legend.position = "none") +
      theme(legend.text = element_text(size = 6))+
      ggtitle(paste('MAF = ',maf))+
      theme(plot.title = element_text(hjust = 0.5))+
      theme(axis.title.x =element_text(size=8), axis.title.y=element_text(size=8)
            ,title =element_text (size=10))+
      scale_x_discrete(name=expression(beta), #x轴坐标名称
                       labels= c("ln1.1","ln1.2","ln1.3","ln1.4","ln1.5"))+ #离散的标签
      scale_y_continuous(name=paste("R",'(%)')) #y轴坐标名称
  }


  return(a)
}
