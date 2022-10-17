## Compare SPREE versions

compare_spree<- function(spree_type1, spree_type2, spree_type3=NULL,spree_type4=NULL){
  # compare_spree<- function(spree_type1, spree_type2, spree_type3=NULL,spree_type4=NULL,spree_type5=NULL){

  point0<- melt(spree_type1$updated_point)
  point0$type<- rep(deparse(substitute(spree_type1)), nrow(point0))
  point1<- melt(spree_type2$updated_point)
  point1$type<- rep(deparse(substitute(spree_type2)), nrow(point1))
  point2<- melt(spree_type3$updated_point)
  point2$type<- rep(deparse(substitute(spree_type3)), nrow(point2))
  point3<- melt(spree_type4$updated_point)
  point3$type<- rep(deparse(substitute(spree_type4)), nrow(point3))
  #  point4<- melt(spree_type5$updated_point)
  #  point4$type<- rep(deparse(substitute(spree_type5)), nrow(point4))

   point<-rbind(point0, point1, point2,point3)



  #MSE

  MSE0<- melt(spree_type1$MSE)
  MSE0$type<- rep(deparse(substitute(spree_type1)), nrow(MSE0))

  MSE1<- melt(spree_type2$MSE)
  MSE1$type<- rep(deparse(substitute(spree_type2)), nrow(MSE1))

  MSE2<- melt(spree_type3$MSE)
  MSE2$type<- rep(deparse(substitute(spree_type3)), nrow(MSE2))

  MSE3<- melt(spree_type4$MSE)
  MSE3$type<- rep(deparse(substitute(spree_type4)), nrow(MSE3))

  #  MSE4<- melt(spree_type5$MSE)
  #  MSE4$type<- rep(deparse(substitute(spree_type5)), nrow(MSE4))


  MSE<-rbind(MSE0, MSE1, MSE2, MSE3)



  #CV


  CV0<- melt(spree_type1$CV)
  CV0$type<- rep(deparse(substitute(spree_type1)), nrow(CV0))
  CV1<- melt(spree_type2$CV)
  CV1$type<- rep(deparse(substitute(spree_type2)), nrow(CV1))
  CV2<- melt(spree_type3$CV)
  CV2$type<- rep(deparse(substitute(spree_type3)), nrow(CV2))
  CV3<- melt(spree_type4$CV)
  CV3$type<- rep(deparse(substitute(spree_type4)), nrow(CV3))
  #  CV4<- melt(spree_type5$CV)
  # CV4$type<- rep(deparse(substitute(spree_type5)), nrow(CV4))

  CV<-rbind(CV0, CV1, CV2, CV3)

  points<-ggplot(point, aes(x=variable, y=value, fill=type)) +
    geom_boxplot(alpha=0.6)+labs(title="Point estimates",
                                 x ="Categories", y = "Counts")+labs(fill = "Type")+
    scale_fill_manual(values=c("#31688e", "#440154", "#35b779", "#fde725","#5ec962"))+
    theme(axis.text.x = element_text(size = 15))+
    theme(axis.text.y = element_text(size = 15))+
    theme(axis.title = element_text(size = 15)) +
    theme(plot.title = element_text(size = 17))+
    theme(
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13)
    )


  MSEs<-ggplot(MSE, aes(x=variable, y=value, fill=type)) +
    geom_boxplot(alpha=0.6)+labs(title="Mean Squared Error",
                                 x ="Categories", y = "MSE")+labs(fill = "Type")+
    scale_fill_manual(values=c("#31688e", "#440154", "#35b779", "#fde725","#5ec962"))+
    theme(axis.text.x = element_text(size = 15))+
    theme(axis.text.y = element_text(size = 15))+
    theme(axis.title = element_text(size = 15)) +
    theme(plot.title = element_text(size = 17))+
    theme(
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13)
    )


  CVs<-ggplot(CV, aes(x=variable, y=value, fill=type)) +
    geom_boxplot(alpha=0.6)+labs(title="Coefficient of variation",
                                 x ="Categories", y = "CV %")+labs(fill = "Type")+
    scale_fill_manual(values=c("#31688e", "#440154", "#35b779", "#fde725","#5ec962"))+
    theme(axis.text.x = element_text(size = 15))+
    theme(axis.text.y = element_text(size = 15))+
    theme(axis.title = element_text(size = 15)) +
    theme(plot.title = element_text(size = 17))+
    theme(
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13)
    )


  return(list(points= points,
              MSEs= MSEs,
              CVs=CVs))
}
