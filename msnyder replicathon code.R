install.packages("dplyr")
install.packages("tidyr")
install.packages("swirl")

#determining the number of cell lines there are in the dataset

cell_line<- summarizedPharmacoData$cellLine

counter<-0
newvec<-vector()
for (i in cell_line){
  if (i %in% newvec==FALSE){
    newvec<-c(newvec, i)
    counter=counter+1
  }
}
print(counter)
