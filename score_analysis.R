#!/usr/bin/env Rscript
#import arguments 
args = commandArgs(trailingOnly=TRUE)

#planed usage
#Rscript scores_analysis.R -Print graphs(TRUE or FALSE) -Print recommended Parameters(TRUE or FALSE) -Location of scores file to import/wd -Project Name 
#Rscript ~/Scripts/github/score_analysis.R TRUE TRUE ~/Scripts/github/Mid1_run Mid1_run


printgraphsT_F<-args[1]
printRecomPara<-args[2]
location<- args[3]
projectname<-args[4]

#set libraries
library(ggplot2)
#Import data
#dataset <- read.csv("~/Documents/OneDrive/Antarctica Files/LS Project/runAssembly_opti/bunassemnly_opti/Mid_3_newbler_scores.txt")

dataset<-read.csv(paste(location,projectname,"_newbler_scores.txt",sep = ""))

#set wd
#setwd("~/Documents/OneDrive/Antarctica Files/LS Project/runAssembly_opti/bunassemnly_opti")

setwd(location)

#MAKE GRAPHS for each parameter and each score
make_graphs<-function(dataset,size,projectname){
  ggplot(data=dataset, aes(x=dataset$Readlength,y=dataset[,5],group=dataset$Readlength))+
    geom_boxplot()+
    geom_smooth(method = "glm", se=FALSE, color="red", aes(group=1))+
    labs(title =paste("Readlength vs", names(dataset)[5]), 
         y= paste(names(dataset)[5],"(# of contigs)"), 
         x= "Readlength (# of bases)")+
    theme_bw()
  filename_x<-paste(projectname,"_Readlength","_",names(dataset)[5],".png")
  filename_x<-gsub("\\s", "",filename_x)
  ggsave(filename=filename_x, width = size, height = size)
  
  ggplot(data=dataset, aes(x=dataset$Overlap,y=dataset[,5],group=dataset$Overlap))+
    geom_boxplot()+
    geom_smooth(method = "glm", se=FALSE, color="red", aes(group=1))+
    labs(title =paste("Overlap vs", names(dataset)[5]), 
         y= paste(names(dataset)[5],"(# of contigs)"), 
         x= "Overlap (#  of bases)")+
    theme_bw()
  filename_x<-paste(projectname,"_Overlap","_",names(dataset)[5],".png")
  filename_x<-gsub("\\s", "",filename_x)
  ggsave(filename=filename_x, width = size, height = size)
  
  ggplot(data=dataset, aes(x=dataset$Id,y=dataset[,5],group=dataset$Id))+
    geom_boxplot()+
    geom_smooth(method = "glm", se=FALSE, color="red", aes(group=1))+
    labs(title =paste("Id% vs", names(dataset)[5]), 
         y= paste(names(dataset)[5],"(# of contigs)"), 
         x= "Id %")+
    theme_bw()
  filename_x<-paste(projectname,"_Id","_",names(dataset)[5],".png")
  filename_x<-gsub("\\s", "",filename_x)
  ggsave(filename=filename_x, width = size, height = size)
  
for(i in 6:9){
ggplot(data=dataset, aes(x=dataset$Readlength,y=dataset[,i],group=dataset$Readlength))+
  geom_boxplot()+
    geom_smooth(method = "glm", se=FALSE, color="red", aes(group=1))+
    labs(title =paste("Readlength vs", names(dataset[i])), 
         y= paste(names(dataset[i]),"(# of bases)"), 
         x= "Readlength (# of bases)")+
  theme_bw()
  filename_x<-paste(projectname,"_Readlength","_",names(dataset[i]),".png")
  filename_x<-gsub("\\s", "",filename_x)
  ggsave(filename=filename_x, width = size, height = size)

ggplot(data=dataset, aes(x=dataset$Overlap,y=dataset[,i],group=dataset$Overlap))+
  geom_boxplot()+
  geom_smooth(method = "glm", se=FALSE, color="red", aes(group=1))+
  labs(title =paste("Overlap vs", names(dataset[i])), 
       y= paste(names(dataset[i]),"(# of bases)"), 
       x= "Overlap (#  of bases)")+
  theme_bw()
  filename_x<-paste(projectname,"_Overlap","_",names(dataset[i]),".png")
  filename_x<-gsub("\\s", "",filename_x)
  ggsave(filename=filename_x, width = size, height = size)

ggplot(data=dataset, aes(x=dataset$Id,y=dataset[,i],group=dataset$Id))+
  geom_boxplot()+
  geom_smooth(method = "glm", se=FALSE, color="red", aes(group=1))+
  labs(title =paste("Id % vs", names(dataset[i])), 
       y= paste(names(dataset[i]),"(# of bases)"), 
       x= "Id %")+
  theme_bw()
  filename_x<-paste(projectname,"_Id","_",names(dataset[i]),".png")
  filename_x<-gsub("\\s", "",filename_x)
  ggsave(filename=filename_x, width = size, height = size)
}
}

#Test variables 
#printgraphsT_F=TRUE
#projectname<-"Mid3"
#End of test variables

if (printgraphsT_F==TRUE){
  print("Building Graphs")
  make_graphs(dataset,5,projectname)
  print("Done Building Graphs")
}

if(printRecomPara==TRUE){
print("Printing Analysis Report")
#Summary Stats
#stats_summary<-summary(dataset[,5:9])

min_numberOfContigs_parameters<-which(dataset$numberOfContigs==min(dataset$numberOfContigs,na.rm=TRUE))
max_numberOfBases_parameters<-which(dataset$numberOfBases==max(dataset$numberOfBases,na.rm=TRUE))
max_avgContigSize_parameters<-which(dataset$avgContigSize==max(dataset$avgContigSize,na.rm=TRUE))
max_N50ContigSize_parameters<-which(dataset$N50ContigSize==max(dataset$N50ContigSize,na.rm=TRUE))
max_largestContigSize_parameters<-which(dataset$largestContigSize==max(dataset$largestContigSize,na.rm=TRUE))

df1<-data.frame(dataset[min_numberOfContigs_parameters,1:4])
df2<-data.frame(dataset[max_numberOfBases_parameters,1:4])
df3<-data.frame(dataset[max_N50ContigSize_parameters,1:4])
df4<-data.frame(dataset[max_largestContigSize_parameters,1:4])
df5<-data.frame(dataset[max_avgContigSize_parameters,1:4])

#Make text output
line1<-paste("Based on lowest number of contigs -->",min(dataset$numberOfContigs,na.rm=TRUE),"(# of contigs)\n", "- The best assembly was:",df1[1,1],"\n",
             "- That assembly was made using the following parameters:","\n", " Read length=",df1[1,2],
             "\n"," Overlap length=",df1[1,3],"\n"," Percent Identity=",df1[1,4])

line2<-paste("Based on largest number of total bases -->",max(dataset$numberOfBases,na.rm=TRUE),"(# of bases)\n", "- The best assembly was:",df2[1,1],"\n",
             "- That assembly was made using the following parameters:","\n", " Read length=",df2[1,2],
             "\n"," Overlap length=",df2[1,3],"\n"," Percent Identity=",df2[1,4])

line3<-paste("Based on largest N50ContigSize -->",max(dataset$N50ContigSize,na.rm=TRUE),"(# of bases)\n", "- The best assembly was:",df3[1,1],"\n",
             "- That assembly was made using the following parameters:","\n", " Read length=",df3[1,2],
             "\n"," Overlap length=",df3[1,3],"\n"," Percent Identity=",df3[1,4])

line4<-paste("Based on the size of the largestContigSize -->",max(dataset$largestContigSize,na.rm=TRUE),"(# of bases)\n","- The best assembly was:",df4[1,1],"\n",
             "- That assembly was made using the following parameters:","\n", " Read length=",df4[1,2],
             "\n"," Overlap length=",df4[1,3],"\n"," Percent Identity=",df4[1,4])

line8<-paste("Based on highest average contig size -->",max(dataset$avgContigSize,na.rm=TRUE),"(# of bases)\n", "- The best assembly was:",df5[1,1],"\n",
             "- That assembly was made using the following parameters:","\n", " Read length=",df5[1,2],
             "\n"," Overlap length=",df5[1,3],"\n"," Percent Identity=",df5[1,4])

line5<-paste("# Based on ",projectname,"_newbler_scores.txt", sep="")
line6<-paste("####################################################################################################")
line7<-paste("# Newbler Parameter Sweep Analysis Report")

lines_to_write<-c(line6,line7,line5,line6,line1,line2,line8,line3,line4)
fileConn<-file(paste(projectname,"_newbler_scores_Analysis_Report.txt",sep = ""))
writeLines(lines_to_write, fileConn)
close(fileConn)
print("Done Printing Analysis Report")
}
