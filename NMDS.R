library(vegan)
library(ggplot2)
library(extrafont)
#font_import() # only one time required when first time use the library extrafont
#y
fonts() 
loadfonts()

#set meta data
MyMeta = data.frame(
  sites = c(2,3,4,5,6,7,8,9,10),
  amt = c("YL", "YL", "YL", "ML", "ML", "ML", "SL", "SL", "SL"),
  row.names = "sites")


#upload file
x<- read.delim("RPKM_original",header=TRUE,row.names=1)
rawdata <- x[,c(4,5,6,7,8,9,10,11,12)]
head(rawdata)



# Remove all gene which has 0 value in all sample
all <- apply(rawdata, 1, function(x) all(x==0) )
newdata <- rawdata[!all,]
dim(newdata)


# remove uninformative genes keep only genes that are expressed in at least 1 count in 2 samples
counts <- newdata[rowSums(newdata > 1) >= 2,]
head (counts)
dim(counts)

#transpose data
t_x <- t(counts)


#log2 conversion
logx <- log2(t_x+3)

#NMDS analysis
sol <- metaMDS(logx, distance='bray', k=2, trymax=100)
head(sol)

#convert NMDS analysis results into dataframe
NMDS = data.frame(MDS1 = sol$points[,1], MDS2 = sol$points[,2])


png("NMDS.png", units="in", family="Times New Roman",  width=2, height=2, res=300, pointsize = 1) #pointsize is font size| increase image size to see the key

#plot the dataframe
ggplot(data = NMDS, aes(MDS1, MDS2)) + 
geom_point(aes(color = MyMeta$amt),size=1)+
labs(color = "Tissue ")+ #change legend title
theme_bw()+ 
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=5, face = "bold"), # remove x-axis labels
        axis.title.y = element_text(size=5, face = "bold"), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())+
theme(legend.direction = 'horizontal', 
        legend.position = 'top',
        legend.key = element_rect(size = 3),
        legend.key.size = unit(1.5, 'lines'))+
theme(legend.text=element_text(size=5))+
# guides(color = guide_legend(nrow = 1))+
theme(plot.title = element_text(size = 5, face = "bold") , legend.title=element_text(size=6) , legend.text=element_text(size=5))
dev.off()












