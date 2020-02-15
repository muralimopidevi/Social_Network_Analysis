
#--------------------------------------INSTALLING REQURIED LIBRARIES--------------------------------------------------------
install.packages('igraph')
install.packages('rgl')
install.packages('plot3D')
install.packages('plot3Drgl')
install.packages('plyr')
install.packages('factoextra')
install.packages('NbClust')



#--------------------------------------LOADING INSTALLED LIBRARIES----------------------------------------------------------
library(rgl)
library(igraph)
library(plot3D)
library(plot3Drgl)
library (plyr)
library(factoextra)
library(NbClust)


#--------------------------------------CREATING A FOR LOOP for GENETATING RANDOM NETWORKS-----------------------------------
graph_list<- list()
random_graphs <-list()
for(i in 1:1000){
    
    s1 <- sample(100:5000,1) #randomly selecting RANGE S1 from 100 to 5000
    s2 <- sample(150:5000,1) #Randomly selecting Range S2 from 150 to 5000
    
    #creating a Random networks using erdos.renyi.game
    g <- erdos.renyi.game(s1, s2, type = "gnm")
    
    #Calculating degree_centrality
    degree_centrality <- degree(g, mode="all", normalized=T)
    
    #Calculating  closeness_centrality
    closeness_cent <- closeness(g, mode="all", weights=NA, normalized=T)
    
    #Calculating  betweenness_centrality
    betweenness_cent <- betweenness(g, directed = F, weights = NA, normalized = T)
    
    #for Storing random graphs in the form of list
    name <- paste('',i,sep='')
    tmp <- list(g)
    random_graphs[[name]] <- tmp

    
    #Co-relation between the centralities to get X, Y, Z three point from whole network.
    cor_1 <- cor(degree_centrality,closeness_cent, method = c("pearson", "kendall", "spearman"))
    cor_2 <- cor(degree_centrality,betweenness_cent, method = c("pearson", "kendall", "spearman"))
    cor_3 <- cor(closeness_cent,betweenness_cent, method = c("pearson", "kendall", "spearman"))
    
    #for Storing x, y, z points from each  random graphs.
    name1 <- paste('',i,sep='')
    tmp <- list(X=cor_1, y=cor_2, Z=cor_3)
    graph_list[[name1]] <- tmp
}


#--------------------------------------CREATING A DATAFRAME USING THE LIST OF x, Y, Z POINTS FROM RANDOM GRAPHS-------------

#viewing all RANDOM NETWOKRKS WHICH IS IN THE FORM OF LISTS.

#---->>>>>>>>>>>> un-comment 'random_graphs' in below line to view all randomly generated networks .

# random_graphs       # <<<---------------------------------------------------

#CREATING A DATAFRAME FROM LIST'S
data_frame <- ldply(graph_list, data.frame)

#removing unwanted data from DATA_FRAME
df = subset(data_frame, select = -c(.id))

#INFORLATION ABOUT DATA_FRAME WHICH IS CREATED FROM LIST's X, Y, Z.
str(df)
dim(df)
summary(df)
names(df)

#plotting from dataframe in 3D
plot3d(df$X,df$y,df$Z, col = c('red','green','blue'), size = 6)

#----------------- standardisation, where the dimensions now have a mean of zero---------------------------------
#mean of X
m_x <- mean(df$X)

#mean of Y
m_y <- mean(df$y)

#mean of Z
m_z <- mean(df$Z)

x1 <- df$X - m_x
x1

#summary of X
summary(x1)

y1 <- df$y - m_y
y1

#summary of Y
summary(y1)

z1 <- df$Z - m_y
z1

#summary of z
summary(z1)

#plotting from dataframe
plot3d(x1,y1,z1, col = c('red','green','blue'), size = 6)

# to calculate the covariance matrix. Covariance measures how dimensions vary with respect to each other and the 
#  covariance matrix contains all covariance measures between all dimensions

cv_xx <- cov(x1, x1)
cv_xy <- cov(x1, y1)
cv_xz <- cov(x1, z1)
cv_yy <- cov(y1, y1)
cv_yx <- cov(y1, x1)
cv_yz <- cov(y1, z1)
cv_zz <- cov(z1, z1)
cv_zx <- cov(z1, x1)
cv_zy <- cov(z1, y1)


m <- matrix(c(cv_xx, cv_xy, cv_xz,cv_yx,cv_yy,cv_yz,cv_zx,cv_zy,cv_zz),
            nrow=3,
            ncol=3,
            byrow=TRUE,
            dimnames=list(c("x","y","z"),c("x","y","z")))

m

# we need to find the eigenvector and eigenvalues of the covariance matrix. An eigenvector 
# is a direction and an eigenvalue is a number that indicates how much variance is in the data in that direction

ev <- eigen(m)
ev

#The largest eigenvalue is the first principal component; we multiply the standardised values to the 
# first eigenvector, which is stored in ev$vectors[,1].

PC1 <- x1 * ev$vectors[1,1] + y1 * ev$vectors[2,1] + z1 * ev$vectors[3,1]
PC1


PC2 <- x1 * ev$vectors[1,2] + y1 * ev$vectors[2,2] +  z1 * ev$vectors[3,2]
PC2

PC3 <- x1 * ev$vectors[1,3] + y1 * ev$vectors[2,3] +  z1 * ev$vectors[3,3]
PC3

#creating PCA DATA FRAME
PCA_data_frame <- data.frame(PC1,PC2,PC3)

PCA_data_frame


#plotting PCA
plot3d(PC1,PC2,PC3, col = c('red','green','blue'), size = 6)


#Now to perform PCA using the prcomp() function.
pca <- prcomp(df)
pca

#names in  PCA
names(pca)

pca$x

#plotting PCA
plot3d(pca$x[,1], pca$x[,2], pca$x[,3], col = c('red','green','blue'), size = 6)

#SCREE PLOT
fviz_eig(pca, addlabels = TRUE, ylim = c(0, 50))

#determining the optimal number of clusters
fviz_nbclust(df, kmeans, method = "gap_stat")



#Kmeans ALGORITHM WITH PCA VALUES like PC1, PC2, PC3
set.seed(123)
k <- kmeans(df,3, nstart=25, iter.max=1000)

new = cbind(df,cluster = k$cluster)
with(new,plot3d(PC1,PC2,PC3, col=k$cluster, size=1, type='s'))

#Kmeans ALGORITHM WITHOUT PCA VALUES
set.seed(123)
k <- kmeans(df,3, nstart=25, iter.max=1000)


#plotting from Clusters using K-means
plot3d(df$X,df$y,df$Z, col = k$cluster, size=7)


#K-Means
set.seed(123)
km.res <- kmeans(df,9, nstart = 25)

fviz_cluster(km.res, data = df,
             ellipse.type = "convex",
             palette = "jco",
             ggtheme = theme_minimal())



# Compute hierarchical clustering
res.hc <- df %>%
    scale() %>%                    # Scale the data
    dist(method = "euclidean") %>% # Compute dissimilarity matrix
    hclust(method = "ward.D2")     # Compute hierachical clustering

# Visualize using factoextra
# Cut in 4 groups and color by groups
fviz_dend(res.hc, k = 3, # Cut in four groups
          cex = 0.5, # label size
          k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE # Add rectangle around groups
)


#Determining the optimal number of clusters
# Compute
set.seed(123)
res.nbclust <- df %>%
    scale() %>%
    NbClust(distance = "euclidean",
            min.nc = 2, max.nc = 10, 
            method = "complete", index ="all") 
# Visualize
fviz_nbclust(res.nbclust, ggtheme = theme_minimal())




#Compute and visualize hierarchical clustering:
set.seed(123)
# Enhanced hierarchical clustering, cut in 3 groups
res.hc <- df %>%
    scale() %>%
    eclust("hclust", k = 3, graph = FALSE)

# Visualize with factoextra
fviz_dend(res.hc, palette = "jco",
          rect = TRUE, show_labels = FALSE)


#Inspect the silhouette plot:
fviz_silhouette(res.hc)



