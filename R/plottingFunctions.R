#' ### Dependencies:
   
  library("igraph")
  library("sf")
  library("geosphere")
  library("dplyr")

 
#' Helper function for use within the larger Dijkstra algorithm. 
#' Takes in two lat, lon pairs and returns distance (curvature accounted) in kilometers using geoSphere's distGeo function. 
#' 
#' 
#' @param Lat1: Latitude of first point.
#' @param Lon1: Longitude of first point.
#' @param Lat2: Latitude of second point.
#' @param Lon2: Longitude of second point. 
#' 
#' @return: Geodesic distance (accounting for curvature) between the two points. 
#' 
  calcdist <- function(lat1,lon1,lat2,lon2){
    return (geosphere::distGeo(c(lat1,lon1),c(lat2,lon2))/1000)
  }
#' 
#' Helper function for merging together columns. 
#' 
#' @param muni: municipality in question. 
#' 
#' @return: Municipality with CO pasted in front. 
#' # Adds "CO" at the beginning of municipality number to facilitate merge. 
  addCO <- function(muni){
    ifelse(nchar(paste(muni)) == 4,paste("CO0",muni,sep = ""),paste("CO",muni,sep = ""))
  }
#' 
#' Helper function to run Dijkstra given parameters and create corresponding distance dataframe. 
#' 
#' @param metric Distance metric being used. Metric 0 corresponds to hops in the graph, metric 1 to crow-flies distance, metric 2 to PCA without roads, and metric 3 to PCA with roads. 
#' @param a Parameter for first principal component.
#' @param b Parameter for second principal component. 
#' @param epicenter_1 First epicenter for violence origin (municipality code), allowing us to create models of outward spread from different origin pairs. 
#' @param epicenter_2 Second epicenter for violence origin (municipality code). 
#'
#' @return Dataframe containing distance, ring and municipality code. 
#' Dataframe containing Municipality code, Year, Share, Ring, and Standard Deviation. 
#'
#'\itemize{
#'    \item \code{total_dist} Distance between given pair of adjacent municipalities under metric. 
#'  }
generate_distances <- function(metric,a,b,epicenter_1,epicenter_2){
  


munigraph <- forcedMigration::munigraph

# Figure out which PCA version to use (roads or no roads)

cross_section_merged <- forcedMigration::cross_section_merged
  
# Figure out which PCA version to use (roads or no roads)
for_pca <- forcedMigration::pca
if(metric == 3){
  for_pca <- forcedMigration::roads_pca
}

# Shift PCA's to positive range. 
for_pca$PC1 <- for_pca$PC1-min(for_pca$PC1)
for_pca$PC2 <- for_pca$PC2-min(for_pca$PC2)
for_pca$PC1 <- for_pca$PC1 / (max(for_pca$PC1))
for_pca$PC2 <- for_pca$PC2 / (max(for_pca$PC2))

for(i in 1:1120){
  for(j in 1:1120){
    # Initially, set total_dist (distance between municipalities i and j) to an effectively infinite number. 
    total_dist<-1000000
    
    # Check that the municipalities are adjacent. 
    if(igraph::are.connected(munigraph,i,j)){
      
      # Calculate crow-flies-distance. 
      d1 <- forcedMigration::calcdist(forcedMigration::AttributeTableFinal$latnum[i],forcedMigration::AttributeTableFinal$lonnum[i],forcedMigration::AttributeTableFinal$latnum[j],forcedMigration::AttributeTableFinal$lonnum[j])
      
      # If we have no geographic covariates for at least one municipality in the pair (and therefore did not calculate a PCA distance for this pair), set total_dist to 1. (This is the "pure graph hops" case). 
      total_dist <- 1
      
      # If using the crow-flies distance metric (or hiking, for cases where we lack elevation covariate), then set total_dist to d1. 
      if(metric == 1 | metric == 4){
        total_dist <- d1
      }    
      
      # For distinct municipalities for which we have all covariates, calculate the distance. (For municipalities lacking geographic covariates, we default to crow distance). 
      
      if(metric>1 & i!= j & i<nrow(forcedMigration::cross_section_merged) & j<nrow(forcedMigration::cross_section_merged)){
        # Retrieve principal components 1 and 2 for municipalities i and j from PCA dataframe.  
        pcode_i <- forcedMigration::AttributeTableFinal$ADM2_PCODE[i]
        pcode_j <- forcedMigration::AttributeTableFinal$ADM2_PCODE[j]
        row <- for_pca[for_pca$muni_i == pcode_i & for_pca$muni_j == pcode_j,]$index
        PCA_1 <- for_pca$PC1[row]
        PCA_2 <- for_pca$PC2[row]
        
        
        # Calculate total distance using our distance formula and save for summary statistics. 
        total_dist <- exp((a*PCA_1+b*PCA_2))
        #summary_db[row] <- total_dist
        
        # If we are using the hiking metric (and both municipalities are among the 1,076 for which we have geographic covariates), calculate and update distance. Where we don't have covariates, the default is d1.  
        if(metric==4 & i<nrow(forcedMigration::cross_section_merged) & j<nrow(forcedMigration::cross_section_merged) ){
          elev_difference <- max((forcedMigration::cross_section_merged$alt_mean[i]-forcedMigration::cross_section_merged$alt_mean[j]),0)
          total_dist <- d1 + 0.6 * elev_difference
          #summary_db[row] <- total_dist
        }
      }
      
      # Update and set edge weight. 
      #print(total_dist)
      edge <- igraph::get.edge.ids(munigraph,c(i,j))
      munigraph <- igraph::set_edge_attr(munigraph,"weight",edge,total_dist)
    }
  }
}

vertex_ids = vector(length = 1120)
for(i in 1:1120){
  vertex_ids[i]<- forcedMigration::AttributeTableFinal$ADM2_PCODE[i]
}

# Calculate Dijkstra distances 
delta_1 = igraph::distances(munigraph,v = which(vertex_ids == epicenter_1),to = igraph::V(munigraph),algorithm = "dijkstra")
delta_2 = igraph::distances(munigraph,v = which(vertex_ids == epicenter_2),to = igraph::V(munigraph),algorithm = "dijkstra")
deltas <- data.frame(matrix(c(delta_1,delta_2),ncol = 2))
colnames(deltas) <- c("delta_1","delta_2");
deltas <- mutate(deltas, delta_min = as.numeric(purrr::map2(delta_1,delta_2,min)))



# Merge our distances (from specified epicenter) into dataframe containing map polygons. 
distances <- as.numeric(deltas[,3])


# Dataset for merging delta and ring_num into other datasets by ADM2_PCODE. 
merge_deltas <- data.frame(delta = distances,ADM2_PCODE = vertex_ids)

# Calculate ring as 10 * decile in distance distribution plus 1.  
merge_deltas <- mutate(merge_deltas, ring_num = as.integer(10*ecdf(merge_deltas$delta)(delta))+1)


# For the single largest delta in the dataset, decile is 10 and so ring is 11. 
# Set ring_num to 10 for this case. 
merge_deltas <- mutate(merge_deltas, ring_num = ifelse(ring_num == 11,10,ring_num))
return(merge_deltas)
}
#' Helper function to create distances given epicenter and metric parameters. 
#' @param merge_deltas Output from generate_distances; contains mapping of municipality to distance and ring under given metric, a, b, and epicenter. 
#' 
#' @return Dataframe containing municipality, ring, share of total violence for ring in given year, and standard deviation of violence across municipalities in a given ring-year. 
#' \itemize{
#' \item \code{FinalWithDeltas} Dataframe containing geographic covariates and final distances. 
#'    \item \code{vdf} Dataframe containing year-ring violence. 
#'    \item \code{sdf} Dataframe containing standard errors of year-ring violence. 
#'    \item \code{sharedf} Dataframe containing final output. 
#' }
#' 
#' 
#' 
get_df <- function(merge_deltas){
FinalWithDeltas <- merge(forcedMigration::cross_section_merged,merge_deltas,by = "ADM2_PCODE")
violence_data <- forcedMigration::violence_data[violence_data$year < 2009,]
violence_set <- merge(violence_data,FinalWithDeltas,by="ADM2_PCODE")      

#print(summary(FinalWithDeltas$delta))

# Take out municipalities with 0 violence. 

violence_subset <- filter(violence_set,year == 2008)
has_violence_list <- filter(violence_subset,cum_victims_UR >0)$ADM2_PCODE
violence_set <- filter(violence_set, ADM2_PCODE %in% has_violence_list)


# For each ring and year, stores total number of deaths in that ring in that year. 
victimsUR_matrix <- matrix(0,nrow = 10,ncol = 13)

# For each ring and year, stores standard deviation of the set {deaths in municipality in year : municipality in ring} . 
stderr_matrix <- matrix(0,nrow = 10,ncol = 13)

# Stores the two above quantities in their respective matrices. 
for(ring in 1:10){
  for(year in 1996:2008){
    victimsUR_matrix[ring,year-1995] <- sum(violence_set[violence_set$year == year & violence_set$ring == ring,]$victims__UR)
    stderr_matrix[ring,year-1995] <- std.error(violence_set[violence_set$year == year & violence_set$ring == ring,]$victims__UR)
  }
}

# Convert matrices to dataframes. 
vdf <- as.data.frame(t(victimsUR_matrix))
sdf <- as.data.frame(stderr_matrix)

# Generate share of total ring deaths that occur in a given year by dividing each ring-year death figure by total ring deaths. 
sharedf <- mutate(vdf,V1 = V1/sum(V1),V2 = V2/sum(V2),V3 = V3/sum(V3),V4 = V4/sum(V4), V5 = V5/sum(V5), V6 = V6/sum(V6), V7 = V7/sum(V7), V8 = V8/sum(V8), V9 = V9/sum(V9), V10 = V10/sum(V10))
sharedf <- as.data.frame(t(sharedf))

cnames = paste(1996:2008)
colnames(sharedf) <- cnames
colnames(sdf) <- cnames
rownames(sharedf) <- 1:10
rownames(sdf) <- 1:10

# Since the rows represent rings, explicitly create a ring variable corresponding to row number. 
sharedf <- cbind("ring" = as.numeric(rownames(sharedf)),sharedf)

# Make dataframes long to allow for plotting. 
sharedf <- pivot_longer(sharedf,cnames,names_to = "year",values_to = "share")
sdf <- cbind("ring" = as.numeric(rownames(sdf)),sdf)
sdf <- pivot_longer(sdf,cnames,names_to = "year",values_to = "sd")
sharedf <- merge(sharedf,sdf,by = c("ring","year"))

return(sharedf)

}  
#'
#' Heatmap plot to visualize distribution of violence in rings over time. 
#' 
#' @param metric Distance metric being used. Metric 0 corresponds to hops in the graph, metric 1 to crow-flies distance, metric 2 to PCA without roads, and metric 3 to PCA with roads. 
#' @param a Parameter for first principal component.
#' @param b Parameter for second principal component. 
#' @param epicenter_1 First epicenter for violence origin (municipality code), allowing us to create models of outward spread from different origin pairs. 
#' @param epicenter_2 Second epicenter for violence origin (municipality code). 
#' 
#' @return NULL (plots corresponding graph). 
#' 
#'
heatmap <- function(metric,a,b,epicenter_1,epicenter_2){
merged_deltas <- forcedMigration::generate_distances(metric,a,b,epicenter_1,epicenter_2)
sharedf <- forcedMigration::get_df(merged_deltas)
ggplot(data=sharedf,mapping=aes(x=year,y=ring,fill=share))+
  geom_tile()+theme_minimal()+scale_fill_gradient(name="Violence Share",low="darkblue",high="red")+
  scale_x_discrete(breaks = seq(1,10,3),label = paste(seq(1,10,3)))
}
#'
#' Star Wars plot (without bars) to visualize violence over time within rings.  
#' 
#' 
#' @param metric Distance metric being used. Metric 0 corresponds to hops in the graph, metric 1 to crow-flies distance, metric 2 to PCA without roads, and metric 3 to PCA with roads. 
#' @param a Parameter for first principal component.
#' @param b Parameter for second principal component. 
#' @param epicenter_1 First epicenter for violence origin (municipality code), allowing us to create models of outward spread from different origin pairs. 
#' @param epicenter_2 Second epicenter for violence origin (municipality code). 
#' 
#' @return NULL (plots corresponding graph). 
#'
starwars <-  function(metric,a,b,epicenter_1,epicenter_2){
merged_deltas <- forcedMigration::generate_distances(metric,a,b,epicenter_1,epicenter_2)
sharedf <- forcedMigration::get_df(merged_deltas)
ggplot(sharedf,aes(x = year,y = share,group = factor(ring)))+geom_point(aes(colour = factor(ring)))+facet_wrap(~ ring)
}
#'
#' Star Wars plot (with bars) to visualize violence over time within rings. 
#' 
#' @param metric Distance metric being used. Metric 0 corresponds to hops in the graph, metric 1 to crow-flies distance, metric 2 to PCA without roads, and metric 3 to PCA with roads. 
#' @param a Parameter for first principal component.
#' @param b Parameter for second principal component. 
#' @param epicenter_1 First epicenter for violence origin (municipality code), allowing us to create models of outward spread from different origin pairs. 
#' @param epicenter_2 Second epicenter for violence origin (municipality code). 
#' 
#' @return NULL (plots corresponding graph). 
#' 
#' 
starwars_bars <- function(metric,a,b,epicenter_1,epicenter_2){
merged_deltas <- forcedMigration::generate_distances(metric,a,b,epicenter_1,epicenter_2)
sharedf <- forcedMigration::get_df(merged_deltas)
ggplot(sharedf,aes(x = year,y = share,group = factor(ring)))+geom_point(aes(colour = factor(ring)))+geom_errorbar(aes(ymin = share-sd, ymax = share+sd,colour = factor(ring)))+facet_wrap(~ ring)
}

#'
#' Spatial map plot to visualize distances. 
#' 
#' 
#' @param metric Distance metric being used. Metric 0 corresponds to hops in the graph, metric 1 to crow-flies distance, metric 2 to PCA without roads, and metric 3 to PCA with roads. 
#' @param a Parameter for first principal component.
#' @param b Parameter for second principal component. 
#' @param epicenter_1 First epicenter for violence origin (municipality code), allowing us to create models of outward spread from different origin pairs. 
#' @param epicenter_2 Second epicenter for violence origin (municipality code). 
#' 
#' @return NULL (plots corresponding graph). 
#'
#'
spatial_map <- function(metric,a,b,epicenter_1,epicenter_2){
merged_deltas <- forcedMigration::generate_distances(metric,a,b,epicenter_1,epicenter_2)
muni_pol <- forcedMigration::municipalities_shapefile.shp
muni_pol <- select(muni_pol,"ADM2_PCODE")
merged_map <- merge(merged_deltas,muni_pol)
ggplot() +
  geom_sf(data = st_as_sf(merged_map), aes(fill=ring_num),color = 'grey34',lwd=.05) +
  scale_fill_stepsn(colors = c("red","gold","darkgreen","blue","violet"),values = NULL,space = "Lab",na.value = "grey50",guide = "coloursteps",aesthetics = "fill",breaks = c(1:10),limits= c(1, 10))
}
#'
#' 
#'
#' 
#'



