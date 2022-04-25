#' ### Dependencies:



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
  
  geographic_covariates <- forcedMigration::geographic_covariates
  
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
        d1 <- forcedMigration::calcdist(forcedMigration::full_municipalities$latnum[i],forcedMigration::full_municipalities$lonnum[i],forcedMigration::full_municipalities$latnum[j],forcedMigration::full_municipalities$lonnum[j])
        
        # If we have no geographic covariates for at least one municipality in the pair (and therefore did not calculate a PCA distance for this pair), set total_dist to 1. (This is the "pure graph hops" case). 
        total_dist <- 1
        
        # If using the crow-flies distance metric (or hiking, for cases where we lack elevation covariate), then set total_dist to d1. 
        if(metric == 1 | metric == 4){
          total_dist <- d1
        }    
        
        # For distinct municipalities for which we have all covariates, calculate the distance. (For municipalities lacking geographic covariates, we default to crow distance). 
        
        if(metric>1 & i!= j & i<nrow(forcedMigration::geographic_covariates) & j<nrow(forcedMigration::geographic_covariates)){
          # Retrieve principal components 1 and 2 for municipalities i and j from PCA dataframe.  
          pcode_i <- forcedMigration::full_municipalities$ADM2_PCODE[i]
          pcode_j <- forcedMigration::full_municipalities$ADM2_PCODE[j]
          row <- for_pca[for_pca$muni_i == pcode_i & for_pca$muni_j == pcode_j,]$index
          PCA_1 <- for_pca$PC1[row]
          PCA_2 <- for_pca$PC2[row]
          
          
          # Calculate total distance using our distance formula and save for summary statistics. 
          total_dist <- exp((a*PCA_1+b*PCA_2))
          #summary_db[row] <- total_dist
          
          # If we are using the hiking metric (and both municipalities are among the 1,076 for which we have geographic covariates), calculate and update distance. Where we don't have covariates, the default is d1.  
          if(metric==4 & i<nrow(forcedMigration::geographic_covariates) & j<nrow(forcedMigration::geographic_covariates) ){
            elev_difference <- max((forcedMigration::geographic_covariates$alt_mean[i]-forcedMigration::geographic_covariates$alt_mean[j]),0)
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
    vertex_ids[i]<- forcedMigration::full_municipalities$ADM2_PCODE[i]
  }
  
  # Calculate Dijkstra distances 
  delta_1 = igraph::distances(munigraph,v = which(vertex_ids == epicenter_1),to = igraph::V(munigraph),algorithm = "dijkstra")
  delta_2 = igraph::distances(munigraph,v = which(vertex_ids == epicenter_2),to = igraph::V(munigraph),algorithm = "dijkstra")
  deltas <- data.frame(matrix(c(delta_1,delta_2),ncol = 2))
  colnames(deltas) <- c("delta_1","delta_2");
  deltas <- dplyr::mutate(deltas, delta_min = as.numeric(purrr::map2(delta_1,delta_2,min)))
  
  
  
  # Merge our distances (from specified epicenter) into dataframe containing map polygons. 
  distances <- as.numeric(deltas[,3])
  
  
  # Dataset for merging delta and ring_num into other datasets by ADM2_PCODE. 
  merge_deltas <- data.frame(delta = distances,ADM2_PCODE = vertex_ids)
  
  # Calculate ring as 10 * decile in distance distribution plus 1.  
  merge_deltas <- dplyr::mutate(merge_deltas, ring_num = as.integer(10*ecdf(merge_deltas$delta)(delta))+1)
  
  
  # For the single largest delta in the dataset, decile is 10 and so ring is 11. 
  # Set ring_num to 10 for this case. 
  merge_deltas <- dplyr::mutate(merge_deltas, ring_num = ifelse(ring_num == 11,10,ring_num))
  return(merge_deltas)
}
#' Helper function to create distances given epicenter and metric parameters. 
#' @param merge_deltas Output from generate_distances; contains mapping of municipality to distance and ring under given metric, a, b, and epicenter. 
#' 
#' @return Dataframe containing municipality, ring, share of total violence for ring in given year, and standard deviation of violence across municipalities in a given ring-year. 
#' \itemize{
#' \item \code{FinalWithDeltas} Dataframe containing geographic and final distances. 
#'    \item \code{vdf} Dataframe containing year-ring violence. 
#'    \item \code{sdf} Dataframe containing standard errors of year-ring violence. 
#'    \item \code{sharedf} Dataframe containing final output. 
#' }
#' 
#' 
#' 
get_df <- function(merge_deltas){
  FinalWithDeltas <- merge(forcedMigration::geographic_covariates,merge_deltas,by = "ADM2_PCODE")
  violence_data <- forcedMigration::violence_data[forcedMigration::violence_data$year < 2009,]
  violence_set <- merge(violence_data,FinalWithDeltas)      
  
  #print(summary(FinalWithDeltas$delta))
  
  # Take out municipalities with 0 violence. 
  
  violence_subset <- violence_set[violence_set$year == 2008,]
  has_violence_list <- violence_subset[violence_subset$cum_victims_UR>0,]$ADM2_PCODE
  violence_set <- violence_set[violence_set$ADM2_PCODE %in% has_violence_list,]
  
  
  # For each ring and year, stores total number of deaths in that ring in that year. 
  victimsUR_matrix <- matrix(0,nrow = 10,ncol = 13)
  
  # For each ring and year, stores standard deviation of the set {deaths in municipality in year : municipality in ring} . 
  stderr_matrix <- matrix(0,nrow = 10,ncol = 13)
  
  # Stores the two above quantities in their respective matrices. 
  for(ring in 1:10){
    for(year in 1996:2008){
      victimsUR_matrix[ring,year-1995] <- sum(violence_set[violence_set$year == year & violence_set$ring == ring,]$victims__UR)
      stderr_matrix[ring,year-1995] <- sd(violence_set[violence_set$year == year & violence_set$ring == ring,]$victims__UR)
    }
  }
  
  # Convert matrices to dataframes. 
  vdf <- as.data.frame(t(victimsUR_matrix))
  sdf <- as.data.frame(stderr_matrix)
  
  # Generate share of total ring deaths that occur in a given year by dividing each ring-year death figure by total ring deaths. 
  sharedf <- dplyr::mutate(vdf,V1 = V1/sum(V1),V2 = V2/sum(V2),V3 = V3/sum(V3),V4 = V4/sum(V4), V5 = V5/sum(V5), V6 = V6/sum(V6), V7 = V7/sum(V7), V8 = V8/sum(V8), V9 = V9/sum(V9), V10 = V10/sum(V10))
  sharedf <- as.data.frame(t(sharedf))
  
  cnames = paste(1996:2008)
  colnames(sharedf) <- cnames
  colnames(sdf) <- cnames
  rownames(sharedf) <- 1:10
  rownames(sdf) <- 1:10
  
  # Since the rows represent rings, explicitly create a ring variable corresponding to row number. 
  sharedf <- cbind("ring" = as.numeric(rownames(sharedf)),sharedf)
  
  # Make dataframes long to allow for plotting. 
  sharedf <- tidyr::pivot_longer(sharedf,cnames,names_to = "year",values_to = "share")
  sdf <- cbind("ring" = as.numeric(rownames(sdf)),sdf)
  sdf <- tidyr::pivot_longer(sdf,cnames,names_to = "year",values_to = "sd")
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
  ggplot2::ggplot(data=sharedf,mapping=ggplot2::aes(x=year,y=ring,fill=share))+
    ggplot2::geom_tile()+ggplot2::theme_minimal()+ggplot2::scale_fill_gradient(name="Violence Share",low="darkblue",high="red")+
    ggplot2::scale_x_discrete(breaks = seq(1,10,3),label = paste(seq(1,10,3)))
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
  ggplot2::ggplot(sharedf,ggplot2::aes(x = year,y = share,group = factor(ring)))+ggplot2::geom_point(ggplot2::aes(colour = factor(ring)))+ggplot2::facet_wrap(~ ring)
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
  ggplot2::ggplot(sharedf,ggplot2::aes(x = year,y = share,group = factor(ring)))+ggplot2::geom_point(ggplot2::aes(colour = factor(ring)))+ggplot2::geom_errorbar(ggplot2::aes(ymin = share-sd, ymax = share+sd,colour = factor(ring)))+ggplot2::facet_wrap(~ ring)
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
distance_map <- function(metric,a,b,epicenter_1,epicenter_2){
  merged_deltas <- forcedMigration::generate_distances(metric,a,b,epicenter_1,epicenter_2)
  muni_pol <- forcedMigration::muni_pol
  merged_map <- merge(merged_deltas,muni_pol)
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = sf::st_as_sf(merged_map), ggplot2::aes(fill=ring_num),color = 'grey34',lwd=.05) +
    ggplot2::scale_fill_stepsn(colors = c("red","gold","darkgreen","blue","violet"),values = NULL,space = "Lab",na.value = "grey50",guide = "coloursteps",aesthetics = "fill",breaks = c(1:10),limits= c(1, 10))
}

geography_map <- function(covariate){
    map_with_attributes <- merge(forcedMigration::muni_pol,forcedMigration::geographic_covariates)
    if(covariate == "road"){
    ggplot2::ggplot() +
    ggplot2::geom_sf(data = sf::st_as_sf(map_with_attributes), ggplot2::aes(fill=has_road),color = 'grey34',lwd=.05) +ggplot2::scale_fill_manual(values = c("blue","red"))
}
    else if(covariate == "ruggedness"){
        ggplot2::ggplot()+ggplot2::geom_sf(data = sf::st_as_sf(map_with_attributes),ggplot2::aes(fill = ruggedness_index),color = 'grey34',lwd=.05)+ggplot2::scale_fill_gradient(high = "red",low = "blue")
    }
    
    else if(covariate == "slope"){
        ggplot2::ggplot()+ggplot2::geom_sf(data = sf::st_as_sf(map_with_attributes),ggplot2::aes(fill = slope_mean),color = 'grey34',lwd=.05)+ggplot2::scale_fill_gradient(high = "red",low = "blue")
    }
    
    else if(covariate == "altitude"){
        ggplot2::ggplot()+ggplot2::geom_sf(data = sf::st_as_sf(map_with_attributes),ggplot2::aes(fill = alt_mean),color = 'grey34',lwd=.05)+ggplot2::scale_fill_gradient(high = "red",low = "blue")
    }
    
    else if(covariate == "elevation_difference"){
        ggplot2::ggplot()+ggplot2::geom_sf(data = sf::st_as_sf(map_with_attributes),ggplot2::aes(fill = elevation_difference),color = 'grey34',lwd=.05)+ggplot2::scale_fill_gradient(high = "red",low = "blue")
    }
    
    else if(covariate == "forest"){
        ggplot2::ggplot()+ggplot2::geom_sf(data = sf::st_as_sf(map_with_attributes),ggplot2::aes(fill = is_forested),color = 'grey34',lwd=.05)+ggplot2::scale_fill_manual(values = c("blue","red"))
    }
    
    else{
        print("Please enter covariate in the following list: road, ruggedness, slope, altitude, elevation_difference, forest.")
    }
}



abandoned_land_map <- function(covariate){
    map_with_attributes <- merge(forcedMigration::muni_pol,forcedMigration::abandoned_land)
    if(covariate == "displaced"){
        ggplot2::ggplot()+ggplot2::geom_sf(data = sf::st_as_sf(map_with_attributes),ggplot2::aes(fill = displaced),color = 'grey34',lwd=.05)+ggplot2::scale_fill_gradient(high = "red",low = "blue")
    }
    
    else if(covariate == "hect_paramilitary"){
        ggplot2::ggplot()+ggplot2::geom_sf(data = sf::st_as_sf(map_with_attributes),ggplot2::aes(fill = hect_abandoned_paramilitary),color = 'grey34',lwd=.05)+ggplot2::scale_fill_gradient(high = "red",low = "blue")
    }
    
    else if(covariate == "hect_other_armed"){
        ggplot2::ggplot()+ggplot2::geom_sf(data = sf::st_as_sf(map_with_attributes),ggplot2::aes(fill = hect_abandoned_other_armed),color = 'grey34',lwd=.05)+ggplot2::scale_fill_gradient(high = "red",low = "blue")
    }
    
    else if(covariate == "total_hect"){
        ggplot2::ggplot()+ggplot2::geom_sf(data = sf::st_as_sf(map_with_attributes),ggplot2::aes(fill = total_hect_abandoned),color = 'grey34',lwd=.05)+ggplot2::scale_fill_gradient(high = "red",low = "blue")
    }
    
    
    else{
        print("Please enter covariate in the following list: displaced, hect_paramilitary, hect_other_armed, total_hect.")
    }
}

violence_map <- function(year,covariate){

   panel <-  forcedMigration::panel[forcedMigration::panel$year == year,]
   map_with_attributes <- merge(forcedMigration::muni_pol,panel,by = "municipality")
   
   if(covariate == "V_cum"){
   ggplot2::ggplot() +
   ggplot2::geom_sf(data = sf::st_as_sf(map_with_attributes), ggplot2::aes(fill=V_cum),color = 'grey34',lwd=.05) +ggplot2::scale_fill_gradient(values = c("blue","red"))
}
   else if(covariate == "V_flow"){
       ggplot2::ggplot()+ggplot2::geom_sf(data = sf::st_as_sf(map_with_attributes),ggplot2::aes(fill = V_flow),color = 'grey34',lwd=.05)+ggplot2::scale_fill_gradient(high = "red",low = "blue")
   }
   
   else if(covariate == "D_AS"){
       ggplot2::ggplot()+ggplot2::geom_sf(data = sf::st_as_sf(map_with_attributes),ggplot2::aes(fill = D_AS),color = 'grey34',lwd=.05)+ggplot2::scale_fill_gradient(high = "red",low = "blue")
   }
   
   else if(covariate == "D_CODHES"){
       ggplot2::ggplot()+ggplot2::geom_sf(data = sf::st_as_sf(map_with_attributes),ggplot2::aes(fill = D_CODHES),color = 'grey34',lwd=.05)+ggplot2::scale_fill_gradient(high = "red",low = "blue")
   }
   
   else if(covariate == "D_RUV"){
       ggplot2::ggplot()+ggplot2::geom_sf(data = sf::st_as_sf(map_with_attributes),ggplot2::aes(fill = D_RUV),color = 'grey34',lwd=.05)+ggplot2::scale_fill_gradient(high = "red",low = "blue")
   }
   
   else if(covariate == "D_CEDE"){
       ggplot2::ggplot()+ggplot2::geom_sf(data = sf::st_as_sf(map_with_attributes),ggplot2::aes(fill = D_CEDE),color = 'grey34',lwd=.05)+ggplot2::scale_fill_gradient(high = "red",low = "blue")
   }
   
   else if(covariate == "D_JYP"){
       ggplot2::ggplot()+ggplot2::geom_sf(data = sf::st_as_sf(map_with_attributes),ggplot2::aes(fill = D_JYP),color = 'grey34',lwd=.05)+ggplot2::scale_fill_gradient(high = "red",low = "blue")
   }
   
   else{
       print("Please enter covariate in the following list: V_flow, V_cum, D_AS, D_CODHES, D_RUV, D_CEDE, D_JYP.")
   }
}
     
 get_neighbor_violence <- function(municipality){
 df <- forcedMigration::generate_distances(0,1,1,municipality,municipality)
 neighbors <- df[df$delta == 1,]$ADM2_PCODE
 municipalities <- as.numeric(unlist(lapply(neighbors,substring,first=3)))
 flows <- forcedMigration::violence_data[forcedMigration::violence_data$municipality %in% municipalities,]    
 flows <- subset(flows,select = c(municipality,year,victims__UR))
 return(flows)
 }
   



#'
#' 
#'
#' 
#'
#'
#'
#'
#'
#'
#'
#'
#'
