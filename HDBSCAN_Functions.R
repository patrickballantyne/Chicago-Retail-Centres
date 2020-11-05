##### Functions for HDBSCAN


### Function to obtain and format euclidean clusters
get_euclidean_clusters <- function(sf, minPts) {
  
  ## Run the clustering
  cl <- dbscan::hdbscan(st_coordinates(sf), minPts = minPts)
  ## Extract cols
  results <- cbind(sf ,paste0(cl$cluster))
  results <- cbind(results, paste0(cl$membership_prob))
  results <- cbind(results, paste0(cl$outlier_scores))
  
  ## Format
  ## FOrmat output    
  results <- dplyr::rename_at(results, "paste0.cl.cluster.",~"clustering")
  results$clustering <- as.numeric(as.character(results$clustering))
  results <- dplyr::rename_at(results, "paste0.cl.membership_prob.",~"membership_prob")
  results$membership_prob <- as.numeric(as.character(results$membership_prob))
  results <- dplyr::rename_at(results, "paste0.cl.outlier_scores.",~"outlier_scores")
  results$outlier_scores <- as.numeric(as.character(results$outlier_scores))
  
  ## Drop noise
  results <- results[results$clustering != 0,]
  
  ## Convert to sf
  results <- st_as_sf(results, coords = c("X", "Y"))
  return(results)
}


### Function for obtaining network clusters
get_network_clusters <- function(sf, dist, minPts) {
  
  ## Run the clustering
  cl <- dbscan::hdbscan(st_coordinates(sf), minPts = minPts, xdist = dist)
  ## Extract cols
  results <- cbind(sf ,paste0(cl$cluster))
  results <- cbind(results, paste0(cl$membership_prob))
  results <- cbind(results, paste0(cl$outlier_scores))
  
  ## Format
  ## FOrmat output    
  results <- dplyr::rename_at(results, "paste0.cl.cluster.",~"clustering")
  results$clustering <- as.numeric(as.character(results$clustering))
  results <- dplyr::rename_at(results, "paste0.cl.membership_prob.",~"membership_prob")
  results$membership_prob <- as.numeric(as.character(results$membership_prob))
  results <- dplyr::rename_at(results, "paste0.cl.outlier_scores.",~"outlier_scores")
  results$outlier_scores <- as.numeric(as.character(results$outlier_scores))
  
  ## Drop noise
  results <- results[results$clustering != 0,]
  
  ## Convert to sf
  results <- st_as_sf(results, coords = c("X", "Y"))
  return(results)
}



### Function for generating buffered hulls based on clustering
get_hulls <- function(sf) {
  sf <- sf %>%
    rename(clustering = county_cluster_id)
  sf <- sf[, c("clustering")]
  polygons <- map(unique(sf$clustering), 
                  ~ concaveman(sf[sf$clustering %in% .,])
  ) %>%
    map2(unique(sf$clustering), ~ mutate(.x, clustering = .y)) %>%
    reduce(rbind)
  poly <- st_make_valid(polygons)
  poly <- st_transform(poly, 32616)
  buffered_boundaries <- st_buffer(poly, dist = 50)
  buffered_boundaries <- st_transform(buffered_boundaries, 4326)
  buffered_boundaries <- buffered_boundaries %>%
    arrange(clustering)
  return(buffered_boundaries)
}

## Function to convert hulls into H3 geoms
get_h3_hulls <- function(sf, h3_res = 14) {
  
  ## Extract H3 geoms from hull polygon
  h3 <- polyfill(sf, res = h3_res, simple = FALSE)
  
  ## Extract list of h3 addresses for cluster
  lookup <- stack(h3$h3_polyfillers)
  lookup$ind <- as.numeric(as.character(lookup$ind))
  lookup$cluster_id <- scales::rescale(lookup$ind, to = c(1, n_distinct(lookup$cluster_id)))
  lookup <- select(lookup, -c(ind))
  lookup <- rename(lookup, h3_address = values)
  
  ## Convert hexagons to polygon
  h3 <- h3_to_polygon(unlist(h3$h3_polyfillers), simple = FALSE)
  
  ## Merge cluster info
  h3_merge <- merge(h3, lookup, by = "h3_address", all.x = TRUE)
  h3_merge <- h3_merge %>%
    arrange(cluster_id)
  
  ## Dissolve into boundary
  dissolved_h3 <- st_dissolve(h3_merge, by = "cluster_id")
  ##  Sort out columns
  h3 <- dissolved_h3 %>%
    rownames_to_column(var = "clustering") %>%
    select(-c(cluster_id))
  return(h3)
  
}

## Function that joins the get_hulls and get_h3_hulls together
get_retail_centre_extent <- function(sf) {
  
  ## Function 1. Get Concave Hulls
  ch <- get_hulls(sf)
  ## Combine adjacent polygons
  ch_u <-
    ch %>%
    st_union(by_feature = FALSE) %>%
    st_as_sf()  %>%
    st_cast(to = "POLYGON") %>%
    rownames_to_column(var = "id") %>%
    st_transform(crs = 32616)
  
  ## Join Points to Polygons and Extract Lookup for Later
  lookup <- st_join(ch_u, sf)
  lookup <- lookup %>%
    as.data.frame() %>%
    select(safegraph_place_id, top_category, sub_category, ldc_aggregation, state_code, county_name, county_code, code_abc, id) %>%
    rename(clustering = id, county_code_abc = code_abc)
  
  ## Function 2.
  h3 <- get_h3_hulls(ch_u)
  h3 <- merge(h3, lookup)
  h3$clustering <- as.numeric(h3$clustering)
  h3 <- h3 %>%
    arrange(clustering)
  

  ## Calculate total units per cluster
  unit_count <- h3 %>%
    as.data.frame() %>%
    select(clustering) %>%
    group_by(clustering) %>%
    summarise(cluster_unit_count = n())
    
  ## Merge on unit count
  polygons <- merge(h3, unit_count, by = "clustering", all.x = TRUE)
  polygons <- polygons %>%
    select(-clustering, everything())
  
  ## Paste county code infront
  polygons$county_cluster_id <- as.character(with(polygons, paste(county_code_abc, clustering, sep = "")))
  polygons <- select(polygons, -c("clustering", "state_code","county_code","county_code_abc"))
  polygons <- polygons %>%
    select(safegraph_place_id, top_category, sub_category, ldc_aggregation, county_cluster_id, county_name, cluster_unit_count)
  return(polygons)
  
}

## Function for retail centre boundaries (no H3 geoms)
get_retail_centre_extent_basic <- function(sf) {
  
  ## Get concave hulls
  ch <- get_hulls(sf)  
  
  ## Join Adjacent Polygons
  ch_u <-
    ch %>%
    st_union(by_feature = FALSE) %>%
    st_as_sf()  %>%
    st_cast(to = "POLYGON") %>%
    rownames_to_column(var = "id") %>%
    st_transform(crs = 32616)
  
  ## Spatial Join Points to Polygons
  ch_u <- st_join(ch_u, sf)
  ch_u <- ch_u %>%
    select(safegraph_place_id, top_category, sub_category, ldc_aggregation, state_code, county_name, county_code, code_abc, id) %>%
    rename(clustering = id, county_code_abc = code_abc)
  
  ## Get unit count
  unit_count <- ch_u %>%
    as.data.frame() %>%
    select(clustering) %>%
    group_by(clustering) %>%
    summarise(cluster_unit_count = n())
  
  ## Merge on unit count
  polygons <- merge(ch_u, unit_count, by = "clustering", all.x = TRUE)
  polygons <- polygons %>%
    select(-clustering, everything())
  
  ## Adding region letter to cluster ID
  polygons$county_cluster_id <- as.character(with(polygons, paste(county_code_abc, clustering, sep = "")))
  polygons <- select(polygons, -c("clustering", "state_code","county_code","county_code_abc"))
  polygons <- polygons %>%
    select(safegraph_place_id, top_category, sub_category, ldc_aggregation, county_cluster_id, county_name, cluster_unit_count)
  return(polygons)
}

get_cook_retail_centre_extent_basic <- function(sf, side = "N") {
  
  ## Get concave hulls
  ch <- get_hulls(sf)  
  
  ## Join Adjacent Polygons
  ch_u <-
    ch %>%
    st_union(by_feature = FALSE) %>%
    st_as_sf()  %>%
    st_cast(to = "POLYGON") %>%
    rownames_to_column(var = "id") %>%
    st_transform(crs = 32616)
  
  ## Spatial Join Points to Polygons
  ch_u <- st_join(ch_u, sf)
  ch_u <- ch_u %>%
    select(safegraph_place_id, top_category, sub_category, ldc_aggregation, state_code, county_name, county_code, code_abc, id) %>%
    rename(clustering = id, county_code_abc = code_abc)
  
  ## Get unit count
  unit_count <- ch_u %>%
    as.data.frame() %>%
    select(clustering) %>%
    group_by(clustering) %>%
    summarise(cluster_unit_count = n())
  
  ## Merge on unit count
  polygons <- merge(ch_u, unit_count, by = "clustering", all.x = TRUE)
  polygons <- polygons %>%
    select(-clustering, everything())
  
  ## Adding region letter to cluster ID
  polygons$county_cluster_id <- as.character(with(polygons, paste(county_code_abc, clustering, sep = "")))
  polygons <- select(polygons, -c("clustering", "state_code","county_code","county_code_abc"))
  polygons <- polygons %>%
    select(safegraph_place_id, top_category, sub_category, ldc_aggregation, county_cluster_id, county_name, cluster_unit_count)
  
  ## Paste identifier
  polygons$county_cluster_id <- as.character(with(polygons, paste(side, county_cluster_id, sep = "")))
  return(polygons)
}


### Get Polygons from Clustering
get_retail_polygons <- function(sf) {
  
  ## Extract columns
  sf[, c("county_cluster_id", "county_name", "cluster_unit_count")]
  ## Extract one polygon per cluster
  sf <- st_difference(sf)
  return(sf)
}

### Get Points from Clustering
get_retail_points <- function(sf, pts) {
  
  ## Extract info from clusters
  sf<- sf %>%
    as.data.frame() %>%
    select(safegraph_place_id, county_cluster_id)
  ## Merge with original pts to extract sfc POINT
  sf <- merge(pts, sf, by = "safegraph_place_id", all.y = TRUE)
  
  ## Format before extraction
  pts <- sf[, c("safegraph_place_id", "top_category", "sub_category", "ldc_aggregation", "state_code", "county_name", "county_cluster_id", "cluster_unit_count")]
  return(pts)
}


## Function for calculating and returning distance to cluster cores for a clustered dataset
get_distance_to_core <- function(sf) {
  
  ## Extract Core Points
  core_pts <- sf %>%
    group_by(clustering) %>%
    filter(membership_prob == max(membership_prob))
  
  ## Format for osrmTable()
  pts <- sf[, c("safegraph_place_id", "clustering")]
  core_pts <- core_pts[, c("safegraph_place_id", "clustering")]
  
  ## Get the ID's set up
  pts <- pts %>%
    arrange(clustering)
  pts$ID <- as.numeric(as.character(rownames(pts)))
  
  ## Calculate the network durations from the cluster cores to all the other points
  durations <- osrmTable(src = core_pts, dst = pts, measure = "duration")
  
  ## Format output
  durations_df <- as.data.frame(durations$durations)
  durations_df$cluster_cores <- rownames(durations_df)
  durations_df <- durations_df %>%
    select(cluster_cores, everything())
  durations_df$cluster_cores <- gsub("X", "", as.character(durations_df$cluster_cores), n)
  durations_df <- gather(durations_df, cluster_cores, value = Mins_To_Core)
  colnames(durations_df) <- c("ID", "Mins_To_Core")
  
  cluster_cores <- as.data.frame(core_pts[, c("clustering")])
  colnames(cluster_cores) <- c("clustering")
  
  ## Joining outputs together
  results <- cbind.fill(durations_df, cluster_cores) 
  colnames(results) <- c("ID", "Duration_to_Core", "Clustering")
  results <- sp::merge(pts, results, by = "ID", all.y = TRUE)
  colnames(results)[6] <- "to_delete"
  results <- select(results, -c("to_delete"))
  
  ## Drop the rows that have distances calculated but not for their cluster assigned
  results <- results %>%
    filter(Clustering == clustering)
  results$ID <- as.numeric(as.factor(results$ID))
  results <- results[, c("safegraph_place_id", "clustering", "Duration_to_Core")]
  results <- distinct(results, safegraph_place_id, .keep_all = TRUE)
  return(results)
  
}


## Function that takes the initial clustering with network distances and refines the clusters, keeping only the points withi
## 3mins walk of the core cluster point and those with 5 points minimum in the cluster
get_refined_network_clusters <- function(sf) {
  
  ## Get the good clusters straight out
  sf_df <- get_distance_to_core(sf)
  ## Modify original data
  sf <- sf %>%
    as.data.frame() %>%
    select(-c(geom))
  good_clusters <- sf_df[sf_df$Duration_to_Core < 3,]
  
  group <- good_clusters %>%
    group_by(clustering) %>%
    summarise(cluster_unit_count = n()) %>%
    as.data.frame() %>%
    select(-c("geometry"))
  
  good_clusters <- merge(good_clusters, group, by = "clustering", all.x = TRUE)
  good_clusters <- good_clusters[good_clusters$cluster_unit_count >= 10,]
  
  clusters <- merge(good_clusters, sf, by = "safegraph_place_id", all.x = TRUE)
  clusters <- clusters[, c("safegraph_place_id", "top_category", "sub_category", "state_code", "county_code", "county_name",
                           "ldc_aggregation", "clustering.x", "membership_prob", "outlier_scores", "Duration_to_Core", "cluster_unit_count")]
  clusters <- clusters %>%
    rename(clustering = clustering.x, duration_to_core = Duration_to_Core)
  return(clusters)
  
  # ## Extract the 'bad' clusters
  # ### Firstly, those where duration to core is more than 3mins
  # long_distance_points <- sf[sf$Duration_to_Core >= 3,]
  # long_distance_points <- select(long_distance_points, -c("clustering", "Duration_to_Core"))
  # long_distance_points <- drop_na(long_distance_points)
  # 
  # ### Secondly, the points with less than 5 pts in cluster as a result of above removal
  # tidy_points <- sf[sf$Duration_to_Core < 3,]
  # any_insignificant <- tidy_points %>%
  #   group_by(clustering) %>%
  #   summarise(count = n()) %>%
  #   mutate(cluster_unit_count = count) %>%
  #   as.data.frame() %>%
  #   select(-c("geometry", "count"))
  # tidy_points <- merge(tidy_points, any_insignificant, by = "clustering", all.x = TRUE)
  # small_points <- tidy_points[tidy_points$cluster_unit_count < 10,]
  # small_points <- select(small_points, -c("clustering", "Duration_to_Core", "cluster_unit_count"))
  # 
  # ## Join the two sets together to have a comprehensive dataset of 'unclustered points'bad' points
  # bad_clusters <- rbind(long_distance_points, small_points)
  # bad_clusters <- drop_na(bad_clusters)
  # 
  # ## Recalculate distance matrices for the 'bad' clusters
  # bad_dist <- osrmTable(loc = bad_clusters, measure = "distance")
  # bad_dist <- bad_dist$distances
  # bad_dist <- as.dist(bad_dist)
  # 
  # ## Re-run clustering
  # additional_clustering <- get_network_clusters(bad_clusters, bad_dist, 10)
  # #return(additional_clustering)
  # 
  # ## Extract clusters that satisfy the two conditions
  # #additional_clustering <- get_distance_to_core(additional_clustering)
  # #additional_clustering <- get_good_clusters(additional_clustering)
  # 
  # 
  # ## Joining the original good clusters with the additional good ones
  # ### Format and Join
  # additional_clustering$clustering <- sapply(0, paste0, additional_clustering$clustering)
  # all_good_clusters <- rbind(good_clusters, additional_clustering)
  # all_good_clusters <- all_good_clusters %>%
  #   arrange(desc(cluster_unit_count))
  # ### Setup new cluster IDs
  # new_cluster_ids <- all_good_clusters %>%
  #   as.data.frame() %>%
  #   select(clustering) %>%
  #   group_by(clustering) %>%
  #   distinct()
  # new_cluster_ids$new_cluster_id <- seq.int(nrow(new_cluster_ids))
  # clusters <- merge(all_good_clusters, new_cluster_ids, by = "clustering", all.x = TRUE)
  # clusters <- select(clusters, - c("clustering"))
  # clusters <- clusters %>% # Format for join
  #   as.data.frame() %>%
  #   select(-c(geometry))
  # 
  # ## Merge back on the original POI data
  #clusters <- merge(clustering, clusters, by = "safegraph_place_id", all.y = TRUE)
  # clusters <- select(clusters, -c("clustering", "membership_prob", "outlier_scores"))
  # #return(clusters)
}

## Function for getting clusters ready for join
clean_clusters <- function(cluster_object = cook_north, cluster_region = "Cook North-Side") {
  
  ## Create column identifier for cluster region
  cluster_object$cluster_region <- cluster_region
  ## Paste Cluster Region in front of clustering
  cluster_object <- cluster_object %>%
    distinct(safegraph_place_id, .keep_all = TRUE)
  return(cluster_object)
}
