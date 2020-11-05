### FUNCTION FOR CALCULATING DISTANCES FROM POINTS TO CLUSTER CORE

distance_to_core <- function(pts) {
  
  ##########################################################################
  
  # PART 1 - FORMATTING DATA 
  
  ##########################################################################
  ## Extract the Cluster 'Cores'
  core_pts <- pts %>%
    group_by(clustering) %>%
    filter(membership_prob == max(membership_prob))
  
  ## Format for osrmTable()
  pts <- pts[, c("safegraph_place_id", "clustering")]
  core_pts <- core_pts[, c("safegraph_place_id", "clustering")]
  
  pts_df <- unlist(st_geometry(pts)) %>% 
    matrix(ncol=2,byrow=TRUE) %>% 
    as_tibble() %>% 
    setNames(c("lon","lat"))
  pts_df <- as.data.frame(bind_cols(pts, pts_df))
  pts_df <- select(pts_df, -c("geometry"))
  pts_df <- pts_df[, c(2, 3:4, 1)]
  
  core_df <- unlist(st_geometry(core_pts)) %>% 
    matrix(ncol=2,byrow=TRUE) %>% 
    as_tibble() %>% 
    setNames(c("lon","lat"))
  core_df <- as.data.frame(bind_cols(core_pts, core_df))
  core_df <- select(core_df, -c("geometry"))
  core_df <- core_df[, c(2, 3:4, 1)]
  
  ## Assigning IDs to the points, to allow them to be returned later
  pts_df <- pts_df[order(pts_df$clustering),]
  pts_df$ID <- seq.int(nrow(pts_df))
  pts_full <- pts_df
  pts_df <- pts_df[, c("ID", "lon", "lat", "safegraph_place_id")]
  
  
  ###############################################################################
  
  # PART 2 - CALCULATING DISTANCES
  
  ################################################################################
  
  ## Run the calculation
  durations <- osrmTable(src = core_df, dst = pts_df, measure = "duration")
  
  ## Format output
  durations_df <- as.data.frame(durations$durations)
  durations_df$cluster_cores <- rownames(durations_df)
  durations_df <- durations_df %>%
    select(cluster_cores, everything())
  durations_df$cluster_cores <- gsub("X", "", as.character(durations_df$cluster_cores), n)
  durations_df <- gather(durations_df, cluster_cores, value = Mins_To_Core)
  colnames(durations_df) <- c("ID", "Mins_To_Core")
  
  cluster_cores <- as.data.frame(core_df[, c("clustering")])
  colnames(cluster_cores) <- c("clustering")
  
  ## Joining outputs together
  results <- cbind.fill(durations_df, cluster_cores) 
  colnames(results) <- c("ID", "Duration_to_Core", "Clustering")
  
  results <- merge(results, pts_full, by = "ID", all.x = TRUE)
  results <- select(results, -c("lon", "lat"))
  
  ## Drop the rows that have distances calculated but not for their cluster assigned
  results <- results %>%
    filter(Clustering == clustering)
  results$ID <- as.numeric(as.factor(results$ID))
  
  results <- results[, c("safegraph_place_id", "clustering", "Duration_to_Core")]
  results <- distinct(results, safegraph_place_id, .keep_all = TRUE)
  return(results)
}



get_untidy_points <- function(sf) {
  
  ## Extract the points not within 3 mins walk of cluster core
  long_distance_points <- sf[sf$Duration_to_Core >= 3,]
  long_distance_points <- select(long_distance_points, -c("cluster_region", "membership_prob", "outlier_scores",
                                                   "clustering", "Duration_to_Core"))
  long_distance_points <- drop_na(long_distance_points)
  long_distance_points <- long_distance_points %>%
    rename(geometry = geom )
  
  ## Extract the points with less than 5 pts in cluster as a result of above removal
  tidy_points <- sf[sf$Duration_to_Core < 3,]
  any_insignificant <- tidy_points %>%
    group_by(clustering) %>%
    summarise(count = n()) %>%
    mutate(cluster_unit_count = count) %>%
    as.data.frame() %>%
    select(-c("geom", "count"))
  
  tidy_points <- merge(tidy_points, any_insignificant, by = "clustering", all.x = TRUE)
  small_points <- tidy_points[tidy_points$cluster_unit_count < 5,]
  small_points <- select(small_points, -c("cluster_region", "membership_prob", "outlier_scores",
                                          "clustering", "Duration_to_Core", "cluster_unit_count"))
  
  ## Join the two sets together to have a comprehensive dataset of unclustered points
  all_untidy_points <- rbind(long_distance_points, small_points)
  all_untidy_points <- drop_na(all_untidy_points)
  return(all_untidy_points)
}
