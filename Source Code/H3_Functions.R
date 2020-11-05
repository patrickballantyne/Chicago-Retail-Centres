## H3 Refining Functions

### Function that refines the clustered points from HDBSCAN 
refine_clusters <- function(cl, h3_res = 11, kring = 1, county_identifier = NULL) {
  
  ## STAGE ONE - Convert all clustered points to individual H3 tiles
  ## #Convert pts to H3
  pts_h3 <- point_to_h3(cl, res = h3_res, simple = FALSE)
  ### Extract list of all H3 addresses in cluster
  hexes_list <- unlist(pts_h3$h3_resolution_11, use.names = TRUE)
  hexes <- h3_to_polygon(hexes_list)
  hexes <- st_as_sf(hexes)
  ### Join data back on
  h3_polygons <- cbind(hexes, cl)
  h3_polygons <- h3_polygons %>%
    select(-c(geom))
  ### Join on list of h3 addresses
  h3_polygons <- cbind(h3_polygons, hexes_list)
  ### Columns
  h3_polygons <- h3_polygons %>%
    rename(h3_address = hexes_list, geom = x)
  
  ## STAGE TWO - Aggregate the clusters to single polygons
  ### Count no. distinct points in polygons
  list <- h3_polygons %>%
    as.data.frame() %>%
    select(h3_address, safegraph_place_id, clustering) %>%
    group_by(h3_address) %>%
    summarise(n_pts = n())
  ### Merge onto main dataset
  db <- h3_polygons %>%
    select(clustering, h3_address)
  db <- merge(db, list, by = "h3_address", all.x = TRUE)
  db <- db %>%
    distinct(.keep_all = TRUE)
  
  ## STAGE THREE - Mini-Function for extracting krings 
  cluster2kring <- function(db, kring = kring) {
    
    ## Get kring addresses
    kring_list <- get_kring(h3_address = db$h3_address, ring_size = kring)
    ## Convert list of addresses to a polygon
    kring_sf <- h3_to_polygon(unlist(kring_list, use.names = TRUE), simple = FALSE)
    return(kring_sf)
  }
  
  ## STAGE FOUR
  ### Apply cluster2kring across the different clustering elements and extract polygon
  out <- db %>%
    group_by(clustering) %>%
    group_map(~ cluster2kring(.x, kring = kring))
  ## Format output - extracting single polygons
  out <- Map(cbind, out, cluster_id = 1:length(out))
  out <- do.call(rbind, out)
  out <- out %>%
    summarise()
  out_poly <- st_cast(out, to = "POLYGON")
  out_poly <- out_poly %>%
    mutate(polygon_id = row_number()) %>%
    select(polygon_id, geometry)
  
  ## STAGE FIVE - Extracting the good clusters
  ### Sort out transformations
  poly <- st_transform(out_poly, crs = 32616)
  pts <- st_transform(cl, crs = 32616)
  # Perform an intersection
  int <- st_intersection(poly, pts)
  # Calculate the number of pts in each distinct polygon
  int_df <- int %>%
    as.data.frame() %>%
    group_by(polygon_id) %>%
    summarise(n_pts = n()) %>%
    filter(n_pts >= 10)
  # Merge on to extract contiguous polygons
  out <- merge(poly, int_df, by = "polygon_id", all.x = TRUE)
  out <- drop_na(out)
  # Remove any holes in polygons
  out <- nngeo::st_remove_holes(out)
  
  ## STAGE SIX - Assigning Cluster ID's
  ### Sort by number of points in cluster and then assign a unique id based on this
  out <- out %>%
    arrange(desc(n_pts)) %>%
    mutate(cluster_id = row_number())
  ### Paste the county identifier in front of them
  out$county_cluster_id <- as.character(with(out, paste(county_identifier, cluster_id, sep = "")))
  ### Final formatting
  out <- out %>%
    select(county_cluster_id, n_pts, geom)
  return(out)
}

## Function that extracts the refined points for the retail centre
get_refined_cluster_points <- function(poly, pts) {
  
  ## Format points
  pts <- pts %>%
    select(-c(clustering, membership_prob, outlier_scores, duration_to_core, cluster_unit_count))
  
  ## Perform intersection
  r_pts <- st_intersection(pts, poly)
  r_pts <- r_pts %>%
    rename(cluster_pts = n_pts)
  return(r_pts)
  
}
