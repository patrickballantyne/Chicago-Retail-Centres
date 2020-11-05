## Typology Functions

## Function to implement PAM across varying seeds
get_pam <- function(df, seed = 123, k = 3) {
  
  ## Set Seed
  set.seed(seed)
  ## Run PAM
  pm <- pam(x = df, k = k)
  
  ## Extract Silhouette Scores
  output <- cbind(df, pm$clustering)
  colnames(output)[12] <- "clustering"
  output <- silhouette(output$clustering, dist(df))
  
  avg_si <- mean(output[, 3])
  return(avg_si)
}

## Function that runs get_pam for every k value and compiles the results
get_silhouette <- function(db) {
  
  ## Extract values across the 8 k values
  k2 <- as.data.frame(cbind("2", get_pam(db, seed = 123, k = 2)))
  k3 <- as.data.frame(cbind("3", get_pam(db, seed = 123, k = 3)))
  k4 <- as.data.frame(cbind("4", get_pam(db, seed = 123, k = 4)))
  k5 <- as.data.frame(cbind("5", get_pam(db, seed = 123, k = 5)))
  k6 <- as.data.frame(cbind("6", get_pam(db, seed = 123, k = 6)))
  k7 <- as.data.frame(cbind("7", get_pam(db, seed = 123, k = 7)))
  k8 <- as.data.frame(cbind("8", get_pam(db, seed = 123, k = 8)))
  k9 <- as.data.frame(cbind("9", get_pam(db, seed = 123, k = 9)))
  k10 <- as.data.frame(cbind("10", get_pam(db, seed = 123, k = 10)))

  ## Compile and Format
  out <- rbind(k2, k3, k4, k5, k6, k7, k8, k9, k10)
  colnames(out) <- c("k", "silhouette_score")
  out
}


## Function that filters dataset to a specific supergroup, plots the supergroup and returns the silhouette scores
identify_groups <- function(db, cluster_id = 1) {
  
  ## Format for S.Scores and Run
  db <- db %>%
    filter(supergroup_id == cluster_id)
  g <- get_silhouette(db[, 3:22])
  
  ## Format for Clustergram and Run
  db_m <- as.matrix(db[, 3:22])
  t <- clustergram(Data = db_m, k.range = 1:10, line.width = 0.004)
  
  ## Show Scores and Clustergram
  t
  g
  
}

## Subgroup clustering
get_pam_groups <- function(db, cluster_id = 1, k = 2) {
  
  ## Subset data to supergroup of interest
  db <- db %>%
    filter(supergroup_id == cluster_id) %>%
    select(county_cluster_id, supergroup_id, everything())
  
  ## Run the clustering
  pm <- pam(x = db[, 3:22], k = k, metric = "euclidean")
  
  ## Format output
  cl <- as.data.frame(pm$clustering)
  colnames(cl) <- "group_id"
  cols <- as.data.frame(db[, 1])
  out <- cbind(cols, cl)
  colnames(out) <- c("county_cluster_id", "group_id")
  db_out <- merge(db, out, by = "county_cluster_id")
  
  ## Final Formatting
  db_out <- db_out %>%
    mutate(group_id = paste(supergroup_id,".", group_id)) %>%
    select(county_cluster_id, supergroup_id, group_id, everything())
  return(db_out)
  
}

## Subgroup medoid values
get_group_medoids <- function(db, cluster_id = 1, k = 2) {
  
  ## Subset data to supergroup of interest
  db <- db %>%
    filter(supergroup_id == cluster_id) %>%
    select(county_cluster_id, supergroup_id, everything())
  
  ## Run the clustering
  pm <- pam(x = db[, 3:22], k = k, metric = "euclidean")
  
  ## Format medoids
  medoids <- as.data.frame(t(pm$medoids))
  medoids <- medoids %>%
    rownames_to_column()
  medoids <- gather(medoids, key = group_id, value = cluster_vals, -1)
  medoids$group_id <- gsub("V", "", paste(medoids$group_id))
  colnames(medoids)[1] <- "variable"
  medoids <- medoids %>%
    mutate(pos = cluster_vals >= 0)
  return(medoids)
}
