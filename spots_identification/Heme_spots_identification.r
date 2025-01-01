library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

args <- commandArgs(TRUE)
sample <- args[1]
n_cluster <- as.numeric(args[2])
cosine_value_threshold <- as.numeric(args[3])
patient_threshold <- as.numeric(args[4])
neighbors_threshold <- as.numeric(args[5])
max_percentage <- as.numeric(args[6])

# sample = "control"
# n_cluster = 7
# cosine_value_threshold = 0.05
# patient_threshold = 10
# neighbors_threshold = 4
# max_percentage = 0.05

#
st_data <- Load10X_Spatial(
    paste0("/Users/haoran/Documents/MGGNN/data/Heme/input/", sample),
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = sample
)
spots <- GetTissueCoordinates(st_data)
spots <- spots[, 1:2]
names(spots)[1:2] <- c("imagerow", "imagecol")

anatomy = read.csv(paste0("/Users/haoran/Documents/MGGNN/data/Heme/input/", sample, "/", sample, "_anatomy.csv"))
anatomy$anatomy[anatomy$anatomy == ""] = "NA"
idx_notnull = anatomy$anatomy != "NA"

spots$anatomy = anatomy$anatomy

spots$anatomy[spots$anatomy == "hypothalamuis"] # = "hypothalamus"



#
marker_genes <- read.csv("/Users/haoran/Documents/MGGNN/data/spots_identification/Heme_marker_genes.csv")

Group1_marker <- marker_genes$Gene[marker_genes$caudate_putamen == 1]
Group2_marker <- marker_genes$Gene[marker_genes$corpus_callosum == 1]
Group3_marker <- marker_genes$Gene[marker_genes$cortex == 1]
Group4_marker <- marker_genes$Gene[marker_genes$globus_pallidus == 1]
Group5_marker <- marker_genes$Gene[marker_genes$hypothalamus == 1]
Group6_marker <- marker_genes$Gene[marker_genes$plexus == 1]
Group7_marker <- marker_genes$Gene[marker_genes$thalamus == 1]

o_1 <- unique(c(Group2_marker, Group3_marker, Group4_marker, Group5_marker, Group6_marker, Group7_marker))
o_2 <- unique(c(Group1_marker, Group3_marker, Group4_marker, Group5_marker, Group6_marker, Group7_marker))
o_3 <- unique(c(Group1_marker, Group2_marker, Group4_marker, Group5_marker, Group6_marker, Group7_marker))
o_4 <- unique(c(Group1_marker, Group2_marker, Group3_marker, Group5_marker, Group6_marker, Group7_marker))
o_5 <- unique(c(Group1_marker, Group2_marker, Group3_marker, Group4_marker, Group6_marker, Group7_marker))
o_6 <- unique(c(Group1_marker, Group2_marker, Group3_marker, Group4_marker, Group5_marker, Group7_marker))
o_7 <- unique(c(Group1_marker, Group2_marker, Group3_marker, Group4_marker, Group5_marker, Group6_marker))

Group1_marker <- setdiff(Group1_marker, o_1)
Group2_marker <- setdiff(Group2_marker, o_2)
Group3_marker <- setdiff(Group3_marker, o_3)
Group4_marker <- setdiff(Group4_marker, o_4)
Group5_marker <- setdiff(Group5_marker, o_5)
Group6_marker <- setdiff(Group6_marker, o_6)
Group7_marker <- setdiff(Group7_marker, o_7)

print(strrep("#", 32))
print(paste0("Layer1_marker_genes(unique): ", length(Group1_marker)))
print(paste0("Layer2_marker_genes(unique): ", length(Group2_marker)))
print(paste0("Layer3_marker_genes(unique): ", length(Group3_marker)))
print(paste0("Layer4_marker_genes(unique): ", length(Group4_marker)))
print(paste0("Layer5_marker_genes(unique): ", length(Group5_marker)))
print(paste0("Layer6_marker_genes(unique): ", length(Group6_marker)))
print(paste0("Layer7_marker_genes(unique): ", length(Group7_marker)))
print(strrep("#", 32))



#
#
#
st_data <- SCTransform(st_data, assay = "Spatial")
st_data <- RunPCA(st_data, assay = "SCT")

count_norm <- st_data[["SCT"]]@data
st_pca <- Embeddings(st_data, reduction = "pca")[, 1:30]
#
#
#



idx <- rownames(count_norm) %in% Group1_marker
count_norm_Group1_marker <- count_norm[idx, , drop = FALSE]

idx <- rownames(count_norm) %in% Group2_marker
count_norm_Group2_marker <- count_norm[idx, , drop = FALSE]

idx <- rownames(count_norm) %in% Group3_marker
count_norm_Group3_marker <- count_norm[idx, , drop = FALSE]

idx <- rownames(count_norm) %in% Group4_marker
count_norm_Group4_marker <- count_norm[idx, , drop = FALSE]

idx <- rownames(count_norm) %in% Group5_marker
count_norm_Group5_marker <- count_norm[idx, , drop = FALSE]

idx <- rownames(count_norm) %in% Group6_marker
count_norm_Group6_marker <- count_norm[idx, , drop = FALSE]

idx <- rownames(count_norm) %in% Group7_marker
count_norm_Group7_marker <- count_norm[idx, , drop = FALSE]



library(lsa) # cosine

min_num <- dim(spots)[1] * 0.01

sample_identification <- function(count_norm, st_pca, spots, cosine_value_threshold, patient_threshold, neighbors_threshold) {
    #
    spot_rank <- apply(count_norm, 1, rank)
    spot_rank <- apply(spot_rank, 1, sum)
    spot_rank <- names(sort(spot_rank, decreasing = TRUE))
    
    cell_idx = 1
    patient = 1
    patient_id = c()
    
    # filter
    for (cell in spot_rank) {
        
        if (cell == spot_rank[1]) {
            d_idx <- which(rownames(spots) == cell)
            d_df <- spots[-d_idx, c("imagerow", "imagecol")]
            d_q <- spots[d_idx, c("imagerow", "imagecol")]
            distance <- sqrt(rowSums((d_df - do.call(rbind,replicate(nrow(d_df), d_q, simplify = FALSE)))**2))
            neighbors <- names(sort(distance)[1:neighbors_threshold])
            neighbors <- st_pca[rownames(st_pca) %in% neighbors, ]
            neighbors_mean <- apply(neighbors, 2, mean)
            mean_pca <- st_pca[rownames(st_pca) == cell, ]
            current_loc <- d_q
            next
        }
        
        cosine_value <- cosine(mean_pca, st_pca[rownames(st_pca) == cell, ])
        
        d_idx <- which(rownames(spots) == cell)
        d_df <- spots[-d_idx, c("imagerow", "imagecol")]
        d_q <- spots[d_idx, c("imagerow", "imagecol")]
        distance <- sqrt(rowSums((d_df - do.call(rbind,replicate(nrow(d_df), d_q, simplify = FALSE)))**2))
        neighbors <- names(sort(distance)[1:neighbors_threshold])
        neighbors <- st_pca[rownames(st_pca) %in% neighbors, ]
        
        cosine_value_2 <- cosine(neighbors_mean, apply(neighbors, 2, mean))
        
        unit_distance <- min(distance)
        cell_distance = sqrt((current_loc$row - d_q$row)**2 + (current_loc$col - d_q$col)**2)
        
        if (cell_idx > min_num & patient > patient_threshold) {break}
        # if (cell_distance > 100*unit_distance) {patient = patient + 1; patient_id = c(patient_id, cell); print("*"); next}
        # if (cosine_value < cosine_value_threshold) {patient = patient + 1; patient_id = c(patient_id, cell); next}
        if (cosine_value_2 < cosine_value_threshold) {patient = patient + 1; patient_id = c(patient_id, cell); next}
        
        cell_idx = cell_idx + 1
        
        mean_pca <- (mean_pca + st_pca[rownames(st_pca) == cell, ]) / 2
    }
    spot_rank <- spot_rank[1:cell_idx]
    spot_rank <- setdiff(spot_rank, patient_id)
    return(spot_rank)
}

Group1_spot_rank <- sample_identification(count_norm_Group1_marker, st_pca, spots, cosine_value_threshold, patient_threshold, neighbors_threshold)
Group2_spot_rank <- sample_identification(count_norm_Group2_marker, st_pca, spots, cosine_value_threshold, patient_threshold, neighbors_threshold)
Group3_spot_rank <- sample_identification(count_norm_Group3_marker, st_pca, spots, cosine_value_threshold, patient_threshold, neighbors_threshold)
Group4_spot_rank <- sample_identification(count_norm_Group4_marker, st_pca, spots, cosine_value_threshold, patient_threshold, neighbors_threshold)
Group5_spot_rank <- sample_identification(count_norm_Group5_marker, st_pca, spots, cosine_value_threshold, patient_threshold, neighbors_threshold)
Group6_spot_rank <- sample_identification(count_norm_Group6_marker, st_pca, spots, cosine_value_threshold, patient_threshold, neighbors_threshold)
Group7_spot_rank <- sample_identification(count_norm_Group7_marker, st_pca, spots, cosine_value_threshold, patient_threshold, neighbors_threshold)

# ceiling
sample_ceiling <- floor(dim(spots)[1] * max_percentage)

if (length(Group1_spot_rank) > sample_ceiling) {Group1_spot_rank = Group1_spot_rank[1:sample_ceiling]}
if (length(Group2_spot_rank) > sample_ceiling) {Group2_spot_rank = Group2_spot_rank[1:sample_ceiling]}
if (length(Group3_spot_rank) > sample_ceiling) {Group3_spot_rank = Group3_spot_rank[1:sample_ceiling]}
if (length(Group4_spot_rank) > sample_ceiling) {Group4_spot_rank = Group4_spot_rank[1:sample_ceiling]}
if (length(Group5_spot_rank) > sample_ceiling) {Group5_spot_rank = Group5_spot_rank[1:sample_ceiling]}
if (length(Group6_spot_rank) > sample_ceiling) {Group6_spot_rank = Group6_spot_rank[1:sample_ceiling]}
if (length(Group7_spot_rank) > sample_ceiling) {Group7_spot_rank = Group7_spot_rank[1:sample_ceiling]}

#
tmp <- c(
sort(Group1_spot_rank), sort(Group2_spot_rank),
sort(Group3_spot_rank), sort(Group4_spot_rank),
sort(Group5_spot_rank), sort(Group6_spot_rank),
sort(Group7_spot_rank)
)
duplicate_spots_idx <- duplicated(tmp)
print(strrep("#", 15))
print(paste0("Duplicated samples: ", sum(duplicated(tmp))))
print(strrep("#", 15))

#
#
#
idx <- rownames(spots) %in% Group1_spot_rank
Group1_marker_spot <- spots[idx, ]
Group1_marker_spot$rank_prediction <- "caudate_putamen"

idx <- rownames(spots) %in% Group2_spot_rank
Group2_marker_spot <- spots[idx, ]
Group2_marker_spot$rank_prediction <- "corpus_callosum"

idx <- rownames(spots) %in% Group3_spot_rank
Group3_marker_spot <- spots[idx, ]
Group3_marker_spot$rank_prediction <- "cortex"

idx <- rownames(spots) %in% Group4_spot_rank
Group4_marker_spot <- spots[idx, ]
Group4_marker_spot$rank_prediction <- "globus_pallidus"

idx <- rownames(spots) %in% Group5_spot_rank
Group5_marker_spot <- spots[idx, ]
Group5_marker_spot$rank_prediction <- "hypothalamus"

idx <- rownames(spots) %in% Group6_spot_rank
Group6_marker_spot <- spots[idx, ]
Group6_marker_spot$rank_prediction <- "plexus"

idx <- rownames(spots) %in% Group7_spot_rank
Group7_marker_spot <- spots[idx, ]
Group7_marker_spot$rank_prediction <- "thalamus"

df <- do.call("rbind", list(Group1_marker_spot, Group2_marker_spot,
Group3_marker_spot, Group4_marker_spot, Group5_marker_spot,
Group6_marker_spot, Group7_marker_spot))
df <- df[!duplicate_spots_idx, ]



#
# plot
#
color_plot <- c(caudate_putamen = "#F0027F", corpus_callosum = "#377EB8",
cortex = "#4DAF4A", globus_pallidus = "#984EA3", hypothalamus = "#FFD700",
plexus = "#FF7F00", thalamus = "#1A1A1A", "NA" = "#666666", other = "#666666")

p_layers <- ggplot(spots, aes(x = imagecol, y = -imagerow, color = anatomy)) +
geom_point(size = rel(0.5)) +
theme_void() +
theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), legend.position = "none") +
scale_color_manual(values = color_plot)
#
#
#



color_plot <- c(color_plot, Right = "green", Wrong = "red")
names(color_plot) <- c("caudate_putamen", "corpus_callosum", "cortex", "globus_pallidus",
"hypothalamus", "plexus", "thalamus", "NA", "Other", "Correct", "Wrong")

colnames(df) = c("imagerow", "imagecol", "layer", "rank_prediction")

for (layer_plot in c("caudate_putamen", "corpus_callosum", "cortex", "globus_pallidus", "hypothalamus", "plexus", "thalamus")) {
    
    print("*")
    
    tmp_df <- df[!is.na(df$layer), ]
    tmp_df$Identification <- "Other"
    
    idx_tmp <- tmp_df$rank_prediction == layer_plot
    tmp_df$Identification[idx_tmp][tmp_df$layer[idx_tmp] == tmp_df$rank_prediction[idx_tmp]] = "_Right"
    tmp_df$Identification[idx_tmp][tmp_df$layer[idx_tmp] != tmp_df$rank_prediction[idx_tmp]] = "Wrong_"
    
    tmp_df <- tmp_df[order(tmp_df$layer), ]
    
    idx_right <- rownames(tmp_df)[tmp_df$Identification == "_Right"]
    idx_wrong <- rownames(tmp_df)[tmp_df$Identification == "Wrong_"]
    
    spots$Identification <- "Other"
    # spots$Identification[is.na(spots$anatomy)] <- "NA"
    spots$Identification[spots$anatomy == layer_plot] <- layer_plot
    spots$Identification[rownames(spots) %in% idx_right] <- "Correct"
    spots$Identification[rownames(spots) %in% idx_wrong] <- "Wrong"
    
    p <- ggplot(spots, aes(x = imagecol, y = -imagerow, color = Identification)) +
    geom_point(size = rel(1)) +
    theme_void() +
    theme(axis.text = element_blank(), legend.title = element_text(size = rel(4)), legend.text = element_text(size = rel(3.5))) +
    scale_color_manual(breaks = c("Correct", "Wrong"), values = color_plot) +
    guides(color = guide_legend(override.aes = list(size = rel(8))))
    
    pdf(paste0("/Users/haoran/Documents/MGGNN/figures/spots_identification/", sample, "_", layer_plot, "_threshold", ".pdf"), width = 14, height = 4)
    print(p_layers + p)
    dev.off()
}

#
tmp <- df[!is.na(df$layer), ]
acc <- sum(tmp$layer == tmp$rank_prediction)/nrow(tmp)
print(paste0("Sample: ", sample, ", ACC: ", acc))

spots <- spots[!is.na(spots$anatomy), ]
print(paste0("Percentage: ", dim(tmp)[1]/dim(spots)[1]))

write.csv(df, paste("/Users/haoran/Documents/MGGNN/data/spots_identification/", sample, "_threshold", ".csv", sep = ""), row.names = TRUE)

df_tmp <- data.frame(sample = sample, percentage = dim(tmp)[1]/dim(spots)[1], acc=acc)
if (file.exists("/Users/haoran/Documents/MGGNN/data/spots_identification/Heme_acc.csv")) {
    write.table(df_tmp, paste("/Users/haoran/Documents/MGGNN/data/spots_identification/Heme_acc.csv", sep = ""), row.names = FALSE, col.names = FALSE, append = TRUE)
} else {
    write.table(df_tmp, paste("/Users/haoran/Documents/MGGNN/data/spots_identification/Heme_acc.csv", sep = ""), row.names = FALSE)
}
