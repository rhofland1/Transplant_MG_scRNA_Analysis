file_paths <- c(
  "C:/Users/hebertlab-admin/AppData/Local/Temp/0b4644ef-c81f-4a26-9ef7-56b74de572d5_JEM_20180247_TablesS1andS2 (1).zip.2d5/Irradiation.csv",
  "C:/Users/hebertlab-admin/AppData/Local/Temp/0b4644ef-c81f-4a26-9ef7-56b74de572d5_JEM_20180247_TablesS1andS2 (1).zip.2d5/B_PLX_LPS.csv",
  "C:/Users/hebertlab-admin/AppData/Local/Temp/0b4644ef-c81f-4a26-9ef7-56b74de572d5_JEM_20180247_TablesS1andS2 (1).zip.2d5/B_PLX_Saline.csv",
  "C:/Users/hebertlab-admin/AppData/Local/Temp/082b2473-3d2c-4fe4-ab0f-bad070e3cdfd_JEM_20180247_TablesS1andS2 (1).zip.dfd/Cx3cr1_Cre.csv",
  "D:/Rohan/Cluster_Data/DAM_DEG_Check/Host.csv",
  "D:/Rohan/Cluster_Data/DAM_DEG_Check/Transplant.csv"
)

# Read CSV files into a list of data frames
data_frames <- lapply(file_paths, read.csv, row.names = 1)

# Convert the first data frame to a matrix
irradiation_matrix <- as.matrix(data_frames[[1]])

# Convert the second data frame to a matrix
b_plx_lps_matrix <- as.matrix(data_frames[[2]])

# Convert the third data frame to a matrix
b_plx_saline_matrix <- as.matrix(data_frames[[3]])

cre_matrix <- as.matrix(data_frames[[4]])

Host_matrix <- as.matrix(data_frames[[5]])

TMG_matrix <- as.matrix(data_frames[[6]])

df1 <- as.data.frame(irradiation_matrix)
df2 <- as.data.frame(b_plx_lps_matrix)
df3 <- as.data.frame(b_plx_saline_matrix)
df4 <- as.data.frame(cre_matrix)
df5 <- as.data.frame(Host_matrix)
df6 <- as.data.frame(TMG_matrix)

df1$EnsemblID <- rownames(df1)
df2$EnsemblID <- rownames(df2)
df3$EnsemblID <- rownames(df3)
df4$EnsemblID <- rownames(df4)
# Merge the data frames by 'EnsemblID' column, keeping only rows present in all dataframes
merged_df <- merge(df1, df2, by = "EnsemblID", all = FALSE)
merged_df <- merge(merged_df, df3, by = "EnsemblID", all = FALSE)
merged_df <- merge(merged_df, df4, by = "EnsemblID", all = FALSE)

# Set 'EnsemblID' back as row names and remove the column
rownames(merged_df) <- merged_df$EnsemblID
merged_df$EnsemblID <- NULL

# Convert the final data frame to a matrix
merged_df <- as.data.frame(merged_df)
cronk_matrix <- as.matrix(merged_df)

# Check the combined matrix
print(head(combined_counts_df))

# Find common row names in both matrices
common_rownames <- intersect(rownames(df5), rownames(df6))

# Subset both matrices to retain only the common row names
aligned_df5 <- Host_matrix[common_rownames, , drop = FALSE]
aligned_df6 <- TMG_matrix[common_rownames, , drop = FALSE]

# Merge the matrices by binding columns side-by-side
combined_matrix <- cbind(aligned_df5, aligned_df6)

combined_df <- as.data.frame(combined_matrix)

# Find the common row names
common_rownames <- intersect(rownames(combined_df), merged_df$Symbol.x)

# Subset both matrices to these common rows
combined_hof_aligned <- combined_matrix[common_rownames, , drop = FALSE]
combined_cronk_aligned <- cronk_matrix[common_rownames, , drop = FALSE]

# Ensure column names are unique before merging to prevent conflicts
colnames(combined_hof_aligned) <- paste0("old_", colnames(combined_hof_aligned))
colnames(combined_cronk_aligned) <- paste0("new_", colnames(combined_cronk_aligned))

# Combine the two matrices into one
final_combined_counts <- cbind(combined_hof_aligned, combined_cronk_aligned)

# Check the resulting combined matrix
print(head(final_combined_counts))

write.csv(final_combined_counts, "D:/Rohan/Counts_Finalv2.csv")