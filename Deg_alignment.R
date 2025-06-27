# Load required libraries
library(dplyr)

# Define file paths
file_paths <- c(
  "C:/Users/hebertlab-admin/AppData/Local/Temp/0b4644ef-c81f-4a26-9ef7-56b74de572d5_JEM_20180247_TablesS1andS2 (1).zip.2d5/Irradiation.csv",
  "C:/Users/hebertlab-admin/AppData/Local/Temp/0b4644ef-c81f-4a26-9ef7-56b74de572d5_JEM_20180247_TablesS1andS2 (1).zip.2d5/B_PLX_LPS.csv",
  "C:/Users/hebertlab-admin/AppData/Local/Temp/0b4644ef-c81f-4a26-9ef7-56b74de572d5_JEM_20180247_TablesS1andS2 (1).zip.2d5/B_PLX_Saline.csv"
)
# Read CSV files into a list of data frames
data_frames <- lapply(file_paths, read.csv, row.names = 1)

# Assuming data_frames contains three data frames from your CSV files
# Extract each data frame and convert to a matrix

# Convert the first data frame to a matrix
irradiation_matrix <- as.matrix(data_frames[[1]])

# Convert the second data frame to a matrix
b_plx_lps_matrix <- as.matrix(data_frames[[2]])

# Convert the third data frame to a matrix
b_plx_saline_matrix <- as.matrix(data_frames[[3]])

# Check the first few rows of each matrix to verify
print(head(irradiation_matrix))
print(head(b_plx_lps_matrix))
print(head(b_plx_saline_matrix))

# Load necessary library
library(biomaRt)

# Initialize the biomaRt for mouse genes
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Function to convert gene symbols to Ensembl IDs
convert_symbols_to_ensembl_ids <- function(data_frames) {
  # Ensure there's a 'Symbol' column
  if (!("Symbol" %in% colnames(data_frames))) {
    stop("The data frame does not have a 'Symbol' column.")
  }
  
  # Retrieve Ensembl IDs
  results <- getBM(
    attributes = c("mgi_symbol", "ensembl_gene_id"),
    filters = "mgi_symbol",
    values = data_frames$Symbol,
    mart = ensembl
  )
  
  # Merge results back using 'Symbol' as key
  merged_data_frames <- merge(data_frames, results, by.x = "Symbol", by.y = "mgi_symbol", all.x = TRUE)
  
  return(merged_data_frames)
}

# Clean up symbols before conversion (e.g., trimming whitespace)
clean_symbols <- function(symbols) {
  symbols <- gsub("[[:space:]]", "", symbols)  # Remove white spaces
  return(symbols)
}

# Apply cleanup and conversion to each data frame
converted_data_frames <- lapply(data_frames, function(data_frames) {
  data_frames$Symbol <- clean_symbols(data_frames$Symbol)  # Clean symbols if necessary
  convert_symbols_to_ensembl_ids(data_frames)  # Convert symbols
})

# Check the updated first data frame
print(head(converted_data_frames[[1]]))
