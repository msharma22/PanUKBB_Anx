library(dplyr)

# Example data frame

# Get unique genes
background_genes <- META_twas_sumstats %>% distinct(gene_name) 


write.csv(background_genes, file = "background_genes.csv", row.names = FALSE)


# Print the list of unique genes
print(unique_genes)



########################################################################## P-value threshold df ##########################################################################
# Load the dplyr package
library(dplyr)

# Define the threshold p-value
threshold <- 3.85e-6

# Filter the data frame using dplyr's filter function
df_filtered <- META_twas_sumstats %>% dplyr::filter(UK_pvalue < threshold)

# Display the filtered data frame
print(df_filtered)

# Assuming your data frame is META_twas_sumstats
# and it contains a column called 'gene'

# Create a new data frame with unique genes
unique_genes_df <- df_filtered %>%
  distinct(gene_name, .keep_all = TRUE)

unique_gene_names <- unique_genes_df %>%
  distinct(gene_name)


# View the result
print(unique_gene_names)

# Save the data frame 'df_filtered' to a CSV file
write.csv(unique_gene_names, file = "filtered_genes.csv", row.names = FALSE)


########################################################################## FDR threshold df ##########################################################################
# Load the dplyr package
library(dplyr)

# Define the threshold p-value
threshold <- 0.05

# Filter the data frame using dplyr's filter function
df_filtered_fdr <- META_twas_sumstats %>% dplyr::filter(UK_pvalue_fdr < threshold)

 
# Assuming your data frame is META_twas_sumstats
# and it contains a column called 'gene'

# Create a new data frame with unique genes
unique_fdr_genes_ <- df_filtered_fdr %>%
  distinct(gene_name, .keep_all = TRUE)

unique_fdr_gene_names <- unique_fdr_genes_ %>%
  distinct(gene_name)

# Save the data frame 'df_filtered' to a CSV file
write.csv(unique_fdr_gene_names, file = "filtered_fdr_genes.csv", row.names = FALSE)


