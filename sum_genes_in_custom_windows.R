library(dplyr)

gene_number<-read_tsv("yourgff3.gff3")

##This code operates under the assumption that your gff3 file is named with 
#CHROM	type	pos_start	pos_end


#filter gff3 to just genes
gene_number<-gene_number %>%
  filter(type=="gene")
#add in a 1 next to each gene for counts
gene_number<-gene_number%>%
  mutate(gene_count = 1)

gene_number$midpoint<-((gene_number$pos_start + gene_number$pos_end)/2)
gene_number$gene_length<-(gene_number$pos_end - gene_number$pos_start)


sum_counts_in_range_by_chrom <- function(df_range, df_values) {
  # Initialize an empty data frame to store results
  result_df <- data.frame()
  
  # Get the list of unique chromosomes
  chromosomes <- unique(df_range$CHROM)
  
  # Loop through each chromosome
  for (chrom in chromosomes) {
    # Subset data frames by current chromosome
    subset_range <- df_range[df_range$CHROM == chrom, ]
    subset_values <- df_values[df_values$CHROM == chrom, ]
    
    # Sum counts within this chromosome range
    sums <- numeric(nrow(subset_range))
    for (i in 1:nrow(subset_range)) {
      subset_df <- subset_values[subset_values$midpoint >= subset_range$kb_start[i] & subset_values$midpoint <= subset_range$kb_end[i], ]
      sums[i] <- sum(subset_df$gene_count, na.rm = TRUE)
    }
    
    # Add the sums as a new column to subset_range
    subset_range$sum_gene_counts <- sums
    
    # Combine the result with the main result_df
    result_df <- rbind(result_df, subset_range)
  }
  
  # Return the result dataframe
  return(result_df)
}

gene_number_recomb<-sum_counts_in_range_by_chrom(your_data_frame_with_custom_windows, gene_number)
