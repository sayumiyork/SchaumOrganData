#' plotAcrossKidney
#'
#' Creates a figure plotting the gene expression using filtered TD-data from the Schaum et al. (2020) kidney dataset across samples from both sexes of mice over time
#'
#'
#' @param x A string containing a valid Mouse Genome Informatics (MGI) ID that is found in the dataset
#' @return A figure plotting the gene expression of filtered TD-data from the Schaum et al. (2020) kidney dataset across samples from both sexes of mice over time
#' @export


plotAcrossKidney<- function( x ) {
  df <- data.frame(counts(kidney, normalized=TRUE)[x,], kidney_metadata$age_sex, row.names = row.names(kidney_metadata))
  colnames(df) <- c("counts", "age_sex")
  df$age_sex <- fct_relevel( factor( df$age_sex), "1_m","3_m","6_m","9_m","12_m", "15_m", "18_m", "21_m", "24_m", "27_m", "1_f","3_f","6_f","9_f","12_f", "15_f", "18_f", "21_f")
  ggplot( df ) +
    geom_bar( aes( age_sex, counts ), stat="identity" ) +
    ggtitle(x)+
    xlab("Mice (Kidney)") +
    ylab("Normalized counts") +
    theme_bw() +
    theme( axis.text.x = element_text( angle=90, vjust=0.5 ) )
}
