rm(list = ls())

## You'll need to install poolr to use the PicMin function
library(poolr)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(PicMin)



## Perform Picmin on the Arabidopsis data from Bohutinska et al


## quick function for getting names

get_names <- function( lineage_df ){
  return( paste( lineage_df$scaff, lineage_df$start,
                 sep = "_") )
}


## quick function for minimising the dataframes

min_lin <- function( lineage_df , lineage_name){
  tmp <- data.frame( emp_p = lineage_df$emp_p,
                     window = lineage_df$name)
  names(tmp) <- c(lineage_name,
                  "window")
  return( tmp )
}


## Read in the data, calculate empirical p-values and use a common naming scheme for all lineages

lin_1 <- read.csv("~/UBC/GEA/pMax/Arabidopsis/Analysis_for_Magdalena/picmin/DIATEA_WS1000_MS10_BPM.txt",
                  sep = "\t")
# Calculate empirical p-values from Hudson's Fst
lin_1$emp_p <- PicMin:::EmpiricalPs(lin_1$FstH, large_i_small_p = T)
# Get locus names
lin_1$name <- get_names( lin_1 )
# Make a lightweight version of the dataframe
lin_1_m <- min_lin( lin_1, "DIATEA" )

# Repeat the above for the other 3 lineages

lin_2 <- read.csv("~/UBC/GEA/pMax/Arabidopsis/Analysis_for_Magdalena/picmin/OSLTEM_WS1000_MS10_BPM.txt",
                  sep = "\t")
lin_2$emp_p <- PicMin:::EmpiricalPs(lin_2$FstH, large_i_small_p = T)

lin_2$name <- get_names( lin_2 )
lin_2_m <- min_lin( lin_2, "OSLTEM" )

lin_3 <- read.csv("~/UBC/GEA/pMax/Arabidopsis/Analysis_for_Magdalena/picmin/STDSCT_WS1000_MS10_BPM.txt",
                  sep = "\t")
lin_3$emp_p <- PicMin:::EmpiricalPs(lin_3$FstH, large_i_small_p = T)
lin_3$name <- get_names( lin_3 )
lin_3_m <- min_lin( lin_3, "STDSCT" )

lin_4 <- read.csv("~/UBC/GEA/pMax/Arabidopsis/Analysis_for_Magdalena/picmin/VLHMOD_WS1000_MS10_BPM.txt",
                  sep = "\t")
lin_4$emp_p <- PicMin:::EmpiricalPs(lin_4$FstH, large_i_small_p = T)
lin_4$name <- get_names( lin_4 )
lin_4_m <- min_lin( lin_4, "VLHMOD" )

## Merge dataframes into 1


#put all data frames into list
df_list <- list(lin_1_m,
                lin_2_m,
                lin_3_m,
                lin_4_m)

#merge all data frames in list
all_lins <- df_list %>% reduce(full_join, by='window')

# Make a dataframe with just the p-values
all_lins_p <- all_lins[ , !(names(all_lins) %in% c("window"))]

# Use locus names as row names
rownames(all_lins_p) <- all_lins$window


# Set the number of lineages
nLins <- 4

# A count variable
count = 0

# An empty list to store the results
results = list()

# Give the number of lineags present that you will analyse
lins_present <- c(3,4)

# NOTE: Adding 2 to the list will likely crash the script
for (n in lins_present){
  count = count + 1
  # Run 400,000 replicate simulations of this situation and build the correlation matrix
  emp_p_null_dat <- t(replicate(400000, PicMin:::GenerateNullData(1.0, a=0.5, b=3, n=n, genes=100000)))
  # Calculate the order statistics p-values for each simulation
  emp_p_null_dat_unscaled <- t(apply(emp_p_null_dat ,1, PicMin:::orderStatsPValues))
  # Use those p-values to construct the correlation matrix
  null_pMax_cor_unscaled <- cor( emp_p_null_dat_unscaled )
  
  # Screen out genes that do not have the appropriate level of missing data
  lins_p_n_screened <-  as.matrix(all_lins_p[rowSums(is.na(all_lins_p)) == nLins-n,])
  
  # Skip to next loop iteration if the above conditions are not met
  if (dim(lins_p_n_screened)[1] ==0){
    next
  }
  # Make containers for the results
  res_p <- rep(-1,
               nrow(lins_p_n_screened))
  res_n <- rep(-1,
               nrow(lins_p_n_screened))
  # For each locus, apply PicMin - Could be parallelised if you were keen
  for (i in seq(nrow(lins_p_n_screened)) ){
    test_result <- PicMin:::PicMin(na.omit(lins_p_n_screened[i,]), null_pMax_cor_unscaled, numReps = 1000000)
    res_p[i] <- test_result$p
    res_n[i] <- test_result$config_est
  }
  # Store the results as a dataframe and add to the overall list
  results[[count]] = data.frame(numLin = n ,
                                p = res_p,
                                q = p.adjust(res_p, method = "fdr"),
                                n_est = res_n,
                                locus = row.names(lins_p_n_screened) )
  
}

# Merge the dataframes
picMin_results <- do.call(rbind, results)

# If there are numerous entries with equal values of p,
# then it could indicate that the number of samples in PoolR was not high enough.
# This is not a huge deal, but can be rectified by tweaking the numReps parameter above
picMin_results[picMin_results$p==min(picMin_results$p),]

# Calculate a single FDR corrected p-value for all analyses performed
picMin_results$pooled_q <- p.adjust(picMin_results$p, method = "fdr")

# How many genome-wide hits are there?
nrow(picMin_results[picMin_results$pooled_q<0.05,])

# Construct the final dataframe
picMin_results <- cbind( picMin_results,
                         read.csv(text=picMin_results$locus, header=FALSE,
                                  sep = "_",
                                  col.names=c('redundan','scaffold','start'))
)

library(ggplot2)
# A colour pallette I happen to like
col_pal <- c("white", "#8ec641", "#897696", "#e93826", "#13a4f5", "#f89b56")

# Drop the null scaffold data:
picMin_results <- picMin_results[picMin_results$scaffold!=-9,]

# Convert the scaffold variable to a factor
picMin_results$scaffold <- factor(picMin_results$scaffold,
                                  labels = paste("Scaffold",1:8))


# Make a Manhattan plot:

manhattanPlot <- ggplot(data = picMin_results,
                        aes(x = start/1e6,
                            y = -log10(pooled_q),
                            fill = factor(n_est)))+
  geom_point(shape = 21,
             size = 4)+
  geom_hline(aes(yintercept = -log10(0.05)),
             lty=2)+
  facet_wrap(~scaffold,
             ncol = 4,
             scales = "free_x")+
  scale_fill_manual(expression(italic(n[est])),values = col_pal)+
  scale_y_continuous(expression(-log[10]*"("*italic("q")*"-value)"))+
  scale_x_continuous("Position in Scaffold (Mbp)")+
  theme_half_open() +
  theme(strip.background = element_blank())+
  background_grid()

# Save the results as a CSV
write.csv(picMin_results,
          file = "~/UBC/GEA/pMax/Arabidopsis/Analysis_for_Magdalena/PicMinResults.csv", row.names = F)
# Save the Plot as a PNG
ggsave(manhattanPlot,
       file = "~/UBC/GEA/pMax/Arabidopsis/Analysis_for_Magdalena/genomeDup_arabidopsisPicMin.png",
       width = 9.50,
       height = 5.00)



