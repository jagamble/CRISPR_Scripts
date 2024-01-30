# top_345_mle_comp.R
# Comparing output of MLE analysis from three separate datasets
# - Venn diagrams
# - volcano plots
# - split barplots

# 18 December 2023
# Edited version 30 January 2024

# Using MEDIAN normalized output
# Selecting hits by FDR only, as is no accepted cutoff for Beta score

library(tidyverse)
library( ggrepel )
library(RColorBrewer)
library(eulerr)

# Part_1: prepare MLE data from Top3, -4 and -5 datasets

# Read in MLE results file
top3_mle <- read_tsv("MLE_Output_Comp/01_Data_Files/top3_mnorm_mle.gene_summary.txt")
head(top3_mle)
top4_mle <- read_tsv("MLE_Output_Comp/01_Data_Files/top4_mnorm_mle.gene_summary.txt")
head(top4_mle)
top5_mle <- read_tsv("MLE_Output_Comp/01_Data_Files/top5_mnorm_mle.gene_summary.txt")
head(top5_mle)

# Select required columns and rename
top3_work <- top3_mle %>%
  select(Gene, contains('|beta'), contains('|fdr'))
head(top3_work)
colnames(top3_work) <- c('Gene', 'beta', 'fdr')
head(top3_work)

top4_work <- top4_mle %>%
  select(Gene, contains('|beta'), contains('|fdr'))
head(top4_work)
colnames(top4_work) <- c('Gene', 'beta', 'fdr')
head(top4_work)

top5_work <- top5_mle %>%
  select(Gene, contains('|beta'), contains('|fdr'))
head(top5_work)
colnames(top5_work) <- c('Gene', 'beta', 'fdr')
head(top5_work)

# Add 'dataset' column to each of above files
top3_work <- top3_work %>% mutate('dataset' = 'top3')
top4_work <- top4_work %>% mutate('dataset' = 'top4')
top5_work <- top5_work %>% mutate('dataset' = 'top5')
head(top3_work)
head(top4_work)
head(top5_work)

# Now join the above into one file
all_mle <- rbind(top3_work, top4_work)
head(all_mle)
all_mle <- rbind(all_mle, top5_work)
head(all_mle)
tail(all_mle)

# Add the 'hit' column
all_mle_hits <- all_mle %>%
  mutate('hit' = case_when(fdr < 0.1 & beta > 0 ~ 'UP',
                           fdr < 0.1 & beta < 0 ~ 'DOWN',
                           fdr < 0.1 & beta == 0 ~ 'NO_CHANGE',
                           TRUE ~ 'NOT'))
head(all_mle_hits)
tail(all_mle_hits)
write_tsv(all_mle_hits, "MLE_Output_Comp/02_Temp/all_mle_hits.txt")
# Gene|beta|fdr|dataset|hit

# Add scoring matrix for Venn diagram
all_mle_hits_scores <- all_mle_hits %>%
  mutate('UP' = case_when(hit == 'UP' ~ 1,
                          TRUE ~ 0),
         'DOWN' = case_when(hit == 'DOWN' ~ 1,
                            TRUE ~ 0))
head(all_mle_hits_scores)
write_tsv(all_mle_hits_scores, "MLE_Output_Comp/02_Temp/all_mle_hits_scores.txt")
# Gene|beta|fdr|dataset|hit|UP|DOWN



# Part_2a: prepare Venn diagram for 'UP' genes

# Extract required columns from previous file
mle_venn_UP <- all_mle_hits_scores %>%
  select(Gene, dataset, UP)
head(mle_venn_UP)
# Gene|dataset|UP

# Pivot wider by dataset
mle_venn_UP_wide <- mle_venn_UP %>%
  pivot_wider(names_from = dataset, values_from = UP)
head(mle_venn_UP_wide)
# Gene|top3|top4|top5

# Remove 'Gene' column and remove any line with 'NA' values
mle_venn_UP_wide_euler <- mle_venn_UP_wide %>%
  select(-Gene)
mle_venn_UP_wide_euler <- na.omit(mle_venn_UP_wide_euler)
head(mle_venn_UP_wide_euler)
# top3|top4|top5

# Create Venn diagram
mle_venn_UP_wide_euler <- euler(mle_venn_UP_wide_euler, shape = 'circle')
# plot(mle_venn_UP_wide_euler, quantities = TRUE, fills = brewer.pal(n = 12, name = "Set3"), alpha = 0.5)
# Increase size of text
plot(mle_venn_UP_wide_euler, quantities = list(cex = 1.5), labels = list(cex = 1.5),
     fills = brewer.pal(n = 12, name = "Set3"), alpha = 0.5)
# Save as venn_top345_mle_mnorm_UP

# Add 'hit_class' data for 'UP' genes
# Start with output from 'pivot_wider' above
# Could probably write this a bit more compactly - link by '%>%'
head(mle_venn_UP_wide)
mle_venn_UP_wide_analysis <- mle_venn_UP_wide %>%
  mutate(hit_count = top3 + top4 + top5)
head(mle_venn_UP_wide_analysis)
mle_venn_UP_wide_analysis <- mle_venn_UP_wide_analysis %>%
  arrange(desc(hit_count))
head(mle_venn_UP_wide_analysis)
mle_venn_UP_wide_analysis_hits <- mle_venn_UP_wide_analysis %>%
  filter(hit_count > 0)
  # Now reduced to 141 entries
# Dropping '_only' suffix in 'mutate()' step - helps with selecting unique hits in the volcano plot (see below)
mle_venn_UP_wide_analysis_hits <- mle_venn_UP_wide_analysis_hits %>%
  mutate(hit_class = case_when(hit_count == 3 ~ 'all',
                               hit_count == 1 & top3 == 1 ~ 'top3',
                               hit_count == 1 & top4 == 1 ~ 'top4',
                               hit_count == 1 & top5 == 1 ~ 'top5',
                               TRUE ~ 'other'))
head(mle_venn_UP_wide_analysis_hits)
write_tsv(mle_venn_UP_wide_analysis_hits, "MLE_Output_Comp/02_Temp/mle_venn_UP_wide_analysis_hits.txt")
# Gene|top3|top4|top5|hit_count|hit_class



# Part_2b: prepare Venn diagram for 'DOWN' genes

# Extract required columns from previous file
mle_venn_DOWN <- all_mle_hits_scores %>%
  select(Gene, dataset, DOWN)
head(mle_venn_DOWN)
# Gene|dataset|DOWN

# Pivot wider by dataset
mle_venn_DOWN_wide <- mle_venn_DOWN %>%
  pivot_wider(names_from = dataset, values_from = DOWN)
head(mle_venn_DOWN_wide)
# Gene|top3|top4|top5

# Remove 'Gene' column and remove any line with 'NA' values
mle_venn_DOWN_wide_euler <- mle_venn_DOWN_wide %>%
  select(-Gene)
mle_venn_DOWN_wide_euler <- na.omit(mle_venn_DOWN_wide_euler)
head(mle_venn_DOWN_wide_euler)
# top3|top4|top5

# Create Venn diagram
mle_venn_DOWN_wide_euler <- euler(mle_venn_DOWN_wide_euler, shape = 'circle')
# plot(mle_venn_DOWN_wide_euler, quantities = TRUE, fills = brewer.pal(n = 12, name = "Set3"), alpha = 0.5)
# Increase size of font
plot(mle_venn_DOWN_wide_euler, quantities = list(cex = 1.5), labels = list(cex = 1.5),
     fills = brewer.pal(n = 12, name = "Set3"), alpha = 0.5)
# Save as venn_top345_mle_mnorm_DOWN

# Add 'hit_class' data for 'DOWN' genes
# Start with output from 'pivot_wider' above
# Could probably write this a bit more compactly - link by '%>%'
head(mle_venn_DOWN_wide)
mle_venn_DOWN_wide_analysis <- mle_venn_DOWN_wide %>%
  mutate(hit_count = top3 + top4 + top5)
head(mle_venn_DOWN_wide_analysis)
mle_venn_DOWN_wide_analysis <- mle_venn_DOWN_wide_analysis %>%
  arrange(desc(hit_count))
head(mle_venn_DOWN_wide_analysis)
mle_venn_DOWN_wide_analysis_hits <- mle_venn_DOWN_wide_analysis %>%
  filter(hit_count > 0)
# Now reduced to 167 entries
# Dropped '_only' suffix in 'mutate()' step - helps with selecting unique hits in the volcano plot (see below)
mle_venn_DOWN_wide_analysis_hits <- mle_venn_DOWN_wide_analysis_hits %>%
  mutate(hit_class = case_when(hit_count == 3 ~ 'all',
                               hit_count == 1 & top3 == 1 ~ 'top3',
                               hit_count == 1 & top4 == 1 ~ 'top4',
                               hit_count == 1 & top5 == 1 ~ 'top5',
                               TRUE ~ 'other'))
head(mle_venn_DOWN_wide_analysis_hits)
write_tsv(mle_venn_DOWN_wide_analysis_hits, "MLE_Output_Comp/02_Temp//mle_venn_DOWN_wide_analysis_hits.txt")
# Gene|top3|top4|top5|hit_count|hit_class



# Part_3: create 'combo' file, with 'direction' column
# Starting files:
# - mle_venn_UP_wide_analysis_hits: Gene|top3|top4|top5|hit_count|hit_class
# - mle_venn_DOWN_wide_analysis_hits: Gene|top3|top4|top5|hit_count|hit_class

# Add a 'direction' column to each file as shown
mle_venn_UP_wide_analysis_hits_merge <- mle_venn_UP_wide_analysis_hits %>%
  mutate(direction = 'UP')
mle_venn_DOWN_wide_analysis_hits_merge <- mle_venn_DOWN_wide_analysis_hits %>%
  mutate(direction = 'DOWN')
head(mle_venn_UP_wide_analysis_hits_merge)
head(mle_venn_DOWN_wide_analysis_hits_merge)
# Both files: Gene|top3|top4|top5|hit_count|hit_class|direction

# Now combine both the above files, one on top of the other
mle_venn_COMBO_wide_analysis_hits <- rbind(mle_venn_UP_wide_analysis_hits_merge, mle_venn_DOWN_wide_analysis_hits_merge)
head(mle_venn_COMBO_wide_analysis_hits)
write_tsv(mle_venn_COMBO_wide_analysis_hits, "MLE_Output_Comp/02_Temp/mle_venn_COMBO_wide_analysis_hits.txt")
# Gene|top3|top4|top5|hit_count|hit_class|direction
# 308 entries



# Part_4: Create input file for Volcano plots
# Starting files:
# - all_mle: Gene|beta|fdr|dataset (53958 entries)
# - mle_venn_COMBO_wide_analysis_hits: Gene|top3|top4|top5|hit_count|hit_class|direction (308 entries)
# NB: In 'Gene' column of 'all_mle', entries are 3x, as file was created by merging three separate datasets
# - causes problems later, when trying to select unique genes - see Notes at end of script
# - solution here may not be the most elegant solution, but it works

head(all_mle)
head(mle_venn_COMBO_wide_analysis_hits)

# Select these columns
mle_venn_COMBO_hit_list <- mle_venn_COMBO_wide_analysis_hits %>%
  select(Gene, hit_class, direction)
head(mle_venn_COMBO_hit_list)
write_tsv(mle_venn_COMBO_hit_list, "MLE_Output_Comp/02_Temp/mle_venn_COMBO_hit_list.txt")
# Gene|hit_class|direction (308 entries)

# Join both files
all_mle_volcano <- left_join(all_mle, mle_venn_COMBO_hit_list)
head(all_mle_volcano)
# Beware! 'Gene' column entries are now 3x

# To identify hits *unique to* a specific dataset, create 'unique_hit' column as shown
# If 'hit_class' == 'dataset', entry is a unique hit, so enter value of '1'. Otherwise, enter '0'
# If 'unique_hit == 1', it is a hit that is unique to one or other of the original datasets
all_mle_volcano <- all_mle_volcano %>%
  mutate('unique_hit' = case_when(hit_class == dataset ~ 1,
                                  TRUE ~ 0))
head(all_mle_volcano)
write_tsv(all_mle_volcano, "MLE_Output_Comp/02_Temp/all_mle_volcano.txt")




# Volcano plot - showing hits UNIQUE to each dataset
MLE_volc <- ggplot(all_mle_volcano, aes(x= beta, y = -log10(fdr+1e-3))) + theme_classic() +
  geom_point(colour = "grey", alpha = 0.5, size = 0.5) +
  ggtitle("MLE hits unique to each dataset") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 20),
        strip.text = element_text(size = 14, face = "bold"),
        strip.background = element_rect(fill = "grey"))
MLE_volc
MLE_volc <- MLE_volc + geom_point(data = subset(all_mle_volcano, fdr < 0.1 & direction == 'UP'), colour = "red", size = 1.5)
MLE_volc
MLE_volc <- MLE_volc + geom_point(data = subset(all_mle_volcano, fdr < 0.1 & direction == 'DOWN'), colour = "blue3", size = 1.5)
MLE_volc
MLE_volc <- MLE_volc + geom_text_repel(data = subset(all_mle_volcano, direction == 'UP' & unique_hit == 1), aes(label = Gene ),
                                                       segment.color =  "red", force = 20, point.padding = 0.6, max.overlaps = 90, fontface = 'italic', colour = "red")
MLE_volc
MLE_volc <- MLE_volc + geom_text_repel(data = subset(all_mle_volcano, direction == 'DOWN' & unique_hit == 1), aes(label = Gene ),
                                                       segment.color =  "blue3", force = 20, point.padding = 0.6, max.overlaps = 90, fontface = 'italic', colour = "blue3")
MLE_volc
MLE_volc <- MLE_volc + ylab( '-log10(FDR+1e-3)' ) + xlab( 'MLE Beta score' ) + geom_hline(yintercept = -log10(0.1), linetype = 2, alpha = 0.5)
MLE_volc
MLE_volc <- MLE_volc + geom_vline(xintercept = c(beta = -0.5, beta =0.5), linetype = 2, alpha = 0.5) +
  facet_wrap(vars(dataset), scales = 'fixed')
MLE_volc
# Filename: volcano_top345_mle_UNIQUE



# Volcano plot - showing hits COMMON TO each dataset
MLE_volc <- ggplot(all_mle_volcano, aes(x= beta, y = -log10(fdr+1e-3))) + theme_classic() +
  geom_point(colour = "grey", alpha = 0.5, size = 0.5) +
  ggtitle("MLE hits common to each dataset") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 20),
        strip.text = element_text(size = 14, face = "bold"),
        strip.background = element_rect(fill = "grey"))
MLE_volc
MLE_volc <- MLE_volc + geom_point(data = subset(all_mle_volcano, fdr < 0.1 & direction == 'UP'), colour = "red", size = 1.5)
MLE_volc
MLE_volc <- MLE_volc + geom_point(data = subset(all_mle_volcano, fdr < 0.1 & direction == 'DOWN'), colour = "blue3", size = 1.5)
MLE_volc
MLE_volc <- MLE_volc + geom_text_repel(data = subset(all_mle_volcano, hit_class == 'all' & direction == 'UP'), aes(label = Gene ),
                                       segment.color =  "red", force = 20, point.padding = 0.6, max.overlaps = 90, fontface = 'italic', colour = "red")
MLE_volc
MLE_volc <- MLE_volc + geom_text_repel(data = subset(all_mle_volcano, hit_class == 'all' & direction == 'DOWN'), aes(label = Gene ),
                                       segment.color =  "blue3", force = 20, point.padding = 0.6, max.overlaps = 90, fontface = 'italic', colour = "blue3")
MLE_volc
MLE_volc <- MLE_volc + ylab( '-log10(FDR+1e-3)' ) + xlab( 'MLE Beta score' ) + geom_hline(yintercept = -log10(0.1), linetype = 2, alpha = 0.5)
MLE_volc
MLE_volc <- MLE_volc + geom_vline(xintercept = c(beta = -0.5, beta =0.5), linetype = 2, alpha = 0.5) +
  facet_wrap(vars(dataset), scales = 'fixed')
MLE_volc
# Filename: volcano_top345_mle_COMMON





# Part_5: Create split barplots - one plot for UNIQUE genes and one for GENES_in_COMMON
# 11 January 2024
# Starting file:  all_mle_volcano

# Filter out all entries w/ FDR < 0.1
# Note that it contains quite a few entries you don't want - filter these later
head(all_mle_volcano)
bar_data_all <- all_mle_volcano %>%
  filter(fdr < 0.1)
head(bar_data_all)
write_tsv(bar_data_all, "MLE_Output_Comp/02_Temp/bar_data_all.txt")
# Filters out 530 entries - 158 unique genes and 372 genes_in_common
# Gene|beta|fdr|dataset|hit_class|direction|unique_hit


# Create file of GENES_in_COMMON
bar_data_common <- bar_data_all %>%
  filter(hit_class == 'all')
head(bar_data_common)
# Gene|beta|fdr|dataset|hit_class|direction|unique_hit (216 entries - confirmed correct)

# Create new column for FDR 'bins' (0.1-0.01 and 0.01-0.001)
bar_data_common <- bar_data_common %>%
  mutate('hit_range' = case_when(fdr > 0.01 ~ 'low_hit',
                                 fdr < 0.01 ~ 'top_hit',
                                 TRUE ~ 'NOT'))
head(bar_data_common)
write_tsv(bar_data_common, "MLE_Output_Comp/02_Temp/bar_data_COMMON.txt")
# Gene|beta|fdr|dataset|hit_class|direction|unique_hit|hit_range

# Create plot for COMMON genes
ggplot(bar_data_common, aes(x = direction, fill = hit_range)) + theme_classic() +
  ggtitle("MLE hits common to Top3, 4 and 5") +
  theme(plot.title = element_text(size = 20),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14, face = "bold"),
        strip.background = element_rect(fill = "grey")) +
  scale_fill_brewer(palette = "Pastel1") +
  geom_bar(stat = "count", colour = "black", position = position_stack(reverse = TRUE)) +
  guides(fill = guide_legend(reverse = TRUE)) +
    facet_wrap(vars(dataset), scales = 'fixed')
# Save as 'barplot_top345_SHARED_hits'



# Create file of UNIQUE genes
head(bar_data_all)
bar_data_unique <- bar_data_all %>%
  filter(unique_hit == 1)
head(bar_data_unique)
# Gene|beta|fdr|dataset|hit_class|direction|unique_hit (158 entries - confirmed correct)

# Create new column for FDR 'bins' (0.1-0.01, and 0.01-0.001)
bar_data_unique <- bar_data_unique %>%
  mutate('hit_range' = case_when(fdr > 0.01 ~ 'low_hit',
                                 fdr < 0.01 ~ 'top_hit',
                                 TRUE ~ 'NOT'))
head(bar_data_unique)
write_tsv(bar_data_unique, "MLE_Output_Comp/02_Temp/bar_data_UNIQUE.txt")
# Gene|beta|fdr|dataset|hit_class|direction|unique_hit|hit_range

# Create plot for UNIQUE genes
ggplot(bar_data_unique, aes(x = direction, fill = hit_range)) + theme_classic() +
  ggtitle("MLE hits unique to Top3, 4 and 5") +
  theme(plot.title = element_text(size = 20),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14, face = "bold"),
        strip.background = element_rect(fill = "grey")) +
  scale_fill_brewer(palette = "Pastel1") +
  geom_bar(stat = "count", colour = "black", position = position_stack(reverse = TRUE)) +
  guides(fill = guide_legend(reverse = TRUE)) +
  facet_wrap(vars(dataset), scales = 'fixed')
# Save as 'barplot_top345_UNIQUE_hits'




# Notes on Script
# ---------------

# Problem with use of 'all_mle' file in Part_4 and how to label unique hits on volcano plot
# The entries in the 'Gene' column of 'all_MLE' are represented 3x, as the file was created by merging data from three datasets
# This created a 'one-to-many' problem when merging to create 'all_mle_volcano'
# Problem is that contents of 'hit_class' and 'direction' columns cannot now be trusted, as only one of them will be correct
# Solution as to:
# - remove the '_only' suffix from entries in the 'hit_class' column
# - this makes them same as the 'dataset' column (ie. 'top3' instead of 'top3_only)
# So... selection for unique hits can now just be IF 'hit_class' = 'dataset', THEN 'unique_hit = 1' (otherwise = '0')
# Makes things much simpler than previous code, which involved splitting 'hit_class' with 'separate()' to create the 'class' columns...etc
# Seems to work just as well

