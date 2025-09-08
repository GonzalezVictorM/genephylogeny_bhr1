# -----------------------------
# Clear environment and load libraries
# -----------------------------
rm(list = ls())
library(tidyverse)
library(ggtree)
library(treeio)
library(ape)
library(phytools)

# -----------------------------
# Config paths
# -----------------------------
datasets <- c("OG0000195","Char_septins")
iterations <- c("","_clean","_clean2")
trimmings <- c("","_80trim")
outgroups <- c("","_GtpA","_Myo2","_hMyo2")
rooting_branches <- c("","Dicdi-GtpA","Sacce1-Myo2","Homsa-Myo2")
names(rooting_branches) <- outgroups


dataset <- datasets[2]
iteration <- iterations[1]
trimming <- trimmings[2]
outgroup <- outgroups[1]
rooting_branch <- rooting_branches[outgroup]
file_height <- 20
file_width <- 20

input_file <- paste0(dataset,iteration,outgroup,"_mafft",trimming)
data_dir <- "local_data"
phyl_dir <- file.path(data_dir, 'phylogeny_analysis')
tree_dir <- file.path(phyl_dir, 'iqtree_files')

tree_file <- file.path(tree_dir, paste0(input_file, ".treefile"))

septin_file <- file.path(phyl_dir, "Septins.csv")
tax_file <- file.path(data_dir, "proteome_list_orthofinder_v2.csv")
tax_outgroup_file <- file.path(phyl_dir, "outgroup_phylogeny.csv")
output_tree <- file.path(tree_dir, paste0(input_file, ".pdf"))
output_node_tree <- file.path(tree_dir, paste0(input_file, "_nodes.pdf"))

# -----------------------------
# Load tree and root
# -----------------------------
stopifnot(file.exists(tree_file), file.exists(septin_file), 
          file.exists(tax_file), file.exists(tax_outgroup_file))

tree <- read.tree(tree_file)
stopifnot(!is.null(tree), length(tree$tip.label) > 0)

# -----------------------------
# Load and process metadata
# -----------------------------
# Load characterized septins
char_septins <- read.csv(septin_file)
stopifnot(all(c("protein_id", "new_id") %in% names(char_septins)))

# Rename tips
label_map <- setNames(tree$tip.label, tree$tip.label)
matches <- char_septins$protein_id %in% tree$tip.label
label_map[char_septins$protein_id[matches]] <- char_septins$new_id[matches]
tree$tip.label <- unname(label_map[tree$tip.label])

# -----------------------------
# Plot
# -----------------------------

p <- ggtree(tree, layout="daylight") +
  geom_tiplab(size = 4, hjust = -0.1) +
  geom_nodelab()

ggsave("Figure2_tree.pdf", plot = p, width = file_width, height = file_height, limitsize = FALSE)
 