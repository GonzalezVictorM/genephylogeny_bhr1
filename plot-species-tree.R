# -----------------------------
# Clear environment and load libraries
# -----------------------------
rm(list = ls())
library(tidyverse)
library(ggtree)
library(treeio)
library(ape)
library(phytools)
library(ggtreeExtra)

# -----------------------------
# Color assignment
# -----------------------------

# Colorblind-safe palette
phylum_palette <- c(
  "Ascomycota" = "#E69F00",
  "Basidiomycota" = "#56B4E9",
  "Mucoromycota" = "#009E73",
  "Glomeromycota" = "#F0E442",
  "Chytridiomycota" = "#0072B2",
  "Blastocladiomycota" = "#D55E00",
  "Evosea" = "#000000",
  "Streptophyta" = "#000000"
)

lifestyle_palette <- c("saprotroph" = "#E69F00",
                       "WRF" = "#E69F00",
                       "BRF" = "#E69F00",
                       "plant" = "#56B4E9",
                       "plant_endophyte" = "#56B4E9",
                       "plant_mycorrhizal"= "#56B4E9",
                       #"plant_oppathogen",
                       "plant_pathogen"= "#56B4E9",
                       "animal"= "#009E73",
                       "animal_pathogen"= "#009E73",
                       #"animal_oppathogen",
                       #"human_pathogen",
                       #"human_oppathogen",
                       #"insect_pathogen",
                       "mycoparasite" = "#F0E442",
                       "entoparasite" = "#F0E442",
                       "symbiont" = "#D55E00")
lifestyle_sets <- c("saprotroph" = "saprotroph",
                    "WRF" = "saprotroph",
                    "BRF" = "saprotroph",
                    "plant" = "plant associated",
                    "plant_endophyte" = "plant associated",
                    "plant_mycorrhizal"= "plant associated",
                    #"plant_oppathogen",
                    "plant_pathogen"= "plant associated",
                    "animal"= "animal associated",
                    "animal_pathogen"= "animal associated",
                    #"animal_oppathogen",
                    #"human_pathogen",
                    #"human_oppathogen",
                    #"insect_pathogen",
                    "mycoparasite" = "ento/myco-associated",
                    "entoparasite" = "ento/myco-associated",
                    "symbiont" = "symbiont")

# -----------------------------
# Config paths
# -----------------------------
rooting_branch <- "Allma1"

input_file <- "species_tree"
data_dir <- "local_data"
phyl_dir <- file.path(data_dir, 'speciestree')
tree_dir <- file.path(phyl_dir, 'astral_results')

tree_file <- file.path(tree_dir, paste0(input_file, ".treefile"))

tax_file <- file.path(data_dir, "proteome_list_orthofinder_v2.csv")
tax_outgroup_file <- file.path(data_dir, "outgroup_phylogeny.csv")
output_tree_rect <- file.path(tree_dir, paste0(input_file, "_rect.pdf"))
output_tree_circ <- file.path(tree_dir, paste0(input_file, "_circ.pdf"))

# -----------------------------
# Load tree and root
# -----------------------------
stopifnot(file.exists(tree_file), 
          file.exists(tax_file), file.exists(tax_outgroup_file))

tree <- read.tree(tree_file)
stopifnot(!is.null(tree), length(tree$tip.label) > 0)

# Root by priority
if (is.na(rooting_branch)) {
  tree <- midpoint.root(tree)
} else {
  tree <- root(tree, outgroup = rooting_branch, resolve.root = TRUE)
}

# According to https://alexknyshov.github.io/R/page3.html modify root branches proportionally
tree$edge.length[which(!(tree$edge[,1] %in% tree$edge[,2]))] <- sum(tree$edge.length[which(!(tree$edge[,1] %in% tree$edge[,2]))])/2

# -----------------------------
# Load and process taxonomy and lifestyle metadata
# -----------------------------

# Load species metadata
tax_data <- read.csv(tax_file) %>%
  full_join(read.csv(tax_outgroup_file)) %>%
  select(portal, name, phylum:genus, primary_lifestyle) %>%
  mutate(portal = as.character(portal))

# Annotate taxonomy tip metadata
protein_metadata <- read.csv(tax_file) %>%
  full_join(read.csv(tax_outgroup_file)) %>%
  select(portal, name, phylum:genus, primary_lifestyle) %>%
  mutate(portal = as.character(portal),
         phylum = factor(phylum, levels = names(phylum_palette)),
         class = factor(class))

# Create binary matrix
lifestyle_matrix <- matrix(0, nrow = nrow(tax_data), ncol = length(names(lifestyle_palette)))
rownames(lifestyle_matrix) <- tax_data$portal
colnames(lifestyle_matrix) <- names(lifestyle_palette)

for (i in 1:nrow(tax_data)) {
  current_lifestyles <- trimws(unlist(str_split(tax_data$primary_lifestyle[i], ";\\s*")))
  lifestyle_matrix[i, names(lifestyle_palette) %in% current_lifestyles] <- 1
}

# Convert to df for ggtree
lifestyle_metadata <- data.frame(lifestyle_matrix) %>%
  rownames_to_column("portal") %>%
  pivot_longer(-portal,names_to = "lifestyle",values_to = "binary") %>%
  mutate(lifestyle_set = lifestyle_sets[lifestyle],
         lifestyle = factor(lifestyle, levels = names(lifestyle_palette)))

lifestyle_metadata$binary[is.na(lifestyle_metadata$binary)]<- 0

# -----------------------------
# Build the final tree
# -----------------------------
# The circular layout tree.
p <- ggtree(tree, size=0.15, open.angle=5) +
  #xlim(0,50) +
  geom_nodelab(size = 0.5, nudge_x = 0.01, nudge_y = 0.5)

p <- p %<+% protein_metadata +
  geom_tippoint(
    aes(color = phylum), 
    size = 0.4,
    show.legend = TRUE
  ) +
  geom_tiplab(
    aes(color = phylum),
    align=TRUE,
    linetype=3,
    size = 1.5,
    hjust = -0.1,
    linesize=0.2,
    show.legend=FALSE,
    offset = 0.15
  ) +
  scale_color_manual(
    name = "Phylum", 
    values = phylum_palette, 
    na.translate = FALSE
  )

p <- p +
  geom_fruit(
    data=lifestyle_metadata, 
    geom=geom_tile,
    mapping=aes(
      y=portal, 
      x=lifestyle, 
      alpha=binary, 
      fill = lifestyle_set),
    offset=0.08,
    axis.params = list(
      axis = "x",  # Display the x-axis
      text.size = 3, # Adjust tick label size
      text.angle = 90, # Rotate tick labels for better readability
      line.size = 0
    )
  ) +
  scale_fill_manual(
    name = "Lifestyle",
    values = unique(unname(lifestyle_palette))
  ) +
  geom_treescale(fontsize=2, linesize=0.3, x=4.9, y=0.1)

file_height <- 30
file_width <- 30

ggsave(output_tree_rect, plot = p, width = file_width, height = file_height, limitsize = FALSE)

# -----------------------------
# Build the final tree
# -----------------------------
# The circular layout tree.
p <- ggtree(tree, layout="fan", size=0.15, open.angle=5) +
  #xlim(0,50) +
  geom_nodelab(size = 0.5, nudge_x = 0.01, nudge_y = 0.5)

p <- p %<+% protein_metadata +
  geom_tippoint(
    aes(color = phylum), 
    size = 0.4,
    show.legend = TRUE
  ) +
  geom_tiplab(
    aes(color = phylum),
    align=TRUE,
    linetype=3,
    size = 1.5,
    hjust = -0.1,
    linesize=0.2,
    show.legend=FALSE,
    offset = 0.15
  ) +
  scale_color_manual(
    name = "Phylum", 
    values = phylum_palette, 
    na.translate = FALSE
  )

p <- p +
  geom_fruit(
    data=lifestyle_metadata, 
    geom=geom_tile,
    mapping=aes(
      y=portal, 
      x=lifestyle, 
      alpha=binary, 
      fill = lifestyle_set),
    offset=0.08,
    axis.params = list(
      axis = "x",  # Display the x-axis
      text.size = 3, # Adjust tick label size
      text.angle = 90, # Rotate tick labels for better readability
      line.size = 0
    )
  ) +
  scale_fill_manual(
    name = "Lifestyle",
    values = unique(unname(lifestyle_palette))
  ) +
  geom_treescale(fontsize=2, linesize=0.3, x=4.9, y=0.1)

file_height <- 30
file_width <- 30

ggsave(output_tree_circ, plot = p, width = file_width, height = file_height, limitsize = FALSE)

