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
datasets <- c("OG0000194","Char_septins")
iterations <- c("","_clean","_clean2","_clean3")
trimmings <- c("","_80trim")
outgroups <- c("","_GtpA","_Myo2","_hMyo2")
rooting_branches <- c("","Dicdi-GtpA","Sacce1-Myo2","Homsa-Myo2")
names(rooting_branches) <- outgroups

dataset <- datasets[1]
iteration <- iterations[4]
trimming <- trimmings[2]
outgroup <- outgroups[3]
rooting_branch <- rooting_branches[outgroup]

input_file <- paste0(dataset,iteration,outgroup,"_mafft",trimming)
data_dir <- "local_data"
phyl_dir <- file.path(data_dir, 'phylogeny_analysis')
tree_dir <- file.path(phyl_dir, 'iqtree_files')

tree_file <- file.path(tree_dir, paste0(input_file, ".treefile"))

char_file <- file.path(data_dir, "char_proteins.csv")
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
# Load and process characterized proteins metadata
# -----------------------------
char_proteins <- read.csv(char_file)
stopifnot(all(c("protein_id", "new_id") %in% names(char_proteins)))

# Rename tips
label_map <- setNames(tree$tip.label, tree$tip.label)
matches <- char_proteins$protein_id %in% tree$tip.label
label_map[char_proteins$protein_id[matches]] <- char_proteins$new_id[matches]
tree$tip.label <- unname(label_map[tree$tip.label])

# -----------------------------
# Load and process taxonomy and lifestyle metadata
# -----------------------------
# Load species metadata
tax_data <- read.csv(tax_file) %>%
  full_join(read.csv(tax_outgroup_file)) %>%
  select(portal, name, phylum:genus, primary_lifestyle) %>%
  mutate(portal = as.character(portal),
         phylum = factor(phylum, levels = names(phylum_palette)),
         class = factor(class))

# Annotate taxonomy tip metadata
protein_metadata <- data.frame(protein_id = tree$tip.label) %>%
  mutate(portal = str_split_i(protein_id, "-", 1)) %>%
  left_join(tax_data, by = "portal") %>%
  select(!portal) %>%
  mutate(font_style = "plain",
         font_size = 1.5)

protein_metadata$font_style[protein_metadata$new_id %in% char_proteins$protein_id] <- "bold"
protein_metadata$font_size[protein_metadata$new_id %in% char_proteins$protein_id] <- 2

# Create binary matrix
lifestyle_matrix <- matrix(0, nrow = nrow(tax_data), ncol = length(names(lifestyle_palette)))
rownames(lifestyle_matrix) <- tax_data$portal
colnames(lifestyle_matrix) <- names(lifestyle_palette)

for (i in 1:nrow(tax_data)) {
  current_lifestyles <- trimws(unlist(str_split(tax_data$primary_lifestyle[i], ";\\s*")))
  lifestyle_matrix[i, names(lifestyle_palette) %in% current_lifestyles] <- 1
}

# Convert to df for ggtree
lifestyle_df <- data.frame(lifestyle_matrix) %>%
  rownames_to_column("portal") %>%
  pivot_longer(-portal,names_to = "lifestyle",values_to = "binary") %>%
  mutate(lifestyle_set = lifestyle_sets[lifestyle],
         lifestyle = factor(lifestyle, levels = names(lifestyle_palette)))

lifestyle_df$binary[is.na(lifestyle_df$binary)]<- 0

lifestyle_metadata <- data.frame(protein_id = tree$tip.label) %>%
  mutate(portal = str_split_i(protein_id, "-", 1)) %>%
  left_join(lifestyle_df, by = "portal", relationship ="many-to-many") %>%
  select(!portal)

# -----------------------------
# Define the nodes to highlight
# -----------------------------

highlight_nodes <- data.frame(
  node = c(1061,1208,908,759,756,1225,1271),
  type = "node")

# -----------------------------
# Identify monophyletic clades: by class
# -----------------------------
annotation_set <- "class"
all_nodes <- 1:max(tree$edge)
clade_labels <- list()

for (node in all_nodes) {
  desc <- getDescendants(tree, node)
  tips <- desc[desc <= length(tree$tip.label)]
  labels <- tree$tip.label[tips]
  vals <- protein_metadata %>%
    filter(protein_id %in% labels) %>% 
    pull(!!sym(annotation_set)) %>% 
    unique()
  if (length(vals) == 1 && !is.na(vals)) {
    clade_labels[[length(clade_labels) + 1]] <- list(node = node, label = vals)
  }
}

# Remove the descendant clades
filtered_set <- list()
for (i in seq_along(clade_labels)) {
  cur <- clade_labels[[i]]
  nested <- any(sapply(clade_labels, function(other) {
    cur$label == other$label &&
      cur$node != other$node &&
      cur$node %in% getDescendants(tree, other$node)
  }))
  if (!nested) filtered_set[[length(filtered_set) + 1]] <- cur
}

set_clades <- as.data.frame(do.call(rbind, filtered_set))
set_clades <- data.frame(node = unlist(set_clades$node),
                         name = unlist(set_clades$label)) %>%
  mutate(offset.text = seq(0.05, 0.3, length.out = n()))

# -----------------------------
# Build the final trees
# -----------------------------
# The rectangular layout tree.
p <- ggtree(tree, 
            size=0.15
) +
  xlim(-0.1,NA) +
  ylim(-0.1,NA) +
  geom_nodelab(size = 2, nudge_x = 0.01, nudge_y = 0.5) +
  geom_highlight(data = highlight_nodes,aes(node=node))

p <- p %<+% protein_metadata +
  geom_tippoint(
    aes(color = phylum), 
    size = 0.4,
    show.legend = TRUE
  ) +
  geom_tiplab(
    aes(color = phylum, fontface = font_style),
    align=TRUE,
    linetype=3,
    size = 3,
    hjust = -0.1,
    linesize= 0.2,
    show.legend=FALSE,
    offset = 0.05
  ) +
  scale_color_manual(
    name = "Phylum", 
    values = phylum_palette, 
    na.translate = FALSE
  )

p <- p +
  geom_cladelab(
    data = set_clades,
    mapping = aes(
      node = node,
      label = name
      # offset.text = offset.text,
    ),
    align = FALSE,
    offset = 10,
    #offset.text = 0.24,
    hjust = "left",
    barsize = 1,
    fontsize = 3,
    #angle = "auto",
    horizontal = FALSE
  )

p <- p +
  geom_fruit(
    data=lifestyle_metadata, 
    geom=geom_tile,
    mapping=aes(
      y=protein_id, 
      x=lifestyle, 
      alpha=binary, 
      fill = lifestyle_set,),
    offset = 0.34,
    # width = 0.1,
    # height = 1,
    axis.params = list(
      axis = "x",  # Display the x-axis
      text.size = 3, # Adjust tick label size
      text.angle = 30, # Rotate tick labels for better readability
      line.size = 0,
      hjust = "right"
    )
  ) +
  scale_fill_manual(
    name = "Lifestyle",
    values = unique(unname(lifestyle_palette))
  ) +
  geom_treescale(fontsize=3, linesize=1, x=0.1, y=-0.05) 

file_height <- 40
file_width <- 30

ggsave(output_tree_rect, 
       plot = p, 
       width = file_width, 
       height = file_height, 
       limitsize = FALSE)

# -----------------------------
# Build the final tree
# -----------------------------
# The circular layout tree.
p <- ggtree(tree, 
            layout="fan", 
            size=0.15, 
            open.angle=10
) +
  xlim(-0.1,NA) +
  #ylim(-0.1,NA) +
  geom_nodelab(size = 3, nudge_x = -0.1, nudge_y = -0.8) +
 #geom_nodelab(aes(label = node),size = 3, nudge_x = 0.01, nudge_y = 0.5) +
  geom_highlight(data = highlight_nodes,aes(node=node))

p <- rotate_tree(p, -90)

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
    size = 4,
    hjust = -0.1,
    linesize= 0.2,
    show.legend=FALSE,
    offset = 0.1
  ) +
  scale_color_manual(
    name = "Phylum", 
    values = phylum_palette, 
    na.translate = FALSE
  )

p <- p +
  geom_cladelab(
    data = set_clades,
    mapping = aes(
      node = node,
      label = name,
      offset.text = offset.text,
    ),
    align = TRUE,
    offset = 1.2,
    #offset.text = 0.24,
    hjust = "center",
    barsize = 1,
    fontsize = 4,
    angle = "auto",
    horizontal = FALSE
  )

p <- p +
  # geom_fruit(
  #   data=lifestyle_metadata, 
  #   geom=geom_tile,
  #   mapping=aes(
  #     y=protein_id, 
  #     x=lifestyle, 
  #     alpha=binary, 
  #     fill = lifestyle_set,),
  #   offset = 0.4,
  #   pwidth = 0.4,
  #   # height = 1,
  #   axis.params = list(
  #     axis = "x",  # Display the x-axis
  #     text.size = 3, # Adjust tick label size
  #     #text.angle = 90, # Rotate tick labels for better readability
  #     line.size = 0,
  #     hjust = "right"
  #   )
  # ) +
  # scale_fill_manual(
  #   name = "Lifestyle",
  #   values = unique(unname(lifestyle_palette))
  # ) +
  geom_treescale(fontsize=3, linesize=1, x=0.5, y=0.1) 

file_height <- 50
file_width <- 50

ggsave(output_tree_circ, 
       plot = p, 
       width = file_width, 
       height = file_height, 
       limitsize = FALSE)