####a seperate file to store functions 

#packages
install.packages("phytools")
library('phytools')
install.packages("devtools")
install.packages('treebalance')
library('treebalance')
library('ape')
install.packages("phyloTop")
library(phyloTop)
library(tidyverse)
install.packages("phangorn")
library(phangorn)
install.packages("sjPlot")
library(sjPlot)

install.packages('PDcalc')
devtools::install_github("davidnipperess/PDcalc",build_vignettes = F )
library(devtools)
install.packages(c("dplyr", "tidyr", "rmarkdown", "knitr"))
library(PDcalc)
install.packages("dispRity")
library(dispRity)
install.packages("TreeTools")
library(TreeTools)
install.packages("ggplot2")

# Bioconductor version
library(BiocManager)
BiocManager::install("YuLab-SMU/treedataverse")

# github version
devtools::install_github("YuLab-SMU/ggtree")
install("ggtree",force = TRUE)
library(ggplot2)
library(ggtree)
library(treeio)

install.packages("project/ggtree-devel/R/ggtree.R", repos = NULL, type = "source")
library(ggtree)

mtree$edge



#1. imbalance of each subtree 
calculate_ewIC <- function(tree) {
  return(ewCollessI(tree))
}


#2. contrast matrix for ancestral state
contrast <- matrix(data = c(
  1, 0, 0, 0, 0, 0, 0, 0,  # 0
  0, 1, 0, 0, 0, 0, 0, 0,  # 1
  0, 0, 1, 0, 0, 0, 0, 0,  # 2
  0, 0, 0, 1, 0, 0, 0, 0,  # 3
  0, 0, 0, 0, 1, 0, 0, 0,  # 4
  0, 0, 0, 0, 0, 1, 0, 0,  # 5
  0, 0, 0, 0, 0, 0, 1, 0,  # 6
  0, 0, 0, 0, 0, 0, 0, 1,  # 7
  1, 1, 0, 0, 0, 0, 0, 0,  # 01 (0 or 1)
  1, 0, 1, 0, 0, 0, 0, 0,  # 02 (0 or 2)
  1, 0, 0, 1, 0, 0, 0, 0,  # 03
  1, 0, 0, 0, 1, 0, 0, 0,  # 04
  1, 0, 0, 0, 0, 1, 0, 0,  # 05
  1, 0, 0, 0, 0, 0, 1, 0,  # 06
  1, 0, 0, 0, 0, 0, 0, 1,  # 07
  0, 1, 1, 0, 0, 0, 0, 0,  # 12 (1 or 2)
  0, 1, 0, 1, 0, 0, 0, 0,  # 13
  0, 0, 1, 1, 0, 0, 0, 0,  # 23 (2 or 3)
  0, 1, 0, 0, 0, 1, 0, 0,  # 15 (1 or 5)
  0, 0, 1, 1, 1, 0, 0, 0,  # (234) 
  1, 1, 1, 0, 0, 0, 0, 0,  #. (012)
  0, 0, 0, 1, 1, 0, 0, 0,  # 34 (3 or 4)
  0, 0, 0, 0, 0, 1, 1, 0,   #67
  1, 1, 1, 1, 1, 1, 1, 1,  # ? (Any state including a gap)
  1, 1, 1, 1, 1, 1, 1, 1.  # - 
), ncol = 8, byrow = TRUE)


# Name the rows and columns of the matrix
dimnames(contrast) <- list(c("0", "1", "2", "3", "4", "5", "6", "7", 
                             "(01)", "(02)", "(03)","(04)", "(05)", "(06)", "(07)",  
                             "(12)", "(13)", "(23)", "(15)", "(234)", "(34)",'(012)','(67)', "?", "-"),
                           c("0", "1", "2", "3", "4", "5", "6", "7"))

# Print the contrast matrix


#3.getting synapomorphy from ancestral states, by characetr 
process_ancestral_data <- function(tree, ances, num_characters, char_phyDat) {
  # Combine tip labels and node labels
  all_labels <- c(tree$tip.label, tree$node.label)
  node_id_to_label <- setNames(all_labels, c(1:length(all_labels)))
  
  # Extract the edge list from the tree
  edges <- tree$edge
  
  # Initialize a list to record changes
  changes <- vector("list", num_characters)
  siteidxs <- unlist(attr(char_phyDat, "index")) # Ensure that char_phyDat is correctly passed and used
  
  for (char_idx in 1:num_characters) {
    # Initialize a vector to record changes for the current character
    char_changes <- c()
    
    # Iterate through each edge and compare states
    for (i in 1:nrow(edges)) {
      parent <- edges[i, 1]
      child <- edges[i, 2]
      
      if (!is.null(ances[[parent]]) && !is.null(ances[[child]])) {
        
        parent_state <- ances[[parent]][siteidxs[char_idx],]
        child_state <- ances[[child]][siteidxs[char_idx],]
        
        if (!any(parent_state != 0 & child_state != 0)) {
          # Get the labels for the parent and child
          parent_label <- node_id_to_label[as.character(parent)]
          child_label <- node_id_to_label[as.character(child)]
          
          # Ensure labels are not NA by using node IDs as fallback
          if (is.na(parent_label)) parent_label <- paste(parent)
          if (is.na(child_label)) child_label <- paste(child)
          
          char_changes <- c(char_changes, sprintf(parent_label))
        }
      }
    }
    
    # Add the changes for the current character to the main list
    changes[[char_idx]] <- char_changes
  }
  
  # Return the list of changes
  return(changes)
}

######
#from synapomorphy per chracter to gfit at branch level:
#1. get the subtrees info from list of subtrees, then convert into dataframe
combined_function <- function(phylo_trees) {
  # Function to extract tree information
  extract_tree_info <- function(phylo_tree) {
    list(
      num_tips = Ntip(phylo_tree),
      num_nodes = Nnode(phylo_tree),
      node_labels = if (!is.null(phylo_tree$node.label)) {
        paste(phylo_tree$node.label, collapse = ", ")
      } else {
        NA
      }
    )
  }
  
  # Convert the list of extracted information to a data frame
  tree_info_list <- lapply(phylo_trees, extract_tree_info)
  df <- do.call(rbind, lapply(tree_info_list, as.data.frame))
  df$subtrees <- seq_len(nrow(df))
  
  # Calculate the imbalance for each phylogenetic tree and add it to the data frame
  df$imbalance <- sapply(phylo_trees, calculate_ewIC)
  
  # Convert the node labels to individual rows
  df_nodes_expanded <- df %>%
    separate_rows(node_labels, sep = ", ") %>%
    mutate(node_labels = as.character(node_labels))
  
  return(df_nodes_expanded)
}


#2. Function to extract nodes, synapomorphy and create a dataframe from 'test'
# datalist is text freviously, move up to test! can be combined
# Initialize an empty list to store the extracted nodes
list_to_list <- function(data_list, tree) {
extracted_list <- list()

# Collect all unique nodes
all_nodes <- unique(unlist(data_list))

# Iterate over the data_list to extract nodes
for (i in seq_along(data_list)) {
  if (!is.null(data_list[[i]])) {
    for (node in data_list[[i]]) {
      if (!is.null(extracted_list[[node]])) {
        extracted_list[[node]] <- c(extracted_list[[node]], i)
      } else {
        extracted_list[[node]] <- i
      }
    }
  }
}

# Ensure all nodes are in the extracted list, set as NA if not present
for (node in all_nodes) {
  if (is.null(extracted_list[[node]])) {
    extracted_list[[node]] <- NA
  }
}

return(extracted_list)
}
# 3.Merge subtree data frame with stats dataframe,
#then calculate the sum and average of Gfit and RI for each subtree
summarise_data <- function(expanded, result_df) {
  
  summary_df <- expanded %>%
    full_join(result_df, by = c("node_labels" = "node")) %>%
    mutate(Gfit = as.numeric(Gfit),
           RI = as.numeric(RI)) %>%
    group_by(subtrees) %>%
    summarise(
      imbalance= unique(imbalance),
      subtrees = unique(subtrees),
      sum_Gfit = sum(Gfit, na.rm = TRUE),
      avg_Gfit = ifelse(any(!is.na(Gfit)), sum(Gfit, na.rm = TRUE) / unique(num_nodes), NA),
      sum_RI = sum(RI, na.rm = TRUE),
      avg_RI = ifelse(any(!is.na(RI)), sum(RI, na.rm = TRUE) / unique(num_nodes), NA)
    )
  
  return(summary_df)
}

#####character completeness
# Function to calculate the proportion of unambiguous states
char_proportion <- function(data) {
  # Initialize a vector to store the proportions
  proportions <- numeric(nrow(data))
  
  # Loop through each row of the data
  for (i in 1:nrow(data)) {
    # Count the total number of characters in the row
    total_characters <- ncol(data)
    
    # Count the number of unambiguous characters 
    unambiguous_characters <- sum(data[i, ] != '-' & data[i, ] != '?')

    proportions <- unambiguous_characters / total_characters
  }
  
  # Return the proportions
  return(proportions)
}


