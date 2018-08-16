#!/usr/bin/env Rscript

# USAGE: Rscript interactive-heatmap.r [emissions .txt file] [species names .yaml file] [# clusters] [output file name] [optional -g to specify group]

library(cba)
library(yaml)
library(tidyr)
library(dplyr)
library(tibble)
library(stringr)
library(heatmaply)
library(pheatmap)

args = commandArgs(trailingOnly=TRUE)
emissionsFileName = args[1]
yamlFileName = args[2]
numClusters = as.numeric(args[3])
outputFileName = args[4]

includeGroups = FALSE
if (length(args) == 5) {
    if (args[5] == "-g") {
        includeGroups = TRUE
    }
}

# PROCESSING EMISSIONS FILE
emissions = read.table(emissionsFileName, header=TRUE, sep="\t")
species_names_yaml = yaml.load_file(yamlFileName)
species_names = as.data.frame(cbind(names(species_names_yaml), matrix(unlist(species_names_yaml), ncol=length(unlist(species_names_yaml[1])), byrow = TRUE)))
rownames(species_names) = NULL
if (includeGroups == TRUE) {
    colnames(species_names) = cbind("genomeName", "commonName", "group", "distanceToHuman")
} else {
    colnames(species_names) = cbind("genomeName", "commonName", "distanceToHuman")
}
numSpecies = (ncol(emissions) - 1) / 2
emissions_matched_replaced = emissions
columns = c()
for (i in c(2:ncol(emissions))) {
    species = sub("_aligned", "", colnames(emissions)[i])
    species = sub("_matched", "", species)
    emission_type = sub(".*_", "", colnames(emissions)[i])
    columns = rbind(columns, c(species, emission_type))
      
    if (emission_type == "matched") {
        emissions_matched_replaced[i] = emissions[,i] * emissions[[paste(species, "aligned", sep="_")]]
    }
}
columns = as.data.frame(columns)
colnames(columns) = c("genomeName", "emission_type")
columns$index = c(2:(numSpecies*2 + 1))
columns = columns %>% mutate(fullName = paste(columns[,1], columns[,2], sep="_"))

emissions_columns_reordered = emissions_matched_replaced
columns = merge(columns, species_names, by="genomeName")
columns$distanceToHuman = as.numeric(as.character(columns$distanceToHuman))
newcolnames = colnames(emissions_matched_replaced)
for (i in c(1:nrow(columns))) {
      if (columns$emission_type[i] == "aligned") {
        emissions_columns_reordered[,(columns$distanceToHuman[i] + 1)] = emissions_matched_replaced[,columns$index[i]]
        newcolnames[columns$distanceToHuman[i] + 1] = as.character(columns$commonName[i])
      } else {
            emissions_columns_reordered[,(columns$distanceToHuman[i] + (numSpecies + 1))] = emissions_matched_replaced[,columns$index[i]]
            newcolnames[columns$distanceToHuman[i] + (numSpecies + 1)] = as.character(columns$commonName[i])
          }
}

OLO_clustering = function(hc, mat) {
    d <- dist(mat)
    hc <- hclust(d)
    co <- order.optimal(d, hc$merge)
    ho <- hc
    ho$merge <- co$merge
    ho$order <- co$order
    hc <- ho
}

columns = columns %>% arrange(emission_type, distanceToHuman)
colnames(emissions_columns_reordered) = c("state", columns$fullName)
if (includeGroups == TRUE) {
    annotation_col = columns %>% select(group)
    rownames(annotation_col) = columns$fullName
}

draw_colnames_45 <- function (coln, ...) {
    m = length(coln)
    x = (1:m)/m - 1/2/m
    grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
        hjust = 1, rot = 45, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}

assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

pdf(paste0(outputFileName, ".pdf"), width = 33, height = 16, onefile = FALSE)
if (includeGroups == TRUE) {
    pheatmap(emissions_columns_reordered[,-1], cluster_rows = TRUE, cluster_cols = FALSE, clustering_callback = OLO_clustering, cellwidth = 10, 
             cellheight = 9, border_color = NA, labels_col = newcolnames[-1], gaps_col = numSpecies, fontsize_col = 7, annotation_col = annotation_col, 
             cutree_rows = numClusters, display_numbers = TRUE, fontsize_number = 3)
} else {
    pheatmap(emissions_columns_reordered[,-1], cluster_rows = TRUE, cluster_cols = FALSE, clustering_callback = OLO_clustering, cellwidth = 10, 
             cellheight = 9, border_color = NA, labels_col = newcolnames[-1], gaps_col = numSpecies, fontsize_col = 7, cutree_rows = numClusters, 
             display_numbers = TRUE, fontsize_number = 3)
}

dev.off()

emissions_columns_reordered = add_column(emissions_columns_reordered, gap = 0.5, .after = (numSpecies+1))

newcolnames[numSpecies+2] = " "
for (i in (numSpecies+3):(2*numSpecies + 2)) {
    newcolnames[i] = paste0(" ", newcolnames[i-(numSpecies+1)])
}

if (includeGroups == TRUE) {
    columns = subset(columns, select = c("fullName", "group"))
    addgap = data.frame(fullName = "gap", group = " ")
    columns = rbind(columns[1:numSpecies,], addgap, columns[(numSpecies+1):(2*numSpecies),])
    columns = na.omit(columns)
    annotation_col = columns %>% select(group)
    rownames(annotation_col) = newcolnames[-1]
    heatmaply(emissions_columns_reordered[,-1], xlab="Species Names\nAligned(left), Matched(right)", ylab="Conservation States", 
              main="Heatmap of Conservation State by Species",  key.title="Emission\nProbability", Colv = FALSE, labCol = newcolnames[-1], 
              file = paste0(outputFileName, ".html"), k_row=numClusters, scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", 
              high = "red", midpoint=0.5), col_side_colors = annotation_col)
} else {
    heatmaply(emissions_columns_reordered[,-1], xlab="Species Names\nAligned(left), Matched(right)", ylab="Conservation States", 
              main="Heatmap of Conservation State by Species",  key.title="Emission\nProbability", Colv = FALSE, labCol = newcolnames[-1], 
              file = paste0(outputFileName, ".html"), k_row=numClusters, scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", 
              high = "red", midpoint=0.5),)
}



