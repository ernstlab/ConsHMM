Post-processing of ConsHMM state parameters
================

The ConsHMM pipeline learns a conservation state model from a multiple-sequence alignment of different genomes. This notebook was used to visualize the parameters of a 100-state ConsHMM model learned from the hg19 Multiz 100-way alignment. It can be used to visualize any ConsHMM model, by editing the input files.

### Required libraries.

``` r
library(cba)
library(yaml)
library(pheatmap)
library(tidyverse)
library(gtable)
library(grid)
library(gridExtra)
```

Input parameters
----------------

This notebook requires an emissions parameter file created by running ChromHMM on a multiple-sequence alignment processed through ConsHMM. Specify the name of this file using the `emissionsFileName` variable below:

``` r
emissionsFileName = "../models/hg19_multiz100way/emissions_100.txt"
genomeName <- unlist(strsplit(emissionsFileName, split='/', fixed=TRUE))[3]
file_name <- unlist(strsplit(emissionsFileName, split='/', fixed=TRUE))[5]
file_name <- substr(file_name, 11, nchar(file_name))
numStates <- as.numeric(unlist(strsplit(file_name, split='.', fixed=TRUE))[1])
```

To help interpret the emission parameters, provide a `.yaml` file that annotates the genome names used in the multiple sequence alignment with a common name (used for readability), a distance to human (used for arranging columns in the order of a phylogenetic tree), and a group (used to group columns into clades). Each species should have an entry in the `.yaml` file such as:

``` yaml
panTro4: 
 commonName: Chimp
 group: Primate
 distanceToHuman: 1
```

Specify the name of the `.yaml` file using the `yamlFileName` variable below. If the yaml file contains group labels for each species then set `groups` to FALSE, otherwise set to TRUE.

``` r
yamlFileName = paste0("../speciesNames/", genomeName, "/speciesNames.yaml")
groups = TRUE
```

The states of the model are clustered hierarchically using optimal leaf ordering. To cut the dendrogram into a set number of clusters specify the number of clusters using the `numClusters` variable below. The resulting heatmap will insert white space between the clusters. Set to 0 for no set number of clusters. The full dendrogram is always displayed.

``` r
numClusters = 6
```

The notebook creates a large heatmap that is best viewed by saving into `.pdf` format. Specify the name of the output file below.

``` r
outputFileName = paste0("../emissionHeatmaps/", genomeName, "_", numStates, "_emissions.pdf")
```

Processing emissions file
-------------------------

Reading in and renaming data.

``` r
emissions = read.table(emissionsFileName, header=TRUE, sep="\t")
species_names_yaml = yaml.load_file(yamlFileName)
species_names = as.data.frame(cbind(names(species_names_yaml), matrix(unlist(species_names_yaml), ncol=length(unlist(species_names_yaml[1])), byrow = TRUE)))
rownames(species_names) = NULL
if (groups) {
    colnames(species_names) = cbind("genomeName", "commonName", "group", "distanceToHuman")
} else {
    colnames(species_names) = cbind("genomeName", "commonName", "distanceToHuman")
}

numSpecies = (ncol(emissions) - 1) / 2 
numStates = nrow(emissions)
```

Due to internal ChromHMM representation, the matching probabilities will not be meaningful when the corresponding aligning probability is low. To correct for this we replace the matching probabilites are replaced with the corresponding align \* match probability.

``` r
emissions_matched_replaced = emissions
columns = c()
for (i in c(2:ncol(emissions))) {
  split_column_name = unlist(strsplit(colnames(emissions)[i], split="_", fixed=TRUE))
  species = paste(split_column_name[1:(length(split_column_name) - 1)], collapse='_')
  emission_type = split_column_name[length(split_column_name)]
  columns = rbind(columns, c(species, emission_type))
  
  if (emission_type == "matched") {
    emissions_matched_replaced[i] = emissions[,i] * emissions[[paste(species, "aligned", sep="_")]]
  }
}
columns = as.data.frame(columns)
colnames(columns) = c("genomeName", "emission_type")
columns$index = c(2:((numSpecies*2) + 1))
columns = columns %>% mutate(fullName = paste(columns[,1], columns[,2], sep="_"))
```

Reordering columns according to distance to human.

``` r
emissions_columns_reordered = emissions_matched_replaced
columns = merge(columns, species_names, by="genomeName")
columns$distanceToHuman = as.numeric(as.character(columns$distanceToHuman))
newcolnames = colnames(emissions_matched_replaced)
for (i in c(1:nrow(columns))) {
  if (columns$emission_type[i] == "aligned") {
    emissions_columns_reordered[,(columns$distanceToHuman[i] + 1)] = emissions_matched_replaced[,columns$index[i]] 
    newcolnames[columns$distanceToHuman[i] + 1] = as.character(columns$commonName[i])
  } else {
    emissions_columns_reordered[,(columns$distanceToHuman[i] + numSpecies + 1)] = emissions_matched_replaced[,columns$index[i]]
    newcolnames[columns$distanceToHuman[i] + numSpecies + 1] = as.character(columns$commonName[i])
  }
}
```

Clustering function used for clutering rows hierarchially according to optimal leaf ordering.

``` r
OLO_clustering = function(hc, mat) {
  d <- dist(mat)
  hc <- hclust(d)
  co <- order.optimal(d, hc$merge)
  ho <- hc
  ho$merge <- co$merge
  ho$order <- co$order
  hc <- ho
}
```

Heatmap plotting
----------------

Pheatmap modifications

``` r
convert_annotations = function(annotation, annotation_colors){
    new = annotation
    for(i in 1:ncol(annotation)){
        a = annotation[, i]
        b = annotation_colors[[colnames(annotation)[i]]]
        if(is.character(a) | is.factor(a)){
            a = as.character(a)
            
            if(length(setdiff(setdiff(a, NA), names(b))) > 0){
                stop(sprintf("Factor levels on variable %s do not match with annotation_colors", colnames(annotation)[i]))
            }
            new[, i] = b[a]
        }
        else{
            a = cut(a, breaks = 100)
            new[, i] = colorRampPalette(b)(100)[a]
        }
    }
    return(as.matrix(new))
}

lo <- function (rown, coln, nrow, ncol, cellheight = NA, cellwidth = NA, 
                treeheight_col, treeheight_row, legend, annotation_row, annotation_col, 
                annotation_colors, annotation_legend, annotation_names_row, 
                annotation_names_col, main1, main2, fontsize, fontsize_row, fontsize_col, 
                angle_col, gaps_row, gaps_col, ...) 
{
    if (!is.null(coln[1]) | (!pheatmap:::is.na2(annotation_row) & annotation_names_row)) {
        if (!is.null(coln[1])) {
            t = coln
        }
        else {
            t = ""
        }
        tw = strwidth(t, units = "in", cex = fontsize_col/fontsize)
        if (annotation_names_row) {
            t = c(t, colnames(annotation_row))
            tw = c(tw, strwidth(colnames(annotation_row), units = "in"))
        }
        longest_coln = which.max(tw)
        gp = list(fontsize = ifelse(longest_coln <= length(coln), 
                                    fontsize_col, fontsize), ...)
        coln_height = unit(1, "grobheight", textGrob(t[longest_coln], 
                                                     rot = angle_col, gp = do.call(gpar, gp))) + unit(10, 
                                                                                                      "bigpts")
    }
    else {
        coln_height = unit(5, "bigpts")
    }
    if (!is.null(rown[1])) {
        t = rown
        tw = strwidth(t, units = "in", cex = fontsize_row/fontsize)
        if (annotation_names_col) {
            t = c(t, colnames(annotation_col))
            tw = c(tw, strwidth(colnames(annotation_col), units = "in"))
        }
        longest_rown = which.max(tw)
        gp = list(fontsize = ifelse(longest_rown <= length(rown), 
                                    fontsize_row, fontsize), ...)
        rown_width = unit(1, "grobwidth", textGrob(t[longest_rown], 
                                                   rot = 0, gp = do.call(gpar, gp)))
    }
    else {
        rown_width = unit(5, "bigpts")
    }
    gp = list(fontsize = fontsize, ...)
    if (!pheatmap:::is.na2(legend)) {
        longest_break = which.max(nchar(names(legend)))
        longest_break = unit(1.1, "grobwidth", 
                             textGrob(as.character(names(legend))[longest_break], 
                                      gp = do.call(gpar, gp)))
        title_length = unit(1.1, "grobwidth", textGrob("Scale", 
                                                       gp = gpar(fontface = "bold", ...)))
        legend_width = unit(12, "bigpts") + longest_break * 1.2
        legend_width = max(title_length, legend_width) + unit(10, "bigpts")
    }
    else {
        legend_width = unit(0, "bigpts")
    }
    if (is.na(main1)) {
        main1_height = unit(0, "npc")
    }
    else {
        main1_height = unit(1.5, "grobheight", textGrob(main1, 
                                                        gp = gpar(fontsize = 1.3 * fontsize, ...)))
    }
    textheight = unit(fontsize, "bigpts")
    if (!pheatmap:::is.na2(annotation_col)) {
        annot_col_height = ncol(annotation_col) * (textheight + 
                                                       unit(2, "bigpts")) + unit(2, "bigpts")
        t = c(as.vector(as.matrix(annotation_col)), colnames(annotation_col))
        annot_col_legend_width = unit(1.2, "grobwidth", textGrob(t[which.max(nchar(t))], 
                                                                 gp = gpar(...))) + unit(12, "bigpts")
        if (!annotation_legend) {
            annot_col_legend_width = unit(0, "npc")
        }
    }
    else {
        annot_col_height = unit(0, "bigpts")
        annot_col_legend_width = unit(0, "bigpts")
    }
    if (!pheatmap:::is.na2(annotation_row)) {
        annot_row_width = ncol(annotation_row) * (textheight + 
                                                      unit(2, "bigpts")) + unit(2, "bigpts")
        t = c(as.vector(as.matrix(annotation_row)), colnames(annotation_row))
        annot_row_legend_width = unit(1.2, "grobwidth", textGrob(t[which.max(nchar(t))], 
                                                                 gp = gpar(...))) + unit(12, "bigpts")
        if (!annotation_legend) {
            annot_row_legend_width = unit(0, "npc")
        }
    }
    else {
        annot_row_width = unit(0, "bigpts")
        annot_row_legend_width = unit(0, "bigpts")
    }
    annot_legend_width = max(annot_row_legend_width, annot_col_legend_width)
    treeheight_col = unit(treeheight_col, "bigpts") + unit(5, 
                                                           "bigpts")
    treeheight_row = unit(treeheight_row, "bigpts") + unit(5, 
                                                           "bigpts")
    if (is.na(cellwidth)) {
        mat_width = unit(1, "npc") - rown_width - legend_width - 
            treeheight_row - annot_row_width - annot_legend_width
    }
    else {
        mat_width = unit(cellwidth * ncol, "bigpts") + length(gaps_col) * 
            unit(4, "bigpts")
    }
    if (is.na(cellheight)) {
        mat_height = unit(1, "npc") - main1_height - coln_height - 
            treeheight_col - annot_col_height
    }
    else {
        mat_height = unit(cellheight * nrow, "bigpts") + length(gaps_row) * 
            unit(4, "bigpts")
    }
    gt = gtable(widths = unit.c(treeheight_row, rown_width, mat_width * 0.5, mat_width * 0.5, treeheight_row, legend_width, annot_legend_width), heights = unit.c(main1_height, treeheight_col, annot_col_height, mat_height, coln_height, treeheight_col), vp = viewport(gp = do.call(gpar, gp)))
    #gt = gtable(widths = unit.c(treeheight_row, rown_width,
    #    mat_width, treeheight_row, legend_width, annot_legend_width), 
    #    heights = unit.c(main_height, treeheight_col, annot_col_height, 
    #        mat_height, coln_height), vp = viewport(gp = do.call(gpar, 
    #        gp)))
    cw = convertWidth(mat_width - (length(gaps_col) * unit(4, 
                                                           "bigpts")), "bigpts", valueOnly = T)/ncol
    ch = convertHeight(mat_height - (length(gaps_row) * unit(4, 
                                                             "bigpts")), "bigpts", valueOnly = T)/nrow
    mindim = min(cw, ch)
    res = list(gt = gt, mindim = mindim)
    return(res)
}

# Modified pheatmap:::draw_rownames      
draw_rownames <- function (rown, gaps, ...) 
{
    coord = pheatmap:::find_coordinates(length(rown), gaps)
    y = unit(1, "npc") - (coord$coord - 0.5 * coord$size)
    res = textGrob(rown, x = unit(-3, "bigpts"), y = y, vjust = 0.5, 
                   hjust = 1, gp = gpar(...))
    return(res)
}

heatmap_motor <- function (matrix, border_color, cellwidth, cellheight, tree_col, 
                           tree_row, treeheight_col, treeheight_row, filename, width, 
                           height, breaks, color, legend, annotation_row, annotation_col, 
                           annotation_colors, annotation_legend, annotation_names_row, 
                           annotation_names_col, main1, main2, fontsize, fontsize_row, fontsize_col, 
                           hjust_col, vjust_col, angle_col, fmat, fontsize_number, number_color, 
                           gaps_col, gaps_row, labels_row, labels_col, ...) 
{
    lo = pheatmap:::lo(coln = labels_col, rown = labels_row, nrow = nrow(matrix), 
                       ncol = ncol(matrix), cellwidth = cellwidth, cellheight = cellheight, 
                       treeheight_col = treeheight_col, treeheight_row = treeheight_row, 
                       legend = legend, annotation_col = annotation_col, annotation_row = annotation_row, 
                       annotation_colors = annotation_colors, annotation_legend = annotation_legend, 
                       annotation_names_row = annotation_names_row, annotation_names_col = annotation_names_col, 
                       main1 = main1, main2 = main2, fontsize = fontsize, fontsize_row = fontsize_row, 
                       fontsize_col = fontsize_col, angle_col = angle_col, gaps_row = gaps_row, 
                       gaps_col = gaps_col, ...)
    res = lo$gt
    mindim = lo$mindim
    if (!is.na(filename)) {
        if (is.na(height)) {
            height = convertHeight(gtable_height(res), "inches", valueOnly = T)
        }
        if (is.na(width)) {
            width = convertWidth(gtable_width(res), "inches", valueOnly = T)
        }
        r = regexpr("\\.[a-zA-Z]*$", filename)
        if (r == -1) 
            stop("Improper filename")
        ending = substr(filename, r + 1, r + attr(r, "match.length"))
        f = switch(ending, pdf = function(x, ...) pdf(x, ...), 
                   png = function(x, ...) png(x, units = "in", res = 300, 
                                              ...), jpeg = function(x, ...) jpeg(x, units = "in", 
                                                                                 res = 300, ...), jpg = function(x, ...) jpeg(x, 
                                                                                                                              units = "in", res = 300, ...), tiff = function(x, 
                                                                                                                                                                             ...) tiff(x, units = "in", res = 300, compression = "lzw", 
                                                                                                                                                                                       ...), bmp = function(x, ...) bmp(x, units = "in", 
                                                                                                                                                                                                                        res = 300, ...), stop("File type should be: pdf, png, bmp, jpg, tiff"))
        f(filename, height = height, width = width)
        gt = heatmap_motor(matrix, cellwidth = cellwidth, cellheight = cellheight, 
                           border_color = border_color, tree_col = tree_col, 
                           tree_row = tree_row, treeheight_col = treeheight_col, 
                           treeheight_row = treeheight_row, breaks = breaks, 
                           color = color, legend = legend, annotation_col = annotation_col, 
                           annotation_row = annotation_row, annotation_colors = annotation_colors, 
                           annotation_legend = annotation_legend, annotation_names_row = annotation_names_row, 
                           annotation_names_col = annotation_names_col, filename = NA, 
                           main1 = main1, main2 = main2, fontsize = fontsize, fontsize_row = fontsize_row, 
                           fontsize_col = fontsize_col, hjust_col = hjust_col, 
                           vjust_col = vjust_col, angle_col = angle_col, fmat = fmat, 
                           fontsize_number = fontsize_number, number_color = number_color, 
                           labels_row = labels_row, labels_col = labels_col, 
                           gaps_col = gaps_col, gaps_row = gaps_row, ...)
        grid.draw(gt)
        dev.off()
        return(gt)
    }
    if (mindim < 3) 
        border_color = NA
    if (!is.na(main1)) {
        elem = pheatmap:::draw_main(main1, fontsize = 2 * fontsize, hjust=1)
        res = gtable_add_grob(res, elem, t = 1, l = 3, name = "main1", 
                              clip = "off")
    }
    if (!is.na(main2)) {
        elem = pheatmap:::draw_main(main2, fontsize = 2 * fontsize, hjust=1)
        res = gtable_add_grob(res, elem, t = 1, l = 4, name = "main2", 
                              clip = "off")
    }
    if (!pheatmap:::is.na2(tree_col) & treeheight_col != 0) {
        elem = pheatmap:::draw_dendrogram(tree_col, gaps_col, horizontal = T)
        res = gtable_add_grob(res, elem, t = 2, l = 3, r = 4, name = "col_tree")
    }
    if (!pheatmap:::is.na2(tree_row) & treeheight_row != 0) {
        elem = pheatmap:::draw_dendrogram(tree_row, gaps_row, horizontal = F)
        res = gtable_add_grob(res, elem, t = 4, l = 1, name = "row_tree")
    }
    elem = pheatmap:::draw_matrix(matrix, border_color, gaps_row, gaps_col, 
                                  fmat, fontsize_number, number_color)
    res = gtable_add_grob(res, elem, t = 4, l = 3, r = 4, clip = "off", 
                          name = "matrix")
    if (length(labels_col) != 0) {
        pars = list(labels_col, gaps = gaps_col, fontsize = fontsize_col, 
                    hjust_col = hjust_col, vjust_col = vjust_col, angle_col = angle_col, 
                    ...)
        elem = do.call(pheatmap:::draw_colnames, pars)
        res = gtable_add_grob(res, elem, t = 5, l = 3, r=4, clip = "off", 
                              name = "col_names")
    }
    if (length(labels_row) != 0) {
        pars = list(labels_row, gaps = gaps_row, fontsize = fontsize_row, 
                    ...)
        elem = do.call(pheatmap:::draw_rownames, pars)
        res = gtable_add_grob(res, elem, t = 4, l = 3, clip = "off", 
                              name = "row_names")
    }
    if (!pheatmap:::is.na2(annotation_col)) {
        converted_annotation = convert_annotations(annotation_col, 
                                                   annotation_colors)
        elem = pheatmap:::draw_annotations(converted_annotation, border_color, 
                                           gaps_col, fontsize, horizontal = T)
        res = gtable_add_grob(res, elem, t = 3, l = 3, r=4, clip = "off", 
                              name = "col_annotation")
        if (annotation_names_col) {
            elem = pheatmap:::draw_annotation_names(annotation_col, fontsize, 
                                                    horizontal = T)
            res = gtable_add_grob(res, elem, t = 3, l = 5, clip = "off", 
                                  name = "col_annotation_names")
        }
    }
    if (!pheatmap:::is.na2(annotation_row)) {
        converted_annotation = convert_annotations(annotation_row, 
                                                   annotation_colors)
        elem = pheatmap:::draw_annotations(converted_annotation, border_color, 
                                           gaps_row, fontsize, horizontal = F)
        res = gtable_add_grob(res, elem, t = 4, l = 2, clip = "off", 
                              name = "row_annotation")
        if (annotation_names_row) {
            elem = pheatmap:::draw_annotation_names(annotation_row, fontsize, 
                                                    horizontal = F, hjust_col = hjust_col, vjust_col = vjust_col, 
                                                    angle_col = angle_col)
            res = gtable_add_grob(res, elem, t = 5, l = 2, clip = "off", 
                                  name = "row_annotation_names")
        }
    }
    annotation = c(annotation_col[length(annotation_col):1], 
                   annotation_row[length(annotation_row):1])
    annotation = annotation[unlist(lapply(annotation, function(x) !pheatmap:::is.na2(x)))]
    if (length(annotation) > 0 & annotation_legend) {
        elem = pheatmap:::draw_annotation_legend(annotation, annotation_colors, 
                                                 border_color, fontsize = fontsize, ...)
        t = ifelse(is.null(labels_row), 4, 3)
        res = gtable_add_grob(res, elem, t = t, l = 7, b = 5,  
                              clip = "off", name = "annotation_legend")
    }
    if (!pheatmap:::is.na2(legend)) {
        elem = pheatmap:::draw_legend(color, breaks, legend, fontsize = fontsize, 
                                      ...)
        t = ifelse(is.null(labels_col), 4, 3)
        res = gtable_add_grob(res, elem, t = t, l = 6, b = 5, 
                              clip = "off", name = "legend")
    }
    return(res)
}

draw_colnames <- function (coln, ...) {
    m = length(coln)
    x = (1:m)/m - 1/2/m
    grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
              hjust = 1, rot = 45, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}

assignInNamespace(x="draw_rownames", value=draw_rownames, ns="pheatmap")
assignInNamespace(x="lo", value=lo, ns="pheatmap")
assignInNamespace(x="heatmap_motor", value=heatmap_motor, ns="pheatmap") 
assignInNamespace(x="draw_colnames", value=draw_colnames, ns="pheatmap")
```

The plotting below was optimized for the 100 state ConsHMM model learned from the hg19 Multiz 100-way alignment. Use the code below to specify plotting parameters to the pheatmap package that best fit the model you are analyzing.

``` r
columns = columns %>% arrange(emission_type, distanceToHuman)
colnames(emissions_columns_reordered) = c("state", columns$fullName)
if (groups) {
    annotation_col = columns %>% select(group)
    rownames(annotation_col) = columns$fullName
} else {
    annotation_col = NA
}

olo_order <- OLO_clustering(oc, emissions_columns_reordered[,-1])$order
state_mapping <- tibble("olo_order" = olo_order, "new_order" = c(1:numStates))

plot_list = list()

pdf(paste0("../emissionHeatmaps/", genomeName, "_emissions_", numStates, "_states.pdf"), height=16, width=34)
pheatmap(emissions_columns_reordered[,-1], main1 = "Align probablities", main2="Match probabilities", cluster_rows = TRUE, cluster_cols = FALSE, clustering_callback = OLO_clustering, cellwidth = 10, cellheight = 9, border_color = NA, labels_col = newcolnames[-1], labels_row = (state_mapping %>% arrange(olo_order))$new_order, gaps_col = numSpecies, fontsize_col = 7, annotation_col = annotation_col, display_numbers = TRUE, fontsize_number = 3, legend_breaks = c(min(emissions_columns_reordered), 0.5, max(emissions_columns_reordered[,-1])), legend_labels = c(0,0.5,1))
dev.off()
```

    ## pdf 
    ##   3

``` r
#plot_list[[1]] = heatmap[[4]]
```

If you would like your downstream analyses to use a segmentation in which the states are named according to your reordered emissions, the code below will output two files which can be used to rename all your model and segmentation files to the OLO ordering.

``` r
# Output state mapping

write_tsv(state_mapping %>% arrange(olo_order), paste0("../stateRenamingFiles/", genomeName, "_state_renaming_", numStates, ".tsv"), col_names = FALSE)
write(colnames(emissions_columns_reordered)[2:length(colnames(emissions_columns_reordered))], paste0("../stateRenamingFiles/", genomeName, "_column_reordering_", numStates, ".tsv"))
```

To reorder your model and segmentation files, use the `Reorder` command in ChromHMM as below: `java -jar ChromHMM.jar Reorder -noimage -f <column reordering file> -o <state renaming file> \ -r <original segmentation file> <reordered_segmentation_file> <original model file> <reordered model file>`
