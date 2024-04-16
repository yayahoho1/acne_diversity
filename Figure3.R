# Figure 3A --------------------------------------------------------------------
# Load required packages and data ----------------------------------------------

library(pheatmap)
#install.packages("magick")
library(magick) #
library(circlize)
library(ComplexHeatmap)

data = read.csv("pathway_complete_with_annotation.csv",header = T,row.names = 1,check.names = F)
groupinfo = read.table("group_info.txt",header = T,row.names = 1,sep = "\t")
pathway_info = read.csv("module_classification_info.csv",header = T,row.names = 1)

# col: Determine whether the sample order is consistent
all(rownames(groupinfo) %in% colnames(data))
all(rownames(groupinfo) == colnames(data))
# Reorder the sample；
data = data[, rownames(groupinfo)]
all(rownames(groupinfo) == colnames(data))

# row: Determine whether the sample order is consistent
all(rownames(pathway_info) %in% rownames(data))
all(rownames(pathway_info) == rownames(data))
# Reorder the pathway；
data = data[rownames(pathway_info),]
all(rownames(pathway_info) == rownames(data))

#
groupinfo$Group = factor(groupinfo$Group,levels = c("NC", "AD", "Acne","ConwillA_2022","ZhouW_2020","JoglekarP_2023","KashafSS_2023"))
pathway_info$Pathway = factor(pathway_info$Pathway,levels=c("Amino acid metabolism", "Carbohydrate metabolism", "Energy metabolism",
                                                            "Lipid metabolism", "Metabolism of cofactors and vitamins",	"Nucleotide metabolism",
                                                            "Biosynthesis of other secondary metabolites", "Biosynthesis of terpenoids and polyketides",
                                                            "Gene set"
))

# color palette
col_fun = colorRamp2(c(0,0.5,1),c("#e7e7e5","#bbded7","#D24735"))

# annotation
hac = HeatmapAnnotation(df = groupinfo,
                        col = list(Group=c("NC" = "#90A6BB", "AD" = "#FFD579","Acne"="#9D4E3F","ConwillA_2022"="#B2DF8AFF","ZhouW_2020"="#FB9A99FF","JoglekarP_2023"="#FFFF99FF","KashafSS_2023"="#60A5DFFF")),
                        which = "col"
)
har = HeatmapAnnotation(df = pathway_info,
                        col = list(Pathway=c("Amino acid metabolism"="#FAEBD7", "Carbohydrate metabolism"="#CDC0B0", 
                                             "Energy metabolism"="#FFA07A",  
                                             "Lipid metabolism"="#CD8C95", "Metabolism of cofactors and vitamins"="#8B5F65",	
                                             "Nucleotide metabolism"="#6E8B3D",
                                             "Biosynthesis of other secondary metabolites"="#C1CDCD", "Biosynthesis of terpenoids and polyketides"="#CDC673",
                                             "Gene set"="#66CDAA")),
                        which = "row"
)

# plot
df=as.matrix(data)
pdf("function_completeness.pdf",width = 15,height = 25)
ComplexHeatmap::Heatmap(df, 
                        col = col_fun,
                        row_names_side = "right",
                        show_row_names = T,
                        column_names_gp = gpar(fontsize = 3),
                        row_names_gp = gpar(fontsize = 13),
                        show_column_names = F,
                        
                        cluster_rows = T,cluster_columns = T,
                        column_split = as.factor(groupinfo$Group),
                        row_split = as.factor(pathway_info$Pathway),
                        cluster_column_slices = F,
                        cluster_row_slices =F,
                        row_names_max_width = unit(2,"cm"),
                        rect_gp = gpar(col = "grey", lwd = 0.000001),
                        top_annotation = hac,
                        left_annotation=har,
                        width = ncol(data)*unit(0.04, "mm"), #
                        height = nrow(t(df))*unit(0.12, "mm"),#
                        heatmap_legend_param = list(title= "Completeness",
                                                    title_position = "topcenter", 
                                                    legend_direction= "horizontal",
                                                    legend_width=unit(5,"cm")
                        )
)



dev.off()


# Figure 3B,C ------------------------------------------------------------------
# Load required packages and data ----------------------------------------------

library(pheatmap)
library(magick) #
library(circlize)
library(ComplexHeatmap)

data=read.table("VF_summary.txt",header = T,row.names = 1,check.names = F)
groupinfo=read.table("group_info.txt",header = T,row.names = 1,sep = "\t")

# Determine whether the sample order is consistent
all(rownames(groupinfo) %in% colnames(data))
all(rownames(groupinfo) == colnames(data))

# Reorder the sample
data = data[, rownames(groupinfo)]
all(rownames(groupinfo) == colnames(data))

#
groupinfo$Group = factor(groupinfo$Group,levels = c("NC", "AD", "Acne","ConwillA_2022","ZhouW_2020","JoglekarP_2023","KashafSS_2023"))

# color palette
col_fun = colorRamp2(c(0,1),c("#DCDCDC","#B73508"))

#annotation
hac = HeatmapAnnotation(df = groupinfo,
                        col = list(Group=c("NC" = "#90A6BB", "AD" = "#FFD579","Acne"="#9D4E3F","ConwillA_2022"="#B2DF8AFF","ZhouW_2020"="#FB9A99FF","JoglekarP_2023"="#FFFF99FF","KashafSS_2023"="#60A5DFFF")),
                        which = "col"
)

# plot
df = as.matrix(data)
ComplexHeatmap::Heatmap(df, 
                          col = col_fun,
                          row_names_side = "right",
                          show_row_names = T,
                          column_names_gp = gpar(fontsize = 10),
                          row_names_gp = gpar(fontsize = 9),
                          show_column_names = F,
                          cluster_rows = T,cluster_columns = T,
                          column_split = as.factor(groupinfo$Group),
                          cluster_column_slices = F,
                          row_names_max_width = unit(10,"cm"),
                          rect_gp = gpar(col = "grey", lwd = 0.000001),
                          top_annotation = hac,
                          width = ncol(data)*unit(0.05, "mm"), #
                          height = nrow(t(df))*unit(0.05, "mm"),#
                          show_heatmap_legend = F,
                          
)

