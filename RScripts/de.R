require(methods)

library(limma)
library(EnhancedVolcano)
library("gplots")
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)
library(dplyr)

#n1<-83
#n2<-66
#Load "normalized_data_matrix" and the "design_matrix".
#normalized_data_matrix: each column represents a sample (this is a unique file derived from the fusion of samples inside the two groups), and each row represents a detected gene.
#design matrix: in our case each column represents a group (i.e. group1, group2), and each row corresponds to a sample in the normalized_data_matrix.

#HERE EXAMPLE OF COMPARISON AMONG DKD AND NORMAL GLOMERULI

matrix_data <- read.delim("matrix_group2_group1.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", row.names = 1)
matrix_design <- read.delim("matrix_design_final.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", row.names = 1)

#Create a Contrast matrix: Then, we must tell limma whom we are going to compare with whom.
cont_matrix <- makeContrasts(group1vsgroup2 = group1-group2, levels = matrix_design)   #levels means samples

#Fit the expression matrix to a linear model
fit <- lmFit(matrix_data, matrix_design)

#Compute contrast
fit_contrast <- contrasts.fit(fit, cont_matrix)

#Bayes statistics of differential expression
fit_contrast1 <- eBayes(fit_contrast)


#Summary of results (number of differentially expressed genes). 
result <- decideTests(fit_contrast1)
summary(result)

#Generate a list of top 'n' differentially expressed genes. The user may choose this 'n' number (If n > of DE_genes --> n = all_DE_genes)
DE_table1 <- topTable(fit_contrast1, adjust="fdr", n=500)

#print the table. Here the user could save this table.
DE_table1

write.csv(DE_table1, file = "DE_genes.txt", row.names = TRUE)
#To prepare heatmap, you select rows/genes of "DE_table1" from "matrix_data".
table1 <- matrix_data[row.names(DE_table1),]
map1 <- data.matrix(table1)


#Prepare ColSideColors of heatmap
x <- 0
z <- NULL
#83 = ncol(group1_file)
while(x<n1) { 	
x <- x+1
z <- c(z,"red3") }

x <- 0
#66 = ncol(group2_file)
while(x<n2) { 	
x <- x+1
z <- c(z,"royalblue") }

#heatmap
svg("Heatmap.svg",width=15,height=30)		#to open in high resolution the plot, the problem is that 'width' and 'height' change on the basis of how much is big the plot!
heatmap.2(map1, col=redgreen(75), key=TRUE, symkey=FALSE, density.info="none", cexRow=0.1, trace = "none", scale = "row",  margins=c(12,15), ColSideColors = z)
legend(x="right", legend=c("group2", "group1"), fill=c("red", "royalblue"))
dev.off()

#VolcanoPlot
svg("VolcanoPlot.svg",width=15,height=8)    		#to open in high resolution the plot, the problem is that 'width' and 'height' change on the basis of how much is big the plot!
#EnhancedVolcano(fit_contrast1, lab = rownames(fit_contrast1), x = "coefficients", y = "p.value", title = 'Normal-vs-group2 Glomeruli', pCutoff = 10e-10, FCcutoff = 1, pointSize = 1.0, labSize = 3.0) #also these sizes may changes...
volcanoplot(fit_contrast1, style = "p-value", highlight = 10, names = rownames(fit_contrast1), hl.col="red", xlab = "Log2 Fold Change", ylab = NULL, pch=1, cex=0.5)
#The user may choose 'pCutoff' and 'FCcutoff' 
dev.off()


#For CELL comparison
#IMPORTANT: delete the column with sample!
matrix_cell <- read.delim("cell_new_final.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".")

x <- 0
for (name in colnames(matrix_cell)) {
x <- x+1
#print(x)
if (x>ncol(matrix_cell)-1){break}
res <- oneway.test(matrix_cell[[name]] ~ Group, data = matrix_cell)
if (res$p.value < 0.05 & res$p.value != "NaN") { 	#p.value chosen by user
print(name)
print(res$p.value)
p <- ggboxplot(matrix_cell, x = "Group", y = name, color = "Group", palette = c("#00AFBB", "#E7B800"), add = "jitter") 
svg(paste0(name, ".svg"),width=5,height=5)
print(p + stat_compare_means(method = "t.test", label.x = 1.5))
dev.off()
			} 
				}

#For GO comparison
#IMPORTANT: delete the column with sample!
matrix_GO <- read.delim("GO_final.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".")

x <- 0
for (name in colnames(matrix_GO)) {
x <- x+1
#print(x)
if (x>ncol(matrix_GO)-1){break}
res <- oneway.test(matrix_GO[[name]] ~ Group, data = matrix_GO)
if (res$p.value < 0.00001 & res$p.value != "NaN") { 	#p.value chosen by user	
print(name)
print(res$p.value)
			} 
				}

