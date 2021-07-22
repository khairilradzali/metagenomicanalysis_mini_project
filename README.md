# metagenomicanalysis_mini_project

# Metagenomic Analysis - Mini Project using Phyloseq R Package
Credit to: Torondel, B., Ensink, J. H., Gundogdu, O., Ijaz, U. Z., Parkhill, J., Abdelahi, F., ... & Quince, C. (2016). Assessment of the influence of intrinsic environmental and geographical factors on the bacterial ecology of pit latrines. Microbial biotechnology, 9(2), 209-223.

## Brief Introduction
The study demonstrated the bacterial ecology and composition of pit latrines together with inherent latrine characteristic at an unprecedented level of details. Physicochemical data are included to observe any predictor variables that drove the bacterial composition.

## The analysis including:
1. Rarefaction Curve
2. Taxonomic Composition
3. Alpha and Beta Diversity
4. Differential Abundance
5. Network Analysis (igraph R package and Cytoscape)

## 1. Package installation
```
install.packages("glue")
install.packages("ggplot2")
install.packages("vegan")
install.packages("dplyr")
install.packages("scales")
install.packages("reshape2")
install.packages("tidyr")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
```

## 2. Load libraries
```
library(glue)
library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(tidyr)
library(DESeq2)
```

## 3. Set plotting theme
```
theme_set(theme_bw())
```

## 4. Read the input files
```
#Read and transpose OTU table
abund_table<-read.csv("All_Good_P2_C03.csv",row.names=1,check.names=FALSE)
abund_table<-t(abund_table)

#Read, split, combine, and define metadata
meta_table<-read.csv("ENV_pitlatrine.csv",row.names=1,check.names=FALSE)
grouping_info<-data.frame(row.names=rownames(meta_table),t(as.data.frame(strsplit(rownames(meta_table),"_"))))                      
colnames(grouping_info)<-c("Countries","Latrine","Depth")
meta_table1 <- cbind(meta_table,grouping_info)

meta_table2 <- meta_table1 %>%
  mutate(Country = case_when(
    endsWith(Countries, "T") ~ "Tanzania",
    endsWith(Countries, "V") ~ "Vietnam"
  ))


#Read taxonomy table
OTU_taxonomy<-read.csv("All_Good_P2_C03_Taxonomy.csv",row.names=1,check.names=FALSE)
```

## 5. Formating the input files into Phyloseq
```
OTU = otu_table(as.matrix(abund_table), taxa_are_rows = FALSE)
TAX = tax_table(as.matrix(OTU_taxonomy))
SAM = sample_data(meta_table2)
```

## 6. Creating phyloseq object
```
physeq<-merge_phyloseq(phyloseq(OTU,TAX,SAM))
physeq

# Make a data frame with a column for the read counts of each sample
sample_sum_df <- data.frame(sum = sample_sums(physeq))

# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

# mean, max and min of sample read counts
smin <- min(sample_sums(physeq))
smean <- mean(sample_sums(physeq))
smax <- max(sample_sums(physeq))

smin
smean
smax
```

# Objective 1: Rarefaction curve
```
rarecurve((otu_table(physeq)), step=50, cex=0.5)
```

# Objective 2: Taxonomic composition (Phylum and Family rank)
```
# melt to long format (for ggploting) 
# prune out phyla below 2% in each sample
physeq_phylum <- physeq %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

# Set colors for plotting
phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
)


# Plot
p<-ggplot(physeq_phylum,aes(Sample,Abundance,fill=Phylum))+
  geom_bar(stat="identity", position = "fill")+
  facet_grid(.~Country, drop=TRUE,scale="free",space="free_x")
p<-p+scale_fill_manual(values = phylum_colors)  
p<-p+theme_bw()+ylab("Proportions")  
p<-p+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.margin = unit(0.3, "lines"))
p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
print(p)

####Family Level
# prune out class below 2% in each sample
physeq_family <- physeq %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Family)                                      # Sort data frame alphabetically by phylum

# Set colors for plotting
family_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "darkorchid4", "darkgreen", "red3", 
  "yellowgreen", "slategray4", "palevioletred1", "lightsalmon2", "indianred4",
  "green3", "khaki2", "darkseagreen2", "cyan2", "gray0",
  "aquamarine4", "chocolate", "coral2", "chartreuse4", "cornflowerblue",
  "darkgoldenrod", "darkmagenta", "firebrick3", "yellow", "tan2",
  "seashell4", "peachpuff3", "navajowhite3", "azure4"
)

# Plot 
pf <-ggplot(physeq_family,aes(Sample,Abundance,fill=Phylum))+
  geom_bar(stat="identity", position = "fill")+
  facet_grid(.~Country, drop=TRUE,scale="free",space="free_x")
pf<-pf+scale_fill_manual(values = family_colors)  
pf<-pf+theme_bw()+ylab("Proportions")  
pf<-pf+ scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.margin = unit(0.3, "lines"))
pf<-pf+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
print(pf)
```

# Objective 3.1: Alpha diversity
```
rich = estimate_richness(physeq)
plot_richness(physeq, x = "Country", color = "Depth", measures = c("Observed", "Chao1", "Shannon", "Simpson"))+ geom_boxplot()
plot_richness(physeq, x = "Country", color = "Depth", measures = c("Observed", "Chao1", "Shannon", "Simpson"))


#Test whether the observed number of OTUs differs significantly between seasons. We make a non-parametric test, the Wilcoxon rank-sum test (Mann-Whitney):
pairwise.wilcox.test(rich$Observed, sample_data(physeq)$Country)
pairwise.wilcox.test(rich$Chao1, sample_data(physeq)$Country)
pairwise.wilcox.test(rich$ACE, sample_data(physeq)$Country)
pairwise.wilcox.test(rich$Shannon, sample_data(physeq)$Country)
pairwise.wilcox.test(rich$Simpson, sample_data(physeq)$Country)
pairwise.wilcox.test(rich$Fisher, sample_data(physeq)$Country)

#Export out range of alpha diversity
rich
write.table(rich, "alpha_diversity.csv")
```

# Objective 3.2: Beta diversity (PCoA)
```
#PCoA Plotting
ordinate(physeq, "PCoA", "bray") %>% 
  plot_ordination(physeq, ., color = "Country", title = "PCoA Bray-Curtis")+
  geom_point(size = 3) +   stat_ellipse()

#PERMANOVA - Distance based multivariate analysis of variance
# Function to run adonis test on a phyloseq object and a variable from metadata
phyloseq_to_adonis <- function(physeq, distmat = NULL, dist = "bray", formula) {
  
  if(!is.null(distmat)) {
    phydist <- distmat
  } else {
    phydist <- phyloseq::distance(physeq, dist)
  }
  
  metadata <- as(sample_data(physeq), "data.frame")
  
  # Adonis test
  f <- reformulate(formula, response = "phydist")
  adonis.test <- adonis(f, data = metadata)
  print(adonis.test)
  
  # Run homogeneity of dispersion test if there is only 1 variable
  if (grepl("\\+", formula)) {
    l <- list(
      dist = phydist, 
      formula = f, 
      adonis = adonis.test
    )
  } else {
    group <- metadata[ ,formula]
    beta <- betadisper(phydist, group)
    disper.test = permutest(beta)
    print(disper.test)
    
    l <- list(
      dist = phydist, 
      formula = f, 
      adonis = adonis.test, 
      disper = disper.test
    )
  }
  return (l)
}

#Run an adonis (PERMANOVA) test to see if they have different centroids
adonis.site <- phyloseq_to_adonis(
  physeq = physeq, 
  dist = "bray", 
  formula = "Country"
)
```

# Objective 4. Differential abundance
```
#DESeq2 conversion and call
diagdds = phyloseq_to_deseq2(physeq, ~ Country)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric") 

#if line 220 appeared as error, execute 225 and you will be able to execute it again (line 228)
#Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : 
#every gene contains at least one zero, cannot compute log geometric means

cts <- counts(diagdds) #Reference: https://support.bioconductor.org/p/72011/
geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
dds <- estimateSizeFactors(diagdds, geoMeans=geoMeans)
diagdds = DESeq(dds, test="Wald", fitType="parametric")

#Investigate test results table
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)

#Phylum Order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))

#Genus Order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

#Plotting the differential abundance using ggplot2
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=1) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
```

# Objective 5. Network analysis (igraph and Cytoscape)
## Load the library (install the package if you have not)
```
library(igraph) 
library(Hmisc)  
library(Matrix)   
library(RCy3)
```

## Read the input files (transpose if needed depending on your attribute correlations)
```
otu.table<-read.csv("SPE_pitlatrine.csv",row.names=1,check.names=FALSE)
otu.table<-t(otu.table)
dim(otu.table)
```

## Filtering out low abundant OTUs
```
otu.table.filter <- otu.table[ ,colSums(otu.table) >= 10]
print(c(ncol(otu.table),"versus",ncol(otu.table.filter))) #to compare before and after the filtering
```

## Calculating p-values, r-values, and set the threshold parameters for the network
```
#Calculating Spearman correlation
otu.cor <- rcorr(as.matrix(otu.table.filter), type="spearman")

#Calculating the p-value of correlation
otu.pval <- forceSymmetric(otu.cor$P)

#Filtering the correlation based on p-values and r-value
p.yes <- otu.pval<0.05
r.val = otu.cor$r

#write.table(r.val,"NetworkAnalysis.csv")

#Only select correlation values based on p-value criterion 
p.yes.r <- r.val*p.yes

#Selecting OTUs by correlation level
p.yes.r <- abs(p.yes.r)>0.7
p.yes.rr <- p.yes.r*r.val

#Creating adjacency matrix
adjm <- as.matrix(p.yes.rr)

#Creating a graph object from adjacency matrix
net.grph=graph.adjacency(adjm,mode="undirected",weighted=TRUE,diag=FALSE)

#To check the details of vertices and edges
#igraph::vcount(net.grph)
#igraph::ecount(net.grph)

#To remove isolated nodes
bad.vs<-V(net.grph)[degree(net.grph) == 0]

#Removing the isolated nodes from the graph object using the function delete.vertices()
net.grph <-delete.vertices(net.grph, bad.vs)
#igraph::vcount(net.grph)
#igraph::ecount(net.grph)

#To calculate the degree (number of interaction)
degAll <- igraph::degree(net.grph, v = igraph::V(net.grph), mode = "all")
#summary(net.grph)

net.grph <- igraph::set.vertex.attribute(net.grph, "degree", index = igraph::V(net.grph), value = degAll)

#To obtain the edge weight based on Spearman rank correlation
edgew<-E(net.grph)$weight
```

## Connect to Cytoscape (make sure to open the application before executing)
```
cytoscapePing()
```

## Feed the calculated network into Cytoscape
```
#Feed the calculated igprah object to Cytoscape
createNetworkFromIgraph(net.grph,"myIgraph")
```
# References
1. Phyloseq R Package for Objective 1-4

[Joey711 Github](http://joey711.github.io/phyloseq-demo/phyloseq-demo.html)

[Joey711 - DESEQ2](https://joey711.github.io/phyloseq-extensions/DESeq2.html)

[Deneflab Github](http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#libraries)

2. igraph R Package for Network Analysis

[Alcala & Poudel](https://kelseyandersen.github.io/NetworksPlantPathology/Microbiome_network_ICPP2018_v2.html)

[How to customize network analysis using Cytoscape](http://cytoscape.org/RCy3/articles/Network-functions-and-visualization.html)

# Datasets are available at:
Torondel, B., Ensink, J. H., Gundogdu, O., Ijaz, U. Z., Parkhill, J., Abdelahi, F., ... & Quince, C. (2016). Assessment of the influence of intrinsic environmental and geographical factors on the bacterial ecology of pit latrines. Microbial biotechnology, 9(2), 209-223.

[DOI:10.1111/1751-7915.12334](DOI:10.1111/1751-7915.12334)

[R Code for Ecological Data Analysis](http://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/ecological.html)
