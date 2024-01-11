# Breast Cancer Diagnosis

# This example aims to assess whether a lump in a breast could be malignant (cancerous) or benign (non-cancerous) from digitized images of a fine-needle aspiration biopsy.The breast cancer database was obtained from the University of Wisconsin Hospitals, Madison, from Dr William H. Wolberg. You can download the data from the UCI Machine learning database.
# 
# The Dataset
# 
# clump_thickness: (1-10). Benign cells tend to be grouped in monolayers, while cancerous cells are often grouped in multilayers. Higher clumpthickness associated with cancer
# 
# cell_size_uniformity: (1-10). Cancer cells tend to vary in size and shape. Loss of size uniformity associated with cancer
# 
# cell_shape_uniformity: (1-10). Uniformity of cell size/shape: Cancer cells tend to vary in size and shape. Loss of size uniformity associated with cancer
# 
# marginal_adhesion: (1-10). Normal cells tend to stick together. Cancer cells tend to lose this ability. So the loss of adhesion is a sign of malignancy.
# 
# single_epithelial_cell_size: (1-10). It is related to the uniformity mentioned above. Epithelial cells that are significantly enlarged may be a malignant cell.
# 
# bare_nuclei: (1-10). This is a term used for nuclei not surrounded by cytoplasm (the rest of the cell). Those are typically seen in benign tumors.
# 
# bland_chromatin: (1-10). Describes a uniform "texture" of the nucleus seen in benign cells. In cancer cells, the chromatin tends to be more coarse.
# 
# normal_nucleoli: (1-10). Nucleoli are small structures seen in the nucleus. In normal cells, the nucleolus is usually very small if visible at all. In cancer cells, the nucleoli become more prominent, and sometimes there are more of them.
# 
# mitoses: (1-10). Cancer is essentially a disease of uncontrolled mitosis. Higher number indicates cancer
# 
# classes: Benign (non-cancerous) or malignant (cancerous) lump in a breast.

# Importing the data an wrangling----

library(tidyverse) # for tidy data analysis
library(readr)     # for fast reading of input files

bc_data0 <-  read.csv(paste0("http://archive.ics.uci.edu/ml/machine-learning-databases/","breast-cancer-wisconsin/breast-cancer-wisconsin.data"), header = FALSE, stringsAsFactors = F)


names (bc_data0) <-  c("sample_code_number", 
                      "clump_thickness", 
                      "uniformity_of_cell_size", 
                      "uniformity_of_cell_shape", 
                      "marginal_adhesion", 
                      "single_epithelial_cell_size", 
                      "bare_nuclei", 
                      "bland_chromatin", 
                      "normal_nucleoli", 
                      "mitosis", 
                      "classes")

str(bc_data0)

# Notice how bare_nuclei does not follow the correct formatting bare_nuclei = as.numeric(bare_nuclei),
bc_data0$bare_nuclei = as.integer(bc_data0$bare_nuclei)
# NA introduced by coercion is just a warning message that lets us know that some values got coerced into NA

bc_data1 <- bc_data0 %>%
  dplyr::mutate(classes = ifelse(classes == "2", "benign",
                          ifelse(classes == "4", "malignant", NA)))

str(bc_data1)

# Exploratory Data Analysis----
summary(bc_data1) # Basic 

# Does each row repesent a unique patient record?
length(bc_data1$sample_code_number)
length(unique(bc_data1$sample_code_number))
# No! we want to filter out unique patients

bc_data2 <- bc_data1 %>% distinct(sample_code_number,.keep_all = TRUE)

# we also don't really care about sample code number, as it is not a biological variable.
row.names(bc_data2) <- bc_data2$sample_code_number

bc_data3 <- bc_data2 %>% select(-sample_code_number)

#Automated EDA
library(DataExplorer)
create_report(bc_data3)


# Missing Data---
# use only complete cases
# Missing values can be imputed we will see this later-- for now we will remove!
bc_data4 <- bc_data3[complete.cases(bc_data3),]
bc_pred <- bc_data4 %>% select(-classes)
bc_class <- bc_data4$classes


# Exploratory Data Analysis
ggplot(bc_data4, aes(x = classes, fill = classes)) + geom_bar()

# Notice that we have more benign than malignant-- there is a class imbalance.

# We will deal with unbalanced datasets later
#theme_set(theme_light(base_size = 18, base_family = "Poppins"))

g <- ggplot(stack(bc_pred), aes(x = ind, y = values)) + geom_boxplot() + 
  labs(title = "Boxplots of columns") + labs(x = "", y = "Values") + 
  scale_y_continuous(breaks = seq(1, 10, by = 1))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust=1),
    legend.position = "none",
    axis.title = element_text(size = 16),
    panel.grid = element_blank()
  )

g +
  geom_jitter(size = 2, alpha = 0.25, width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 5)



gather(bc_data4, x, y, clump_thickness:mitosis) %>%  # selecting data pairs
  ggplot(aes(x = y, color = classes, fill = classes)) +
  geom_density(alpha = 0.3) +
  facet_wrap( ~ x, scales = "free", ncol = 3)

# Create correlation matrices
# Benign
co_mat_benign <- filter(bc_data4, classes == "benign") %>%
  select(-classes) %>%
  cor()
# Malignant
co_mat_malignant <- filter(bc_data4, classes == "malignant") %>%
  select(-classes) %>%
  cor()

library(igraph)
g_benign <- graph.adjacency(co_mat_benign,
                            weighted = TRUE,
                            diag = FALSE,
                            mode = "upper")

g_malignant <- graph.adjacency(co_mat_malignant,
                               weighted = TRUE,
                               diag = FALSE,
                               mode = "upper")


# http://kateto.net/networks-r-igraph

cut.off_b <- mean(E(g_benign)$weight)
cut.off_m <- mean(E(g_malignant)$weight)

g_benign_2 <- delete_edges(g_benign, E(g_benign)[weight < cut.off_b])
g_malignant_2 <- delete_edges(g_malignant, E(g_malignant)[weight < cut.off_m])

c_g_benign_2 <- cluster_fast_greedy(g_benign_2) # implements network clustering methods 
c_g_malignant_2 <- cluster_fast_greedy(g_malignant_2) 

par(mfrow = c(1,2))

plot(c_g_benign_2, g_benign_2,
     vertex.size = colSums(co_mat_benign) * 10, # the larger the vertex/node the more correlated that vertex is with other features
     vertex.frame.color = NA, 
     vertex.label.color = "black", 
     vertex.label.cex = 0.8,
     edge.width = E(g_benign_2)$weight * 15,
     layout = layout_with_fr(g_benign_2),
     main = "Benign tumors")

plot(c_g_malignant_2, g_malignant_2,
     vertex.size = colSums(co_mat_malignant) * 10,
     vertex.frame.color = NA, 
     vertex.label.color = "black", 
     vertex.label.cex = 0.8,
     edge.width = E(g_malignant_2)$weight * 15,
     layout = layout_with_fr(g_malignant_2),
     main = "Malignant tumors")

# The nodes in the graph represent each feature, and edge between the two nodes indicates that the features are correlated


# Principal Component Analysis
library(ellipse)

# perform pca and extract scores
pcaOutput <- prcomp(as.matrix(select(bc_data4, -classes)), scale = TRUE, center = TRUE)
pcaOutput2 <- as.data.frame(pcaOutput$x)
PoV <- pcaOutput$sdev^2/sum(pcaOutput$sdev^2)
PoV
# define groups for plotting
pcaOutput2$groups <- bc_data4$classes

centroids <- aggregate(cbind(PC1, PC2) ~ groups, pcaOutput2, mean)

conf.rgn  <- do.call(rbind, lapply(unique(pcaOutput2$groups), function(t)
  data.frame(groups = as.character(t),
             ellipse(cov(pcaOutput2[pcaOutput2$groups == t, 1:2]),
                     centre = as.matrix(centroids[centroids$groups == t, 2:3]),
                     level = 0.95),
             stringsAsFactors = FALSE)))

g1 <- ggplot(data = pcaOutput2, aes(x = PC1, y = PC2, group = groups, color = groups)) + 
  geom_polygon(data = conf.rgn, aes(fill = groups), alpha = 0.2) +
  geom_point(size = 2, alpha = 0.6) + 
  labs(color = "",
       fill = "") 

g1
# Can also make output more informative
g2 <- g1+
  labs(x = paste0("PC1: ", round(PoV[1], digits = 2) * 100, "% variance"),
       y = paste0("PC2: ", round(PoV[2], digits = 2) * 100, "% variance"))

g2


# There is a package to help with ggplot grammar
install.packages("esquisse")
esquisse::esquisser()

# Basic summary stats
psych::describeBy(bc_data4)
psych::describeBy(bc_pred, bc_class)

# a better looking table
arsenal::tableby(classes~., data = bc_data4, total= TRUE) %>% summary(text = TRUE)
