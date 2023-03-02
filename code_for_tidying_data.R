library(tidyverse)
library(xlsx)

## BRINGING THE COUNTS MATRIX AND METADATA CLOSER TOGETHER

# Writing in the counts matrix
mRNA <- read.table("HMCL_mRNA.txt", header = TRUE, sep = "\t")
mRNA <- as_tibble(mRNA)

# Writing in the meta data
md <- read.xlsx("meta_data.xlsx", 1, endRow = 93)
md <- as_tibble(md)

# Removing duplicate cell line from mRNA matrix (Karpas929_ECACC_p10 and _p15 both available; delete _p10)
mRNA <- mRNA[,-23]

# Ensuring cell names in mRNA (i.e. colnames) match names in md (i.e md$Public.Name)
colnames(mRNA)[3:68] <- gsub("_p[0-9]*", "", colnames(mRNA)[3:68])
colnames(mRNA)[3:68] <- gsub("_PLB", "", colnames(mRNA)[3:68])

md$Keats_Lab_Name <- gsub("_PLB", "", md$Keats_Lab_Name)
md$Keats_Lab_Name <- gsub("JCRBsus", "JCRB_Sus", md$Keats_Lab_Name)
md$Keats_Lab_Name <- gsub("JCRBadh", "JCRB_Adh", md$Keats_Lab_Name)
md$Keats_Lab_Name <- gsub("RIKEN", "Riken", md$Keats_Lab_Name)

# Changing the order of the md rows (MMM1 features before MM1R in md)
md <- arrange(md, Keats_Lab_Name)

# Removing unsequenced cell lines from md
md <- filter(md, Keats_Lab_Name %in% colnames(mRNA))

# Ensuring that cell lines are the row names in md
md <- as.data.frame(md)
rownames(md) <- md$Keats_Lab_Name

# Arranging mRNA as a true counts matrix
mRNA <- as.data.frame(mRNA)
rownames(mRNA) <- mRNA[,1]
mRNA1 <- mRNA
mRNA <- mRNA[,-c(1,2)]

## TIDYING UP MD

md$Canonical_Translocations <- gsub(":", ";", md$Canonical_Translocations)
md_tibble <- as_tibble(md)

# Finding all the different translocations
tl <- unique(md_tibble$Canonical_Translocations)
tl <- gsub(":", ";", tl)
tl <- gsub("\\*\\*\\*", "", tl)
tl <- unlist(strsplit(tl, "\\+"))
tl <- gsub(" ", "", tl)
tl <- unique(tl)
tl <- tl[-9]

# Creating translocation columns in md

for(i in tl) md_tibble <- md_tibble %>% mutate("{i}" := ifelse(grepl(i, Canonical_Translocations, fixed = TRUE), 1, 0))

# Getting rid of extraneous columns

md_tibble <- select(md_tibble, !c(Canonical_Translocations, NA.:NA..10))
