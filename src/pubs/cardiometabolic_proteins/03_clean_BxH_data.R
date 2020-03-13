library(impute)
library(data.table)

data.dir <- "data/BxH/raw/"

# Load in gene expression
exprSet <- list(
  adipose = read.table(
    file.path(data.dir, "adipose_gene_expression.txt"), header=T, 
    sep="\t", row.names=1
  ),  
  brain = read.table(
    file.path(data.dir, "brain_gene_expression.txt"), header=T, 
    sep="\t", row.names=1
  ),  
  liver = read.table(
    file.path(data.dir, "liver_gene_expression.txt"), header=T, 
    sep="\t", row.names=1
  ),  
  muscle = read.table(
    file.path(data.dir, "muscle_gene_expression.txt"), header=T, 
    sep="\t", row.names=1
  )
)

# Load in the probe annotations for the custom Agilent array
geneInfo <- read.table(
  file.path(data.dir, "array_annotation.txt"), sep="\t", stringsAsFactors=FALSE,
  header=TRUE
)
geneInfo <- as.data.table(geneInfo)
# map column names to old reporter file so scripts don't break.
setnames(geneInfo, "Probe.ID", "reporter_id")
setnames(geneInfo, "Entrez.ID", "entrez_id")
setnames(geneInfo, "Gene.Symbol", "gene_sym")
setkey(geneInfo, reporter_id)

# Load in the phenotype data
phen.desc <- fread(file.path(data.dir, "trait_description.txt"))
phen.desc[, trait_id := as.character(trait_id)]
setkey(phen.desc, trait_id)

phen <- read.table(file.path(data.dir, "phenotype.txt"))
phen <- t(phen)
phen <- as.data.table(phen, keep.rownames=TRUE)
setnames(phen, c("id", phen.desc[colnames(phen)[-1], trait_name]))


#----------------------------
# Quality control
#----------------------------

# Fix class
exprSet <- lapply(exprSet, as.matrix)

# Impute remaining missing data
exprSet <- lapply(exprSet, impute.knn)

# Add a column to the phenotype data for sex
sex <- fread(file.path(data.dir, "gender_info.txt"))
setnames(sex, c("id", "sex"))
sex[, id := paste0("X", id)]
sex[, sex := factor(sex, levels=c("MALE", "FEMALE"))]
phen <- data.table:::merge.data.table(phen, sex, by="id")

# now remove traits with un-normalised equivalents
normed <- colnames(phen)[grep("(_log)|(_sqrt)", colnames(phen))]
non_normed <- gsub("(_log)|(_sqrt)", "", normed)
non_normed <- gsub("MCP-1_phys", "MCP-1", non_normed)
non_normed_cols <- match(non_normed, colnames(phen))
phen <- phen[,-non_normed_cols, with=FALSE]
phen.desc <- phen.desc[!(trait_name %in% non_normed)]

# Add long form names
phen.names <- c(
  "Length", "Weight", "Abdominal fat", "Other fat", "Total fat", 
  "100 x fat / weight", "Total cholesterol", "Unesterified cholesterol", 
  "Free fatty acids", "Glucose", "LDL+VLDL", "Aortic Lesion", "Insulin",
  "MCP-1 (CCL2)", "Triglyceride", "HDL", "Adiponectin", "Leptin", 
  "Aneurysm severity", "Medial aortic calcification", "Lateral aortic Calfication",
  "Bone density (all limbs)", "Bone density (femur)", "Glucose / Insulin",
  "HDL / LDL+VLDL", "Leptin / Weight", "Leptin / Abdominal fat", 
  "Leptin / Total fat", "Insulin / Weight", "Insulin / Abdominal Fat",
  "Insulin / Total fat", "Apolipoprotein A1", "Apolipoprotein A2"
)
names(phen.names) <- colnames(phen)[-c(1, ncol(phen))]
phen.desc[,long_name := phen.names[trait_name]] 

# Remove redundant columns
phen.desc[, timepoint := NULL]
phen.desc[, adjustment := NULL]
phen.desc[,trait_id := NULL]

# Fix some descriptions
phen.desc[transform == "normal", transform := "none"]

# Add missingness
phen.desc[, missingness := phen[,sapply(.SD, function(x) sum(is.na(x))/.N)][trait_name]]

# Add non-missing N
phen.desc[, sample_size := phen[,sapply(.SD, function(x) sum(!is.na(x)))][trait_name]]

#--------------------------
# Write out cleaned datasets
#--------------------------

for (ii in seq_along(exprSet)) {
  dt <- as.data.table(exprSet[[ii]]$data, keep.rownames="array_id")
  fwrite(dt, sep="\t", quote=FALSE, file=sprintf("data/BxH/cleaned/%s_expression.txt", names(exprSet)[[ii]]))
}

fwrite(phen, sep="\t", quote=FALSE, file="data/BxH/cleaned/phenotype.txt")
fwrite(phen.desc, sep="\t", quote=FALSE, file="data/BxH/cleaned/phenotype_info.txt")
fwrite(geneInfo, sep="\t", quote=FALSE, file="data/BxH/cleaned/expression_info.txt")
