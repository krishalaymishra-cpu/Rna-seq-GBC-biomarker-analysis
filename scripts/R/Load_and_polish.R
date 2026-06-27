#THE FEATURE COUNT GENERATED WAS LOADED IN R

count<-read.delim("final_feature_count.txt",header = FALSE,sep = "\t")
View(count)
count$V2=NULL
count$V3=NULL
count$V4=NULL
count$V5=NULL
count$V6=NULL
count<-count[-1,(1:41)]

new_colnames <- as.character(unlist(count[1, ]))
colnames(count) <- new_colnames

# Remove that first row from the data
count <- count[-1, ]

# Optional: strip BAM paths down to SRR IDs
colnames(count) <- gsub(".*/(SRR[0-9]+).*", "\\1", colnames(count))

# Check result
head(colnames(count))


#sorting columns
library(gtools)

# Reorder columns 2:41 by SRR ID (numeric ascending)
sorted_cols <- mixedsort(colnames(count)[2:41])
count <- count[, c("Geneid", sorted_cols)]


###MAPPING THE ENTRENZ IDS WITH ITS GENE SYMBOLS
library(biomaRt)

# Load your count matrix file
#count_data <- read.csv("new count(202479).csv", stringsAsFactors = FALSE)

# Extract original gene identifiers
original_ids <- as.character(count[[1]])

# Connect to Ensembl Biomart (Human)
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", 
                   dataset = "hsapiens_gene_ensembl", 
                   host = "https://useast.ensembl.org")

# Entrez IDs
mapping1 <- getBM(
  attributes = c("entrezgene_id", "hgnc_symbol"),
  filters = "entrezgene_id",
  values = original_ids,
  mart = ensembl
)
colnames(mapping1) <- c("original_id", "gene_symbol")

# Ensembl Gene IDs
mapping2 <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = original_ids,
  mart = ensembl
)
colnames(mapping2) <- c("original_id", "gene_symbol")

# Already HGNC symbols
mapping3 <- getBM(
  attributes = c("hgnc_symbol"),
  filters = "hgnc_symbol",
  values = original_ids,
  mart = ensembl
)
colnames(mapping3) <- "gene_symbol"
mapping3$original_id <- mapping3$gene_symbol
mapping3 <- mapping3[, c("original_id", "gene_symbol")]

# Combine all mappings
all_mappings <- unique(rbind(mapping1, mapping2, mapping3))
# Replace column A name to "original_id" for merging
colnames(count)[1] <- "original_id"

# Merge mapping with expression data
merged <- merge(all_mappings, count, by = "original_id", all.y = TRUE)

# Replace blank gene symbols with NA
merged$gene_symbol[merged$gene_symbol == ""] <- NA

# Remove rows where gene_symbol is NA
filtered <- merged[!is.na(merged$gene_symbol), ]

# Final output with gene_symbol as first column
final_count <- filtered[, c("gene_symbol", setdiff(colnames(filtered), c("gene_symbol", "original_id")))]
View(final_count)

final_count[,1] <- make.unique(as.character(final_count[,1]))
#shifting the gene symbol as rownames
rownames(final_count) <- final_count[,1]
count <- final_count[, -1]
View(count)
str(count)


##FIlter the noise

# Load the library
library(edgeR)

dge<-DGEList(counts = count)
library(zoo)
count[]<-t(na.aggregate(t(count)))
count
#filtering lowly expressed genes
#keeping genes with counts per million above 1 in atleast 2 samples
keep_cpm <- filterByExpr(dge)
keep_raw <- rowSums(dge$counts) >= 10
keep <- keep_cpm & keep_raw

dge_filtered <- dge[keep,, keep.lib.sizes = FALSE]

#recalculate lib size
dge_filtered<-calcNormFactors(dge_filtered)
filtered_counts<-dge_filtered
# Add gene symbols back as a column
filtered_counts <- data.frame(GeneSymbol = rownames(filtered_counts), filtered_counts)
View(filtered_counts)
rownames(filtered_counts) <- filtered_counts[,1]
count <- filtered_counts[, -1]
View(count)


###REPLACING SRR WITH ITS GSM

rownames(count)=count[,1]
count$X=NULL
colnames(count)[1:40]=paste0("GSM40995",32:71)
