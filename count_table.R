ct.dir <- c("/Users/cigom/metagenomics/phyloseq_in/count_table")
ct.files <- list.files(ct.dir, full.names = FALSE)


# read the list of count tables
count.tbl <- lapply(ct.files, read.csv, sep = "\t")

# ..then sum the columns per sample
list <- list()
for (i in 1:length(ct.files)) {
    list[[i]] <- colSums(count.tbl[[i]][-c(1:2)])
}

# tidy the object list and remove duplicate values:
df <- do.call(rbind, list)
df <- t(df[!duplicated(df), ])


df_names <- strsplit(rownames(df), "_")
df_names <- sapply(df_names, "[", c(3))

init.contigs <- read.csv("../initial.contigs", sep = "\t" )

plotdat <- cbind(df_names, init.contigs, df)
colnames(plotdat) <- c("Estacion", "Total", "Buena Calidad", "Alineamiento", "Clasificacion")
rownames(plotdat) <- df_names

plotdat <- reshape::melt(plotdat)


# And the plot:
library(ggpubr)

ggbarplot(plotdat, x = "Estacion", y = "value",
          fill = "variable",           # change fill color by mpg_level
          color = "variable",            # Set bar border colors to white
          palette = "Paired",            # jco journal color palett. see ?ggpar
          x.text.angle = 90,          # Rotate vertically x axis texts
          merge = TRUE,
          yscale = "none",
          ylab = "Numero de secuencias",
          legend.title = "Secuencias",
          rotate = TRUE,
          ggtheme = theme_minimal()
          #position = position_fill(vjust = 0.05)
          )

R1 = 5, 639, 021
R2 = 5, 643, 518

Total Contigs: 5,630,570 <- averiguar por muestra!
Good Quality Contigs (gqc):5,032,850 <- Buena calidad 
gqc Aligned (barcodes):5,024,954 <-- Alineamiento
non-quimeric barcodes:5,017,815 <-- clasificados
gqc euk-Unique (gqcU):207,349
gqcU aligned (barcodes):204,011
Barcodes made unique:139,527
Pre-clustering barcodes:44,242
unique non-quimeric barcodes:42,138
OTUs with 0.03 % similitud:19,263
Phylotypes:1349
Single-tones:10329
Double-tones:3173
Tripe-tones:1498
Tetratones:831

