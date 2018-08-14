otus = matrix(sample(0:100, 100, replace = TRUE), nrow = 10, ncol = 10)
otus[9,] <- sample(0:1, 10, replace = TRUE)
otus[3,] <- sample(0:2, 10, replace = TRUE)
otus[2,] <- sample(0:2, 10, replace = TRUE)
otus[5,] <- 0
otus[5,8] <- 100
otus[5,4] <- 1
otus[5,2] <- 2
rownames(otus) <- paste0("OTU", 1:nrow(otus)) 
colnames(otus) <- paste0("S", 1:ncol(otus))
taxmat = matrix(sample(letters, 70, replace = TRUE), 
                nrow = nrow(otus), 
                ncol = 7)
rownames(taxmat) <- rownames(otus)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
OTU = otu_table(otus, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
physeq = phyloseq(OTU, TAX)
print(physeq)


# suma de otus dentro de UNA de las muestras
colSums(otus)
# suma de UN otu a lo largo de las muestras
rowSums(otus) # 
# para tratar los k-tones, (sea caso de singletones) esperas que por minimo, rowSums(otus) >= phyloseq::nsamples(physeq) +1
# trabajamos sobre los rowsums (taxa_sums) como a continuacion:
# condicionamos que almenos algun-solo otu se encuentre dos veces en cada uno de las muestras. de este modo tratamos

otu_table(prune_taxa(
  taxa_sums(physeq) >= phyloseq::nsamples(physeq) +1, # aqui podria ir algun valor estadistico! varianza??
  physeq))

# change to percent

physeq <- transform_sample_counts(physeq, function(x) (x / sum (x) )) 
physeq <- transform_sample_counts(physeq, function(x) x * 100 ) 

otu_table(filter_taxa(physeq, function(x) sum(x > 2) > (0.1*length(x)), TRUE))
otu_table(filter_taxa(physeq, function(x) x > 2, TRUE))


f1 <- filterfun_sample(topk(2))
wh1 <- genefilter_sample(physeq, f1, A=2)
otu_table(prune_taxa(wh1,physeq))


### top10

f1 <- filterfun_sample(topp(0.1))
wh1 <- genefilter_sample(Animalia, f1, A=2)
top10 <- prune_taxa(wh1,Animalia)
top10 <- transform_sample_counts(top10, function(x) (x / sum (x) ))
top10 <- transform_sample_counts(top10, function(x) x * 100 )

Animalia_glom <-tax_glom(top10, taxrank="Family")

h <- plot_heatmap(top10, low="#66CCFF", high="#000033",
                  taxa.label="Phylum",
                  method = "NMDS", distance = "bray"
)

h + facet_grid(Class ~ ., scales = "free", space = "free") 


# por tanto removamos esos singletones

Ap <- prune_taxa(
  taxa_sums(Animalia) >= phyloseq::nsamples(Animalia) +1, 
  Animalia)

y2 <- plyr::ldply(taxa_sums(Ap), rbind, use.names = TRUE)

Animalia_1 <- otu_table(prune_taxa(taxa_sums(Animalia) >= 2, Animalia))

# 1) remove k-tones criteria OTUsize > 2 & sumOTUlocation > 4 # vemos otus que estan una sola vez a lo largo de las muestras, singletones
# ie. remove doubletones and remove an OTU that does not appear in at least two samples 
# 2) remove NA into  y rank ??? ex. phylum <- subset_taxa(phyloseq, Phylum!="NA")
# 3) .... 

# singletones
# otu_table(prune_taxa(taxa_sums(Animalia) <= 1, Animalia))

# Remove taxa not seen more than 2 times in at least 20% of the samples. This protects against an OTU with small mean & trivially large C.V.
# 

Aprune <- filter_taxa(Animalia, function(x) sum(x)  > 2, TRUE)
Aprune <- filter_taxa(Animalia, function(x) sum(x > 3) > (1*length(x)), TRUE)


# function to find the most abundant taxa. author: smajor; in https://github.com/joey711/phyloseq/issues/847
# Goes through a phyloseq object, picks out the most abundant taxa and gives the abundance for each
# and identifies which taxa is most abundant for which sample
# 

top <- find.top.taxa(rr,"Phylum")


Animalia <- subset_taxa(phyloseq, Kingdom == "Animalia")
Family <-tax_glom(Animalia, taxrank="Family")
#
taxaSums <- plyr::ldply(taxa_sums(Family), rbind)
colnames(taxaSums) <- c("Taxa", "Size")
taxaSums <- taxaSums[order(-taxaSums[,2]),]
#
sampleSums <- plyr::ldply(sample_sums(Family), rbind)
colnames(sampleSums) <- c("Sample", "Abundance")
sampleSums <- sampleSums[order(-sampleSums[,2]),]
#


find.top.taxa <- function(x, taxa){
  keepTaxa = apply(X = as(otu_table(x), "matrix") > 2L,
                   MARGIN = 1, FUN = sum) >= 2L ## Remove OTUs not k greater than k (2L)
  phy = prune_taxa(keepTaxa, x) # get abundance in %
  phy <- transform_sample_counts(phy, function(y) y/sum(y)) # agglomerate taxa
  glom <- tax_glom(phy, taxrank = taxa) # create dataframe from phyloseq object
  dat <- psmelt(glom) # convert Phylum to a character vector from a factor because R
  dat[,taxa] <- as.character(dat[,taxa]) # group dataframe by Phylum, calculate median rel. abundance
  dat[, median := median(Abundance, na.rm = TRUE), 
      by = taxa] # Change name to remainder of Phylum less than 1%
  dat[(median <= 0.01), taxa := "Remainder"]
  return(dat)
}


find.taxa <- function(x,taxa, find){
  # select find.taxa(phyloseq-obj, "Rank-level", find either, "max" or "min" Abundance)
  # Ex. find.taxa(phyloseq, "Phylum", "min")
  wfind <- as.name(paste0("which.",find))
  require(phyloseq)
  x <- subset_taxa(x, Kingdom == "Animalia")
  top.taxa <- tax_glom(x, taxa)
  otu <- otu_table(t(top.taxa)) # remove the transformation if using a merge_sample object
  tax <- tax_table(top.taxa)
  j<-apply(otu,1, wfind)
  k <- j[!duplicated(j)]
  l <- as.matrix(tax[k,])
  m <- as.matrix(otu[,k])
  s <- as.name(taxa)
  colnames(m) = l[,taxa]
  n <- colnames(m)[apply(m,1, wfind)]
  m <- as.data.frame(cbind(m,n))
  return(m)
}


#xglom <- subset_taxa(Animalia, Family!="NA")
#zglom <- tax_glom(xglom, taxrank = 'Family')


f1 <- filterfun_sample(topp(0.1)) # fraction of the most abundant taxa to be kept
wh1 <- genefilter_sample(Family, f1, A=42) # in all of the samples in which a taxa/OTU passed the filter
ts <- prune_taxa(wh1, Family)
ts
dp <- psmelt(ts)

ggrare2 <- function (physeq_object, step = 10, label = NULL, color = NULL, 
          plot = TRUE, parallel = FALSE, se = TRUE) 
{
  x <- methods::as(phyloseq::otu_table(physeq_object), "matrix")
  if (phyloseq::taxa_are_rows(physeq_object)) {
    x <- t(x)
  }
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  rarefun <- function(i) {
    #cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- vegan::rarefy(x[i, , drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    }
    else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- parallel::mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  }
  else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)
  if (!is.null(phyloseq::sample_data(physeq_object, FALSE))) {
    sdf <- methods::as(phyloseq::sample_data(physeq_object), 
                       "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  if (length(color) > 1) {
    data$color <- color
    names(data)[names(data) == "color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  if (length(label) > 1) {
    labels$label <- label
    names(labels)[names(labels) == "label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  p <- ggplot2::ggplot(data = data, ggplot2::aes_string(x = "Size", 
                                                        y = ".S", group = "Sample", color = color))
  p <- p + ggplot2::labs(x = "Sequence Sample Size", y = "Species Richness")
  if (!is.null(label)) {
    p <- p + ggplot2::geom_text(data = labels, ggplot2::aes_string(x = "x", 
                                                                   y = "y", label = label, color = color), size = 4, 
                                hjust = 0)
  }
  p <- p + ggplot2::geom_line()
  if (se) {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes_string(ymin = ".S - .se", 
                                                      ymax = ".S + .se", color = NULL, fill = color), alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}

rrplot <- ggrare2(rr, step = 100, label = "Estacion", 
                 color="Transecto", 
                 se = FALSE, plot = FALSE) 

