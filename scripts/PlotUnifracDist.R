library("ape")
library(ggfortify)
library(ggrepel)
library(fs)


file_list <- dir_ls("MAGsGeneration/Results/UnifracDist", regexp = "\\.diss$")


#ggbiplot(fit, labels =  rownames(USArrests))par(mfrow = c(4, 4))

for (i in seq_along(file_list)) {
  filename <- file_list[[i]]
  a <- paste("tail -n 1 ", filename, " | wc -w")

  num_cols <- as.numeric(system(a, intern = TRUE))

  mat <- read.table(
    filename,
    fill = TRUE,
    col.names = rep('', num_cols + 1),
    row.names = 1
  )
  mat <- mat[-1, ]
  rownames(mat) <- gsub("_1", "", rownames(mat))
  distUni <- as.dist(as(mat, "matrix"))

  res <- pcoa(distUni)
  map1 <- map[match(rownames(mat), map$Samples), ]
  # rownames(map1) <- map1$GoodSampleID

  list <- unique(map1$Category)
  names(colors) <- list
  filenameOut <- gsub(".diss", ".pdf", filename)

  autoplot(prcomp(distUni, scale. = TRUE), label = FALSE, main = filename) +
    geom_text_repel(
      label = labels(distUni),
      size = 3,
      aes(color = map1$Category)
    ) +
    geom_point() +
    theme_bw()
  ggsave(filenameOut, width = 10, height = 7)
}
