#EMT ref
cos.sim <- function(scaled.matrix, input.mat) {
  dot.mat <- t(scaled.matrix) %*% input.mat
  ref.l2 <- apply(scaled.matrix, 2, function(xx) sqrt(sum(xx^2)))
  query.l2 <- apply(input.mat, 2, function(xx) sqrt(sum(xx^2)))
  dot.mat / (matrix(ref.l2, ncol = 1) %*% matrix(query.l2, nrow=1))
}


EMT <- function(input.mat, ref.db, expression.cutoff=0.25) {
  message(glue::glue("Loading {ref.db} ..."))
  seu <- qs::qread(ref.db)
  EMT_factor <- unique(seu$Treatment)
  input.mat <- as.matrix(input.mat)
  normalized.matrix <- as.matrix(seu[["RNA"]]@data) # [note] Seurat V5 conflict
  avg.expr <- rowMeans(normalized.matrix)
  genes.use <- names(avg.expr)[avg.expr > expression.cutoff]
  genes.use <- intersect(genes.use, rownames(input.mat))
  input.mat <- input.mat[genes.use, , drop=F]
  normalized.matrix <- normalized.matrix[genes.use, ]
  scaled.matrix <- t(scale(t(normalized.matrix)))
  message("Calculating cosine similarity ...")
  cos.mat <- cos.sim(scaled.matrix, input.mat)
  ## wilcox.test
  message("Calculating p value ...")
  res.list <- pbapply::pblapply(1:ncol(cos.mat), function(ii) {
    pbs.cells <- rownames(subset(seu@meta.data, Time == "0d"))
    s2 <- cos.mat[pbs.cells, ii]
    res <- lapply(EMT_factor, function(ct) {
      treatment.cells <- rownames(subset(seu@meta.data, Treatment == ct & Time != "0d"))
      s1 <- cos.mat[treatment.cells, ii]
      es <- mean(s1) - mean(s2)
      p.value <- wilcox.test(s1, s2)$p.value
      data.frame(
        cytokine = ct,
        es = es,
        p.value = p.value
      )
    })
    res <- do.call(rbind, res)
    rownames(res) <- EMT_factor
    res$p.adjust <- p.adjust(res$p.value, method = "fdr")
    res
  })
  names(res.list) <- colnames(cos.mat)
  return(res.list)
}

EMT_drug <- function(input.mat, ref.db, factor, expression.cutoff=0.25) {
  message(glue::glue("Loading {ref.db} ..."))
  seu <- qs::qread(ref.db)
  EMT_factor <- unique(seu$Drug)
  input.mat <- as.matrix(input.mat)
  normalized.matrix <- as.matrix(seu[["RNA"]]@data) # [note] Seurat V5 conflict
  avg.expr <- rowMeans(normalized.matrix)
  genes.use <- names(avg.expr)[avg.expr > expression.cutoff]
  genes.use <- intersect(genes.use, rownames(input.mat))
  input.mat <- input.mat[genes.use, , drop=F]
  normalized.matrix <- normalized.matrix[genes.use, ]
  scaled.matrix <- t(scale(t(normalized.matrix)))
  message("Calculating cosine similarity ...")
  cos.mat <- cos.sim(scaled.matrix, input.mat)
  ## wilcox.test
  message("Calculating p value ...")
  res.list <- pbapply::pblapply(1:ncol(cos.mat), function(ii) {
    # 首先检查PBS组细胞
    pbs.cells <- rownames(subset(seu@meta.data, Treatment == "Untreated"))
    
    # 确保PBS组有数据
    if(length(pbs.cells) == 0) {
      warning("No PBS cells found")
      return(NULL)
    }
    
    s2 <- cos.mat[pbs.cells, ii]
    
    res <- lapply(EMT_factor, function(ct) {
      # 检查TNF处理组细胞
      treatment.cells <- rownames(subset(seu@meta.data, Drug == ct & Treatment == factor))
      
      # 如果处理组或对照组没有足够数据，返回NA
      if(length(treatment.cells) == 0) {
        warning(paste("No cells found for", ct, "treatment group"))
        return(data.frame(
          drug = ct,
          es = NA,
          p.value = NA
        ))
      }
      
      # 提取对应的数据
      s1 <- cos.mat[treatment.cells, ii]
      
      # 检查是否有足够的非缺失数据
      s1 <- s1[!is.na(s1)]
      s2 <- s2[!is.na(s2)]
      
      # 如果仍然没有足够数据，返回NA
      if(length(s1) < 2 || length(s2) < 2) {
        warning(paste("Not enough observations for", ct))
        return(data.frame(
          drug = ct,
          es = NA,
          p.value = NA
        ))
      }
      
      # 执行Wilcoxon检验
      test_result <- tryCatch(
        wilcox.test(s1, s2),
        error = function(e) NULL
      )
      
      # 如果检验失败，返回NA
      if(is.null(test_result)) {
        return(data.frame(
          drug = ct,
          es = NA,
          p.value = NA
        ))
      }
      
      # 计算效应大小和p值
      es <- mean(s1) - mean(s2)
      p.value <- test_result$p.value
      
      data.frame(
        drug = ct,
        es = es,
        p.value = p.value
      )
    })
    
    # 合并结果
    res <- do.call(rbind, res)
    rownames(res) <- EMT_factor
    
    # 调整p值，忽略NA
    res$p.adjust <- p.adjust(res$p.value, method = "fdr")
    
    return(res)
  })
  
  names(res.list) <- colnames(cos.mat)
  return(res.list)
}

EMT.plot <- function(irea.res.df, color.breaks = c(-0.2,0,0.2)) {
  data.plot <- irea.res.df %>% arrange(p.value) %>% mutate(cytokine = factor(cytokine, levels=cytokine))
  abs.values <- color.breaks
  rel.values <- sapply(abs.values, function(xx) sum(data.plot$es < xx) / nrow(data.plot) )
  max.y <- max(-log10(data.plot$p.value))
  
  n.gap <- nrow(data.plot)
  delta.angle <- 360 / n.gap
  angles <- 90 + 0:(nrow(data.plot)-1) * -delta.angle
  hjust <- ifelse(1:n.gap > n.gap/2, 1, 0)
  angles[hjust == 1] <- angles[hjust == 1] + 180
  
  data.plot$is.receptor.exp <- FALSE
  data.plot$is.receptor.exp[sample(1:nrow(data.plot), 3)] = TRUE
  
  ggplot(data.plot, aes(cytokine, -log10(p.value))) + 
    geom_segment(aes(y=0, xend = cytokine, yend=max.y), color = "black", linewidth = 0.3, 
                 arrow = arrow(angle = 15, length = unit(0.1, "inches"), ends = , type = "closed")) + 
    geom_bar(aes(fill = es), stat="identity", color="white", size=.2) + 
    geom_text(aes(y = max.y*1.1, label = cytokine), angle = angles, hjust = hjust) + 
    # geom_tile(inherit.aes = F, data = subset(data.plot, is.receptor.exp),
    #           aes(cytokine, y = 1.05*max.y), fill = "grey", color = "black", linewidth = 0.5, width=0.8, height=0.9) +
    # geom_tile(inherit.aes = F, data = subset(data.plot, !is.receptor.exp),
    #           aes(cytokine, y = 1.05*max.y), fill = "white", color = "black", linewidth = 0.5, width=0.8, height=0.9) +
    scale_fill_gradientn(colors = c("blue", "white", "red"), values = rel.values) + 
    guides(fill = guide_colorbar(title = "ES")) + 
    coord_polar() + 
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank())
}

drug.plot <- function(irea.res.df, color.breaks = c(-0.2,0,0.2)) {
  data.plot <- irea.res.df %>% arrange(p.value) %>% mutate(drug = factor(drug, levels=drug)) %>% na.omit()
  abs.values <- color.breaks
  rel.values <- sapply(abs.values, function(xx) sum(data.plot$es < xx) / nrow(data.plot) )
  max.y <- max(-log10(data.plot$p.value))
  
  n.gap <- nrow(data.plot)
  delta.angle <- 360 / n.gap
  angles <- 90 + 0:(nrow(data.plot)-1) * -delta.angle
  hjust <- ifelse(1:n.gap > n.gap/2, 1, 0)
  angles[hjust == 1] <- angles[hjust == 1] + 180
  
  data.plot$is.receptor.exp <- FALSE
  data.plot$is.receptor.exp[sample(1:nrow(data.plot), 3)] = TRUE
  
  ggplot(data.plot, aes(drug, -log10(p.value))) + 
    geom_segment(aes(y=0, xend = drug, yend=max.y), color = "black", linewidth = 0.3, 
                 arrow = arrow(angle = 15, length = unit(0.1, "inches"), ends = , type = "closed")) + 
    geom_bar(aes(fill = es), stat="identity", color="white", size=.2) + 
    geom_text(aes(y = max.y*1.1, label = drug), angle = angles, hjust = hjust) + 
    # geom_tile(inherit.aes = F, data = subset(data.plot, is.receptor.exp),
    #           aes(drug, y = 1.05*max.y), fill = "grey", color = "black", linewidth = 0.5, width=0.8, height=0.9) +
    # geom_tile(inherit.aes = F, data = subset(data.plot, !is.receptor.exp),
    #           aes(drug, y = 1.05*max.y), fill = "white", color = "black", linewidth = 0.5, width=0.8, height=0.9) +
    scale_fill_gradientn(colors = c("blue", "white", "red"), values = rel.values) + 
    guides(fill = guide_colorbar(title = "ES")) + 
    coord_polar() + 
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank())
}
