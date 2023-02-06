
plot_violin_bar_srt <- function(input, title = "", gene, color_by, log_scale = F, colors = NULL, 
                                facet_by = NULL, spread = NULL, jitter_pts = T, plot_mean = T, 
                                size = 1, sig = 3, number_labels = T, text_sizes = c(20, 10, 7, 15, 10, 10), alpha = 0.5, theme = "classic", plot_legend = TRUE, patient_meta_col = "Patient", plot_patient_pos_frac = TRUE) {
  if (class(input) == "ExpressionSet") {
    df <- pData(input)[, colnames(pData(input)) %in% c(gene, 
                                                       color_by, facet_by), drop = F]
    df <- cbind(df, raw = exprs(input)[gene, ])
    colnames(df) <- gsub("-", "", colnames(df))
    gene <- gsub("-", "", gene)
    
    colnames(df)[which(colnames(df) == "raw")] <- paste(gene)
    
    if (any(!is.null(spread))) {
      others <- setdiff(unique(df[, spread[1]]), spread[2])
      ind <- which(df[, spread[1]] == spread[2])
      rmdf <- df[ind, ]
      df <- df[-ind, ]
      for (i in 1:length(others)) {
        rmdf[, spread[1]] <- others[i]
        df <- rbind(df, rmdf)
      }
    }
  }
  
  if (class(input) == "Seurat") {
    require(stringr)
    
    meta <- input@meta.data
    
    if (any(Assays(input) == "SCT") & nrow(input@assays$SCT@data) > 0) {
      normalized_cpm <- t(input@assays$SCT@data)
    }
    
    #Use RNA assay normalized counts if SCT-transform not used:
    if (eval(any(Assays(input) == "SCT") & nrow(input@assays$SCT@data) > 0) == FALSE) {
      normalized_cpm <- t(input@assays$RNA@data)
    }
    
    pca <- as.data.frame(input@reductions$pca@cell.embeddings)
    
    meta <- cbind(meta, pca)
    
    #find all factor columns:
    factor.columns <- which(sapply(meta, FUN = class) == "factor")
    meta.data.factors <- meta[,factor.columns]
    
    #Temporarily convert all to characters:
    meta[,factor.columns] <- sapply(meta[,factor.columns], FUN = as.character)
    
    if (any(str_detect(string = colnames(meta), pattern = gene)) == TRUE ) {
      
      df <- meta[,which(colnames(meta) %in% c(gene, color_by, facet_by, patient_meta_col))] 
      
    }
    
    if (any(str_detect(string = colnames(normalized_cpm), pattern = gene)) == TRUE ) {
      
      df <- as.data.frame(cbind(meta[,which(colnames(meta) %in% c(gene, color_by, facet_by, patient_meta_col))], normalized_cpm[,which(colnames(normalized_cpm) %in% c(gene, color_by, facet_by))]))
      
      colnames(df)[1:length(colnames(meta)[which(colnames(meta) %in% c(color_by, facet_by, patient_meta_col))])] <- colnames(meta)[which(colnames(meta) %in% c(color_by, facet_by,patient_meta_col))]
      colnames(df)[-which(colnames(df) %in% c(color_by, facet_by, patient_meta_col))] <- paste(gene)
      df[,-which(colnames(df) %in% c(color_by, facet_by, patient_meta_col))] <- sapply(df[,-which(colnames(df) %in% c(color_by, facet_by, patient_meta_col))], FUN = as.numeric)
      
      # Convert plotting df to factors where needed:
      if (any(colnames(df) %in% colnames(meta[,factor.columns]))) {
        
        columns.to.convert <- colnames(df)[which(colnames(df) %in% colnames(meta[,factor.columns]))]
        
        for (i in 1:length(columns.to.convert)) {
          col.name <- columns.to.convert[i]
          
          df[,which(colnames(df) == col.name)] <- factor(as.character(df[,which(colnames(df) == col.name)]),
                                                         levels = levels(meta.data.factors[,which(colnames(meta.data.factors) == col.name)]))
        }
      }
    }
    
  }
  
  #Potential Log Scale:
  if (log_scale == T) {
    df <- df %>% mutate_at(.vars = colnames(df)[which(colnames(df) %in% gene)], function(x) log2(x+1))
  } 
  
  
  #Colors:
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  color_by_vals <- df %>% dplyr::select(color_by) %>% unlist() %>% as.character %>% unique()
  cols <- gg_color_hue(length(unique(color_by_vals)))
  
  #Plot:
  g <- ggplot(df)
  if (all(!is.null(colors))) {
    g <- g + scale_color_manual(values = c(colors))
    g <- g + scale_fill_manual(values = c(colors))
  }
  if (theme == "bw") {
    g <- g + theme_bw()
  }
  else {
    g <- g + theme_classic()
  }
  if (title == "") {
    title <- gene
  }
  g <- g + labs(title = title, y = "Individual CPM")
  
  if (log_scale == T) {
    g <- g + labs(title = title, y = "Individual log2 CPM")
  } 
  g <- g + theme(plot.title = element_text(size = text_sizes[1]), 
                 axis.title = element_text(size = text_sizes[2]), axis.text = element_text(size = text_sizes[3]), 
                 legend.title = element_text(size = text_sizes[4]), legend.text = element_text(size = text_sizes[5]))
  g <- g + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5), 
                 axis.title.x = element_blank(), axis.text.x = element_blank(), 
                 axis.ticks.x = element_blank())
  
  if (jitter_pts == T) {
    #note: always plot zeros first because otherwise it might drop factor levels by plotting nonzero dots first.
    g <- g + geom_jitter(data = df[which(unlist(df[,gene]) == 0),], aes_string(x = color_by, y = gene, 
                                                                               col = color_by), width = 0.2, height = 0, size = size)
    
    g <- g + geom_jitter(data = df[which(unlist(df[,gene]) > 0),], aes_string(x = color_by, y = gene, 
                                                                              col = color_by), width = 0.2, height = 0.05*(max(df[,gene])-min(df[,gene])) ,size = size)
    g <- g + geom_violin(aes_string(x = color_by, y = gene, 
                                    fill = color_by), col = "black", trim = T, scale = "width", 
                         alpha = alpha)
    #plot mean positive expression line:  
    g <- g + stat_summary(fun = function(x) mean(x[which(x >0)]), color = "black", geom ="crossbar", aes_string(x = color_by, y = gene), size = 0.5, width = 0.5)
  }
  
  if (length(facet_by) == 1) {
    g <- g + facet_grid(facets = reformulate(facet_by))
  } else if (length(facet_by) == 2) {
    stop("Ever since I added positive fraction barplots to this function, it can only facet by 1 variable at a time.")
  } else if (length(facet_by) > 2) {
    stop("Parameter facet_by needs to be a string with equal or less than two variables.")
  }
  
  #Store legend to plot separately at bottom:
  legend <- ggpubr::as_ggplot(ggpubr::get_legend(g))
  g <- g + theme(legend.position = "none")
  
  #Barplot of positive fractions and cell count strip text:
  if (is.null(facet_by)) {
    df.stats <- df %>% group_by_at(color_by) %>% summarize(cellcount = n(),
                                                           positive_cells = sum(eval(as.symbol(gene))>0))
    
    df.stats <- df.stats %>% mutate(positive_frac = round(positive_cells/cellcount, digits = 3))
    
    df.patient.stats <- df %>% group_by_at(c(color_by, patient_meta_col)) %>% summarize(cellcount = n(),
                                                                                        positive_cells = sum(eval(as.symbol(gene))>0))
    
    df.patient.stats <- df.patient.stats %>% mutate(positive_frac = round(positive_cells/cellcount, digits = 3))
    
    f <- ggplot(df.stats, aes(x = eval(as.symbol(color_by)), y = positive_frac, fill = eval(as.symbol(color_by)))) + 
      theme_classic() +
      geom_col()
    
    if (plot_patient_pos_frac == TRUE) {
      f <- f + geom_jitter(data = df.patient.stats[which(df.patient.stats$positive_frac == 0),], aes(x = eval(as.symbol(color_by)), y = positive_frac), width = 0.2, height = 0, size = 0.75, alpha = 0.5)
      f <- f + geom_jitter(data = df.patient.stats[which(df.patient.stats$positive_frac > 0),], aes(x = eval(as.symbol(color_by)), y = positive_frac), width = 0.1,size = 0.75, alpha = 0.5)
    }
    
    if (all(!is.null(colors))) {
      f <- f + scale_color_manual(values = c(colors))
      f <- f + scale_fill_manual(values = c(colors))
    }
    
    f <- f + labs(x = paste(color_by), y = "Positive Fraction")
    f <- f + theme(legend.position = "none",
                   axis.title.x = element_blank())
    
    
    c <- ggpubr::ggsummarytable(data = df.stats, x = paste(color_by) , y = "cellcount", size = 3) + ggpubr::clean_table_theme()
    
  }
  
  if (length(facet_by) == 1) {
    df.stats <- df %>% group_by_at(c(color_by, facet_by)) %>% summarize(cellcount = n(),
                                                                        positive_cells = sum(eval(as.symbol(gene))>0))
    
    df.stats <- df.stats %>% mutate(positive_frac = round(positive_cells/cellcount, digits = 3))
    
    df.patient.stats <- df %>% group_by_at(c(color_by, facet_by, patient_meta_col)) %>% summarize(cellcount = n(),
                                                                                                  positive_cells = sum(eval(as.symbol(gene))>0))
    
    df.patient.stats <- df.patient.stats %>% mutate(positive_frac = round(positive_cells/cellcount, digits = 3))
    
    
    f <- ggplot(df.stats, aes(x = eval(as.symbol(color_by)), y = positive_frac, fill = eval(as.symbol(color_by)))) + 
      theme_classic() +
      geom_col() +
      facet_grid(facets = reformulate(facet_by))    
    if (plot_patient_pos_frac == TRUE) {
      f <- f + geom_jitter(data = df.patient.stats[which(df.patient.stats$positive_frac == 0),], aes(x = eval(as.symbol(color_by)), y = positive_frac), width = 0.2, height = 0, size = 0.75, alpha = 0.5)
      f <- f + geom_jitter(data = df.patient.stats[which(df.patient.stats$positive_frac > 0),], aes(x = eval(as.symbol(color_by)), y = positive_frac), width = 0.1,size = 0.75, alpha = 0.5)
    }
    
    if (all(!is.null(colors))) {
      f <- f + scale_color_manual(values = c(colors))
      f <- f + scale_fill_manual(values = c(colors))
    }
    
    f <- f + labs(x = paste(color_by), y = "Positive Fraction")
    f <- f + theme(legend.position = "none", 
                   axis.text.x = element_blank(),
                   axis.title.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   strip.text.x = element_blank())
    
    c <- ggpubr::ggsummarytable(data = df.stats, x = paste(color_by), color = paste(color_by),y = "cellcount", size = 2) + ggpubr::clean_table_theme() + theme(legend.position = "none") +
      facet_grid(facets = reformulate(facet_by)) + theme(strip.text.x = element_blank())
  }
  
  
  if (!is.null(facet_by)) {
    g <- g + theme(strip.text.x = element_text(size = text_sizes[6]))
    
  }
  
  finished.plot <- ggpubr::ggarrange(g, f,c, ncol = 1, align = "v",
                                     heights = c(0.5, 0.3, 0.1, 0.15))
  
  finished.plot <- ggpubr::ggarrange(finished.plot,legend, ncol = 1, align = "v",
                                     heights = c(0.9, 0.1))
  return(finished.plot)
}