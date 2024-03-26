### Authors/Moderators: Dr Felicity Lamrock, Miss Courtney Bull and Dr Sudhir Malla ----------------

dualGSEA <- function(data, group_data, group_colname, geneset_list) {
  
  ## Environment setup ----
  set.seed(127)
  
  ## Input data checkup ----
  
  ### check if the data is dataframe or matrix
  if(!all(sapply(data, is.numeric))) 
    stop("The input expression data needs to be all numeric. Genes as rows and samples as columns!",
         call. = F)
  
  ### check if the groups are n=2 and more or less
  ### NOTE: we will have this as dataframe with group column and rows as sammples (we require sample id to align with expression)
  if(!is.data.frame(group_data))
    stop("The group data input must be in a dataframe format with samples as rows.",
         call. = F)
  
  ### check if the samples are the same size in both
  if(ncol(data) !=  nrow(group_data))
    stop("The number of samples must be same size between expression and group data.",
         call. = F)
  
  ### extract and check group colname
  if(!group_colname %in% colnames(group_data))
    stop("The group name was not found in the group data.",
         call. = F)
  ### next, check if the group = 2. 
  if(is.character(group_data[ , group_colname %in% names(group_data)]))
    if(length(unique(group_data[ , group_colname %in% names(group_data)])) != 2)
      stop("There must be exact 2 groups.",
           call. = F)
  if(is.data.frame(group_data[ , group_colname %in% names(group_data)]))
    if(length(unique(group_data[ , group_colname %in% names(group_data)][[group_colname]])) != 2)
      stop("There must be exact 2 groups.",
           call. = F)
  
  
  ### next check if the sample ids matches to each other and aligns (order the group data first)
  group_data <- group_data %>% arrange(as.factor(.[[group_colname]]))
  ### first check if the sample is exists
  if(!all(colnames(data) %in% rownames(group_data)))
    stop("Please make sure the rownames in group data is your sample IDs")
  ### then check if they are ordered the same between expression and group data
  if(!all(colnames(data) == rownames(group_data)))
    data = data %>% dplyr::select(all_of(rownames(group_data)))
  
  
  ### checking geneset list
  #### check and remove duplicates
  geneset_list <- lapply(geneset_list, function(gene) {gene[!duplicated(gene)] }) 
  #### check and remove NA or blanks
  geneset_list <- lapply(geneset_list, function(gene) { gene[!is.na(gene) |
                                                               gene == ""] } )
  #### check number of matches with the genes in the data
  gene_match_number <- sum( unlist(lapply(geneset_list, function(gene) { 
    sum(gene %in% rownames(data))/length(gene)*100 
  } )) > 80)
  total_geneset <- length(geneset_list)
  gene_match_per <- (gene_match_number/total_geneset)*100
  message(sprintf("NOTICE: %1.0f%% (%s/%s) of geneset(s) have >80%% gene matches with the data.",
                  gene_match_per, gene_match_number, total_geneset))
  
  
  ## Differential Expression Analysis ----
  groupf <- as.factor(group_data[[group_colname]])
  group1 <- levels(groupf)[1]
  group2 <- levels(groupf)[2]
  
  design <- model.matrix(~ 0+groupf)
  ### Note: because the R will take number as its name, if its integer, we need to add letter infront
  if(is.integer(unique(group_data[[group_colname]]))) {
    colnames(design) <- c(paste0("G",group1), paste0("G",group2))
  } else {
    colnames(design) <- c(group1, group2)
  }
  
  fit <- limma::lmFit(data, design)
  contr <- list(paste0(colnames(design)[1],"-", colnames(design)[2]), levels=design)
  cont.matrix <- do.call(limma::makeContrasts, contr)
  
  fit2 <- limma::contrasts.fit(fit, cont.matrix)
  fit2 <- limma::eBayes(fit2)
  degs <- limma::topTable(fit2, adjust.method = "fdr", sort.by = "B", number = Inf)
  degs <- tibble::rownames_to_column(degs, var = "symbol")
  degs$combined_value <- degs$logFC * -log10(degs$P.Value)
  
  message("Differential Expression Analysis finished.")
  
  
  ## Pair-Wise Analysis (fGSEA) ----
  
  message("Running Pairwise GSEA! It may take a while...")
  
  gene_rank_list <- as.vector(degs$t)
  names(gene_rank_list) <- degs$symbol
  gene_rank_list <- sort(gene_rank_list, decreasing = TRUE)
  
  
  fgseaRes <- fgsea::fgsea(geneset_list, 
                    stats = gene_rank_list , 
                    nPermSimple = 1000,
                    minSize = 1, maxSize = Inf, nproc = 4, gseaParam = 1, BPPARAM = NULL)
  gsea.dt <- as.data.frame(fgseaRes)
  
  
  ## Bar plot fGSEA ----
  if(all(gsea.dt$NES > 0) == TRUE) {
    
    p0 <- ggplot(gsea.dt, aes(x = reorder(pathway, NES),y = NES,
                              fill = ifelse(padj < 0.05, "Highlighted", "Normal")))+
      geom_bar(stat = "identity") +
      coord_flip(clip = "off") + # to flip the bar
      theme(panel.background=element_blank(),
            panel.border=element_rect(fill = NA,colour="gray50"),
            axis.title.y=element_blank(), axis.title.x = element_text(size = 14),
            axis.text = element_text(size = 10),
            plot.title = element_text(hjust = 0.5, size = 15),
            legend.position= c(0.75,0.16),
            legend.key = element_rect(fill = "white"),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 12)) +
      scale_y_continuous(limits = c(0, max(gsea.dt$NES)+0.05)) +
      scale_x_discrete(limits=rev(gsea.dt$NAME)) +
      scale_fill_manual(name = "padj",
                        labels = c("< 0.05", "> 0.05"), 
                        values = c("red", "black"), guide = guide_legend(reverse = TRUE))+
      ylab("Normalized enrichment score (NES)") +
      ggtitle("GSEA Barplot") +
      annotation_custom(
      grob = grid::textGrob(
        label = group1,
        x = 0.8, y = 1.02, just = "center", gp = grid::gpar(fontsize = 10, face = "bold", col = "#1B80AD"))
    ) 
  } else if(all(gsea.dt$NES < 0) == TRUE) {
    
    p0 <- ggplot(gsea.dt, aes(x = reorder(pathway, NES),y = NES,
                              fill = ifelse(padj < 0.05, "Highlighted", "Normal")))+
      geom_bar(stat = "identity") +
      coord_flip(clip = "off") + # to flip the bar
      theme(panel.background=element_blank(),
            panel.border=element_rect(fill = NA,colour="gray50"),
            axis.title.y=element_blank(), axis.title.x = element_text(size = 14),
            axis.text = element_text(size = 10),
            plot.title = element_text(hjust = 0.5, size = 15),
            legend.position= c(0.75,0.16),
            legend.key = element_rect(fill = "white"),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 12)) +
      scale_y_continuous(limits = c(min(gsea.dt$NES)-0.05, 0)) +
      scale_x_discrete(limits=rev(gsea.dt$NAME)) +
      scale_fill_manual(name = "padj",
                        labels = c("< 0.05", "> 0.05"), 
                        values = c("red", "black"), guide = guide_legend(reverse = TRUE))+
      ylab("Normalized enrichment score (NES)") +
      ggtitle("GSEA Barplot") +
      annotation_custom(
        grob = grid::textGrob(
          label = group2,
          x = 0.2, y = 1.02, just = "center", gp = grid::gpar(fontsize = 10, face = "bold", col = "#EF7923"))
      ) 
    
  } else {
    
    p0 <- ggplot(gsea.dt, aes(x = reorder(pathway, NES),y = NES,
                              fill = ifelse(padj < 0.05, "Highlighted", "Normal")))+
      geom_bar(stat = "identity") +
      coord_flip(clip = "off") + # to flip the bar
      theme(panel.background=element_blank(),
            panel.border=element_rect(fill = NA,colour="gray50"),
            axis.title.y=element_blank(), axis.title.x = element_text(size = 14),
            axis.text = element_text(size = 10),
            plot.title = element_text(hjust = 0.5, size = 15),
            legend.position= c(0.75,0.16),
            legend.key = element_rect(fill = "white"),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 12)) +
      scale_y_continuous(limits = c(min(gsea.dt$NES)-0.05, max(gsea.dt$NES)+0.05)) +
      scale_x_discrete(limits=rev(gsea.dt$NAME)) +
      scale_fill_manual(name = "padj",
                        labels = c("< 0.05", "> 0.05"), 
                        values = c("red", "black"), guide = guide_legend(reverse = TRUE))+
      ylab("Normalized enrichment score (NES)") +
      ggtitle("GSEA Barplot")
    
    ### adding groups in the barplot
    p0 <- p0 + annotation_custom(
      grob = grid::textGrob(
        label = group1,
        x = 0.8, y = 1.02, just = "center", gp = grid::gpar(fontsize = 10, face = "bold", col = "#1B80AD"))
    ) 
    p0 <- p0 + annotation_custom(
      grob = grid::textGrob(
        label = group2,
        x = 0.2, y = 1.02, just = "center", gp = grid::gpar(fontsize = 10, face = "bold", col = "#EF7923"))
    ) 

  }
  
  #print(p0) ## Barplot
  
  
  ## Enrichment plots fGSEA ----
  
  ### Create pathways list
  pathways_l <- list()
  for (geneset_name in names(geneset_list)) {
    pathways_l[[geneset_name]] <- as.character(geneset_list[[geneset_name]])
  }
  
  ### Run fGSEA for each Hallmark
  enrichmentplot_list <- list()
  for (geneset_name in names(geneset_list)) {
    
    p1 <- fgsea::plotEnrichment(geneset_list[[geneset_name]], gene_rank_list, gseaParam = 1, ticksSize = 0.1) +
      labs(title=geneset_name, y = "Enrichment score") +
      theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), 
            panel.background = element_blank(),
            axis.title.y = element_text(size = 10, face = "bold"),
            axis.text.y = element_text(size = 8),
            panel.border = element_blank(), axis.title.x = element_blank(),
            axis.text.x = element_blank(), axis.ticks.x = element_blank())
    
    #### Fetch NES and padj values
    nes_value <- gsea.dt$NES[match(geneset_name, gsea.dt$pathway)]
    padj_value <- gsea.dt$padj[match(geneset_name, gsea.dt$pathway)]
    
    #### Add NES and padj values as separate grid text elements
    p1 <- p1 +
      annotation_custom(
        grob = grid::textGrob(
          label = ifelse(is.na(nes_value), "NES: NA", sprintf("NES: %.2f", nes_value)),
          x = 0.4, y = 0.98, just = "center", gp = grid::gpar(fontsize = 10)
        )) +
      annotation_custom(
        grob = grid::textGrob(
          label = sprintf("padj: %.2e", padj_value),
          x = 0.6, y = 0.98, just = "center", gp = grid::gpar(fontsize = 10)
        )) +
      annotation_custom(
        grob = grid::textGrob(
          label = group1,
          x = 0.2, y = 0.98, just = "center", gp = grid::gpar(fontsize = 10, face = "bold", col = "#1B80AD")
        )) +
      annotation_custom(
        grob = grid::textGrob(
          label = group2,
          x = 0.8, y = 0.98, just = "center", gp = grid::gpar(fontsize = 10, face = "bold", col = "#EF7923")
        ))
    
    print(p1) ## Enrichment plot
    
    enrichmentplot_list[[geneset_name]] <- recordPlot()
    invisible(dev.off())
    
  }
  
  message("Pairwise GSEA finished!")
  

  ## Single-Sample Analysis (SSGSEA) ----
  
  message("Running single sample GSEA (ssGSEA)...")
  
  ssgsea_scores_1 <- GSVA::gsva(
    as.matrix(data),
    geneset_list,
    method = "ssgsea",
    kcdf = "Gaussian",
    min.sz = 1,
    max.sz = Inf,
    mx.diff = TRUE,
    verbose = TRUE
  )
  
  message("ssGSEA finished!")
  
  ssgsea_data <- as.data.frame(ssgsea_scores_1)
  ssgsea_data <- ssgsea_data %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "ID")
  
  GR <- data.frame(ID = row.names(group_data), Group=group_data[[group_colname]] )
  SSGSEA <- merge(ssgsea_data, GR, by = "ID") %>% dplyr::select("ID", "Group", everything())
  SSGSEA$Group <- factor(SSGSEA$Group, levels = c(levels(groupf)))
  
  
  ## SSGSEA Plots ----
  
  ### Density Plots ----
  
  message("Generating ssGSEA: Density Plots!")
  
  densityplot_list <- list()
  for (geneset_name in names(geneset_list)) {
    
    #### plot
    p2 <- ggplot(SSGSEA, aes_string(x = geneset_name, y = "Group", fill = "Group")) +
      ggridges::geom_density_ridges(rel_min_height = 0.01,
                          jittered_points = TRUE,
                          position = ggridges::position_points_jitter(width = 0.05, height = 0), 
                          point_shape = '|', point_size = 4, point_alpha = 1, alpha = 0.7) +
      scale_color_manual(values = c("#1B80AD", "#EF7923")) +
      scale_fill_manual(values = c("#1B80AD", "#EF7923")) +
      ggridges::theme_ridges() +
      theme(legend.position="bottom",
            axis.line = element_line(),
            plot.title = element_text(hjust = 0.5, size=16, face="bold"),
            axis.title.y = element_text(size = 14,  colour="black", face="bold", hjust = 0.5),
            axis.text.y = element_text(size = 12,  colour="black", face="bold"),
            axis.title.x = element_text(size = 14,  colour="black", face="bold", hjust = 0.5) ,
            axis.text.x = element_text(size = 12,  colour="black", face="bold"),
            panel.background=element_blank(),
            panel.grid = element_blank(),
            panel.spacing = unit(0.1, "lines"),
            strip.background = element_blank()) +
      xlab(paste(geneset_name, "ssGSEA score")) +
      ylab("") +
      ggtitle(geneset_name)
    
    suppressMessages(print(p2)) ## Density plot
    densityplot_list[[geneset_name]] <- recordPlot()
    invisible(dev.off())
    
  }
  
  ### Histograms ----
  
  message("Generating ssGSEA: Histograms!")
  
  histogram_list <- list()
  for (geneset_name in names(geneset_list)) {
    p3 <- easyGgplot2::ggplot2.histogram(data=SSGSEA, xName=geneset_name, groupName='Group',
                            legendPosition="top",
                            groupColors=c("#1B80AD","#EF7923"), alpha=0.5,
                            addDensity=TRUE,
                            removePanelGrid=TRUE,removePanelBorder=TRUE,
                            backgroundColor="white",
                            axisLine=c(0.5, "solid", "black")) +
      theme(legend.position="bottom",
            axis.line = element_line(),
            plot.title = element_text(hjust = 0.5, size=16, face="bold"),
            axis.title.y = element_text(size = 14,  colour="black", face="bold"),
            axis.text.y = element_text(size = 12,  colour="black", face="bold"),
            axis.title.x = element_text(size = 14,  colour="black", face="bold") ,
            axis.text.x = element_text(size = 12,  colour="black", face="bold"),
            panel.background=element_blank(),
            panel.spacing = unit(0.1, "lines"),
            strip.background = element_blank()) +
      xlab(paste(geneset_name, "ssGSEA score")) +
      ylab("Assigned Probability (%)")+
      ggtitle(geneset_name)
    
    suppressMessages(print(p3)) ## histogram
    histogram_list[[geneset_name]] <- recordPlot()
    invisible(dev.off())
  }
  
  
  
  ### ROC plots ----
  
  message("Generating ssGSEA: ROC!")
  
  rocplot_list <- list()
  for (geneset_name in names(geneset_list)) {
    ssgsea_scores <- SSGSEA[[geneset_name]]
    mean_group1 <- mean(ssgsea_scores[SSGSEA$Group == group1])
    mean_group2 <- mean(ssgsea_scores[SSGSEA$Group == group2])
    
    ### If the mean is higher for group 2, we will leave the label as it is, so that it will 
    ### tell us whether the signature is truly predictive for calling group 2
    if (mean_group2 >= mean_group1){
      pred <- ROCR::prediction(SSGSEA[[geneset_name]], SSGSEA$Group)
      perf <- ROCR::performance(pred, "tpr", "fpr")
      auc <- ROCR::performance(pred, "auc")
      ROCR::plot(perf, colorize = TRUE,
           main = paste0(geneset_name, "\nPrediction for ", sprintf("%s", group2), ": AUC = ", round(auc@y.values[[1]], 2))
           )
      }
    else{
      ### The switch of group labels is made in this instance so that signature is predicting FOR the 
      ### group-of-interest (positive), and positively is the signature selecting for that group!
      SSGSEA$Group_switched <- SSGSEA$Group
      SSGSEA$Group_switched <- case_when(SSGSEA$Group_switched == group1 ~ group2,
                                         SSGSEA$Group_switched == group2 ~ group1)
      
      pred <- ROCR::prediction(SSGSEA[[geneset_name]], SSGSEA$Group_switched)
      perf <- ROCR::performance(pred, "tpr", "fpr")
      auc <- ROCR::performance(pred, "auc")
      ROCR::plot(perf, colorize = TRUE, 
           main = paste0(geneset_name, "\nPrediction for ", sprintf("%s", group1), ": AUC = ", round(auc@y.values[[1]], 2))
           )
    }
    
    rocplot_list[[geneset_name]] <- recordPlot()
    invisible(dev.off())
    
  }
  
  ### Waterfall ----
  
  message("Generating ssGSEA: Waterfall Plots! It may take a while...")
  
  #### empty plot list
  waterfallplot_list <- list()
  #### empty dataframe - this is to save the labels
  data_labels <- matrix(nrow = nrow(SSGSEA))
  rownames(data_labels) <- SSGSEA$ID
  data_labels <- as.data.frame(data_labels)
  
  
  for (geneset_name in names(geneset_list)) {
    
    #### Find optimum cutpoint
    set.seed(127)
    opt_cut <- suppressMessages(cutpointr::cutpointr(x = SSGSEA[[geneset_name]], 
                                          class = SSGSEA[["Group"]], 
                                          pos_class = group1, direction = ">=", boot_runs = 1000
                                          , allowParallel = T
    ))
    cutpoint <- opt_cut$optimal_cutpoint
    
    
    if(is.infinite(cutpoint)) {
      message(sprintf("The optimal cutpoint for %s was found Inf. The waterfallplot will not be recorded.",
                      geneset_name))
      next
    }
    
    #### make new column for above and below cutpoint
    SSGSEA[["cutoff_labels"]] <- ""
    #### create label
    SSGSEA$cutoff_labels[SSGSEA[[geneset_name]] >= cutpoint ] <- "Above"
    SSGSEA$cutoff_labels[SSGSEA[[geneset_name]] < cutpoint ] <- "Below"
    
    
    #### add new column for pred_test
    SSGSEA[["roc_centered"]] <- SSGSEA[[geneset_name]] - cutpoint 
    
    
    #### order the roc-centered column
    SSGSEA <- SSGSEA[order(SSGSEA$roc_centered, decreasing = T),]
    #### find midpoint for the plot
    midpoint_line <- SSGSEA[SSGSEA$roc_centered == 0 ,]$ID
    
    
    #### plot
    p4 <- ggplot(SSGSEA, 
                 aes(x=reorder(ID, roc_centered), y=roc_centered, fill=Group, color=Group)) +
      facet_wrap(~Group) +
      geom_bar(stat="identity", width=0.5, position = position_dodge(width=0.4)) +
      geom_vline(xintercept = midpoint_line, color = "red", linetype = "dotted") +
      scale_fill_manual(values = c("#1B80AD", "#EF7923")) +
      scale_colour_manual(values = c("#1B80AD", "#EF7923")) +
      labs(x = "Samples", y = geneset_name, title = geneset_name) +
      theme(axis.line = element_line(),
            legend.position = "none",
            legend.text = element_text(size=10,colour="black", face="bold"),
            legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5, size=10, face="bold"),
            axis.title.y = element_text(size = 8,  colour="black", face="bold"),
            axis.text.y = element_blank(),
            axis.title.x = element_text(size = 8,  colour="black", face="bold") ,
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            panel.background=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border = element_blank(),
            strip.background = element_blank())
    
    print(p4) ## waterfallplot
    waterfallplot_list[[geneset_name]] <- recordPlot()
    invisible(dev.off())
    
    ### into the empty dataframe
    data_labels[[geneset_name]] <- SSGSEA$cutoff_labels
    
  }
  
  ## remove the first column of the data label
  data_labels[,1] <- NULL
  
  
  message("GSEA has finished running!")
  
  
  ## Outputs ----
  return(list("Differential_GeneRanking" = degs, 
              "Pairwise_ResultTable" = gsea.dt,
              "Pairwise_BarPlot" = p0,
              "Pairwise_EnrichmentPlots" = enrichmentplot_list,
              "SingleSample_ResultTable" = ssgsea_data,
              "SingleSample_DensityPlots" = densityplot_list,
              "SingleSample_Histograms" = histogram_list,
              "SingleSample_ROCPlots" = rocplot_list,
              "SingleSample_OptimumCutoff_Labels" = data_labels,
              "SingleSample_WaterfallPlots" = waterfallplot_list))
  
}
