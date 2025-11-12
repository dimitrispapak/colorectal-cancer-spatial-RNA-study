library(umap)
library(parallel)
library(dplyr)
library(ggrepel) 
library(fpc)
library(pheatmap)

run_umap <- function(dataframe,three_d=FALSE){
	custom_umap <- umap::umap.defaults
	custom_umap$random_state <- 42
	if (three_d == FALSE){
		umap_out <- umap(t(log2(assayDataElement(dataframe , elt = "q_norm"))),  config = custom_umap,n_neighbors=15)
		pData(dataframe)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]
		dataframe
	}else{

		custom_umap$n_components <- 3
		umap_out <- umap(t(log2(assayDataElement(dataframe , elt = "q_norm"))),  config = custom_umap,n_neighbors=15)
		pData(dataframe)[, c("UMAP1", "UMAP2","UMAP3")] <- umap_out$layout[, c(1,2,3)]
		dataframe
	}
}
umap_plot <- function(dataframe,title,label, colors = FALSE,label2=FALSE,show_score=FALSE){
	if (show_score == TRUE){
		metrics <- cluster_metrics(dataframe,label)
		score = paste("Score: ",round(metrics[['silhouette']],2))    
	}else{
		score = ""
	}    
	if (label2 == FALSE){
		plot <- ggplot(pData(dataframe), aes(x = UMAP1, y = UMAP2, color = pData(dataframe)[[label]]))
	} else {
		plot <- ggplot(pData(dataframe), aes(x = UMAP1, y = UMAP2, color = pData(dataframe)[[label]], shape  = pData(dataframe)[[label2]]))
	}
	if (colors != FALSE){
		plot <- plot + scale_color_manual(values = colors)
	} 
	plot <- plot +
	geom_point(size = 3) +
	ggtitle(title) +
	coord_cartesian(clip = 'off')  +
	annotate("text", x = Inf, y = Inf, hjust = 0, vjust = 1, label = score) +
	theme(plot.title = element_text(size=22),
	      legend.title=element_blank(),
	      panel.grid.major = element_blank(),
	      panel.grid.minor = element_blank(),
	      panel.border = element_blank(),
	      panel.background = element_blank(),
	      axis.text.x=element_blank(),
	      axis.text.y=element_blank(),
	      axis.ticks=element_blank(),
	      axis.title.x=element_blank(),
	      axis.title.y=element_blank(),
	      legend.text=element_text(size=20)
	      ) +
	scale_shape_manual(values=c(15,5,16,17))
	return(plot)
}
umap_plot2 <- function(dataset,title,cat1,cat2,cat3,shape){
	colors <- c('#00BEC4','#FFFFFF','#F8766C')
	names(colors) <- c('primitive','no_color','metastasis')
	df <- pData(dataset)
	df$new <- ifelse(df[,cat2] == cat3,df[,cat1],'no_color')
	plot <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = df$new)) +
	geom_point(size = 3,shape=shape) +
	ggtitle(title) +
	coord_cartesian(clip = 'off')  +
	theme(plot.title = element_text(size=22),
	      legend.title = element_blank(),
	      panel.grid.major = element_blank(),
	      panel.grid.minor = element_blank(),
	      panel.border = element_blank(),
	      panel.background = element_blank(),
	      axis.text.x = element_blank(),
	      axis.text.y = element_blank(),
	      axis.ticks = element_blank(),
	      axis.title.x = element_blank(),
	      axis.title.y = element_blank(),
	      legend.text = element_text(size=20)) + 
		scale_color_manual(values = colors,breaks=c('metastasis','primitive'))
	return(plot)
}
cluster_metrics <- function(S4_obj,column){
	data <-pData(S4_obj)
	#print(head(data))
	# Calculate the Euclidean distance matrix
	dist_matrix <- dist(data[, c("UMAP1", "UMAP2")])
	# Obtain the clustering vector from the 'label' column
	labels <- as.factor(data[,column])
	clustering_vector <- as.integer(labels)
	# Calculate clustering statistics
	clustering_stats <- cluster.stats(dist_matrix, clustering_vector)
	# Silhouette score
	silhouette_score <- mean(clustering_stats[['avg.silwidth']])
	# Calinski-Harabasz index
	calinski_harabasz_index <- clustering_stats[['ch']]
	result <- list("silhouette" = silhouette_score, "harabasz" = calinski_harabasz_index)
	return(result)
}


convertMouseGeneList <- function(x){
	require("biomaRt")
	human = useMart(biomart="ensembl", dataset = "hsapiens_gene_ensembl", verbose = TRUE, host = "dec2021.archive.ensembl.org")
	mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", verbose = TRUE, host = "dec2021.archive.ensembl.org")
	genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
	humanx <- unique(genesV2[, 2])
	# Print the first 6 genes found to the screen

	print(head(humanx))
	return(humanx)
}
volcano_plot <- function(df,contrast,title=NULL,goi1 = NULL,goi2 = NULL){
	goi1_name <- deparse(substitute(goi1))
	goi2_name <- deparse(substitute(goi2))
	top_g <- c()
	df <- subset(df, Contrast == contrast)
	top_g <- c(top_g,df[order(df$invert_P),]$Gene[1:15],df[order(-df$invert_P),]$Gene[1:15])
	top_g <- unique(top_g)
	df <- df[, -1*ncol(df)] # remove invert_P from matrix
	splitted = strsplit(contrast,split=" - ")

	left = splitted[[1]][2]
	right = splitted[[1]][1]

	x =  paste("Enriched in ",left," <- log2(FC) -> ", right,sep="")
	# Graph results
	p <- ggplot(df,
		    aes(x = Estimate, y = -log10(`Pr(>|t|)`),
			color = Color, label = Gene)) +
geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
geom_hline(yintercept = -log10(0.05), lty = "dashed") +
geom_point()

if (is.null(goi1)){
	p <- p + geom_text_repel(data = subset(df, Gene %in% top_g & FDR < 0.05 & (Estimate>0.5| Estimate < -0.5)),
				 size = 3, point.padding = 0.15, color = "black",
				 min.segment.length = .1, box.padding = .2, lwd = 2,
				 max.overlaps = 50) 
} else {
	p <- p + geom_text_repel(data = subset(df, Gene %in% goi1 ),
				 size = 3, point.padding = 0.15, color = "black",
				 min.segment.length = .1, box.padding = .2, lwd = 2,
				 max.overlaps = 50)
	if (!is.null(goi2)){
		p <- p + geom_text_repel(data = subset(df, Gene %in% goi2 ),
					 size = 3, point.padding = 0.15, color = "black",
					 min.segment.length = .1, box.padding = .2, lwd = 2,
					 max.overlaps = 50) 
	}

}

p <- p +
labs(x = x,
     y = "Significance, -log10(P)",
     color = "Significance" ,
     title = title) +
scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue",
			      `FDR < 0.05` = "lightblue",
			      `P < 0.05` = "orange2",
			      `NS or FC < 0.5` = "gray"),
		   guide = guide_legend(override.aes = list(size = 3))) +
scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
theme_bw(base_size = 13) +
theme(legend.position = "bottom") 

print(p)
}

# Random slope is included in the LMM model when comparing features that co-exist in a given tissue section 
#	(e.g. glomeruli vs tubules in DKD kidneys),
# LMM model does not require a random slope when comparing features that are mutually exclusive in a given tissue section 
#	(healthy glomeruli versus DKD glomeruli)
lmm_diff_expression <- function(df,study_column,between=TRUE, slide = TRUE){
	results <- c()
	uniq_values = unique(pData(df)[,study_column])
	print(uniq_values)
	print(paste("categories of column:",study_column,uniq_values))
	if (length(uniq_values) < 2){
		return(FALSE)
	}
	pData(df)$testRegion <- factor(pData(df)[[study_column]],uniq_values )
	pData(df)[["slide"]] <- factor(pData(df)[["slide name"]])
	pData(df)[["patient"]] <- factor(pData(df)[["patient_id"]])
	if (between == TRUE){
		if (slide == TRUE){
			formula = ~ testRegion + (1 | slide)
		} else {
			formula = ~ testRegion + (1 | patient)
		}

	} else {
		if (slide == TRUE){
			formula = ~ testRegion + (1 + testRegion | slide)
		} else {
			formula = ~ testRegion + (1 + testRegion | patient)
		}
	}
	assayDataElement(object = df, elt = "log_q") <- assayDataApply(df, 2, FUN = log, base = 2, elt = "q_norm")

	print(dim(assayDataElement(df, "log_q") ))
	mixedOutmc <- mixedModelDE(df,
				   elt = "log_q",
				   modelFormula = formula,
				   groupVar = "testRegion",
				   nCores = 10,
				   multiCore = F)

	# format results as data.frame
	r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
	tests <- rownames(r_test)
	r_test <- as.data.frame(r_test)
	r_test$Contrast <- tests

	# use lapply in case you have multiple levels of your test factor to
	# correctly associate gene name with it's row in the results table
	r_test$Gene <- unlist(lapply(colnames(mixedOutmc), rep, nrow(mixedOutmc["lsmeans", ][[1]])))
	r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
	r_test <- r_test[, c("Gene", "Contrast", "Estimate", 
			     "Pr(>|t|)", "FDR")]
		results <- rbind(results, r_test)

		kable(subset(head(results,20)), digits = 3,
		      caption = "DE results for Genes of Interest",
		      align = "lc", row.names = FALSE)


		results$Color <- "NS or FC < 0.5"
		results$Color[results$`Pr(>|t|)` < 0.05] <- "P < 0.05"
		results$Color[results$FDR < 0.05] <- "FDR < 0.05"
		results$Color[results$FDR < 0.001] <- "FDR < 0.001"
		results$Color[abs(results$Estimate) < 0.5] <- "NS or FC < 0.5"
		results$Color <- factor(results$Color,
					levels = c("NS or FC < 0.5", "P < 0.05",
						   "FDR < 0.05", "FDR < 0.001"))
		results$invert_P <- (-log10(results$`Pr(>|t|)`)) * sign(results$Estimate)
		return(results)
}

plot_ontology <- function(genes,title){
	# Convert gene symbols to Entrez IDs
	entrez_ids <- mapIds(org.Hs.eg.db, 
			     keys = genes, 
			     keytype = "SYMBOL", 
			     column = "ENTREZID")
	# Perform the GO enrichment analysis
	go_enrichment <- enrichGO(gene         = entrez_ids,
				  OrgDb        = org.Hs.eg.db,  # Human
				  ont          = "ALL",  # Biological process ontology
				  pAdjustMethod = "BH",  # Benjamini & Hochberg adjustment method
				  qvalueCutoff = 0.05)  # q-value cutoff

	# Print the results
	results_length <- nrow(as.data.frame(go_enrichment))
	if (results_length > 0){
		p <- plot(barplot(go_enrichment, showCategory = 10,label_format=50, title=title)) + theme(axis.text.y=element_text(face="bold", size=8))
	} else {
		return(NULL)
	}

}


heatmap <- function(df,results,title) {

	GOI = head(subset(results[order(results$`Pr(>|t|)`, decreasing = FALSE),], `Pr(>|t|)` < 0.05)$Gene, n = 100)
	#GOI <- unique(subset(results, `Pr(>|t|)` < 0.05)$Gene)
	print(GOI)
	p <- pheatmap(log2(assayDataElement(df[GOI, ], elt = "q_norm")),
		      main = title,
		      scale = "row", 
		      fontsize_row = 4,
		      show_rownames = TRUE, show_colnames = FALSE,
		      border_color = NA,
		      clustering_method = "average",
		      clustering_distance_rows = "correlation",
		      clustering_distance_cols = "correlation",
		      treeheight_row = 0, treeheight_col = 0,
		      cutree_cols = 2, cutree_rows = 2,
		      breaks = seq(-3, 3, 0.05),
		      color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
		      annotation_col = pData(df)["primitive_metastasis"],
		      annotation_names_col=FALSE)
	p
}

positions_of_change  <- function(vec){
	res <- c()
	for (i in 1:length(vec)){
		if (i+1 <= length(vec)){
			if (vec[i] != vec[i + 1]){
				res <- c(res,i)
			}
		}
	}
	return(res)
}

deg_heatmap_prim_meta <- function(data,results_df,features,title) {
	GOI = head(subset(results[order(results_df$`Pr(>|t|)`, decreasing = FALSE),], `Pr(>|t|)` < 0.05 & Estimate > 0.5 | Estimate < -0.5)$Gene, n = 100)
	feature <- features[1]
	print(feature)
	meta <- pData(data)[order(pData(data)[feature]), ]
	print(colnames(meta))
	data <- as.data.frame(log2(assayDataElement(data[GOI,], elt = "q_norm")))
	data <- select(data,all_of(rownames(meta)))
	indices <- positions_of_change(meta[[feature]])
	print(indices)
	print(class(indices))
	#GOI <- unique(subset(results, `Pr(>|t|)` < 0.05)$Gene)
	p <- pheatmap(data,
		      main = title,
		      scale = "row", 
		      fontsize =8,
		      fontsize_row = 3,
		      show_rownames = TRUE, show_colnames = FALSE,
		      border_color = NA,
		      cluster_cols = F,
		      #clustering_method = "average",
		      #clustering_distance_rows = "correlation",
		      #clustering_distance_cols = "correlation",
		      treeheight_row = 0, treeheight_col = 0,
		      #cutree_cols = 2, cutree_rows = 2,
		      #breaks = seq(-3, 3, 0.05),
		      gaps_col = c(indices),
		      color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
		      annotation_col = meta[features],
		      annotation_names_col=FALSE)
	return(p)

}

extreme_genes <- function(results,n){
	positive <- subset(results[order(results$`Pr(>|t|)`, decreasing = FALSE),], `Pr(>|t|)` < 0.05 & Estimate > 0.5 )
	i <- 1
	res <- NULL
	while (i >= 0.5) {
		i <- i - 0.01
		x_quartiles <- quantile(positive$Estimate, probs = c(i))
		y_quartiles <- quantile(positive$invert_P, probs = c(i))
		extrem_pos <-	subset(positive[,c('Gene','Estimate','invert_P')], Estimate > x_quartiles & invert_P > y_quartiles)$Gene # Both Log Fold change and p-value
		#print(length(extrem_pos))
		if (n < length(extrem_pos)){
			res <- extrem_pos
			break
		}
	}
	negative <- subset(results[order(results$`Pr(>|t|)`, decreasing = FALSE),], `Pr(>|t|)` < 0.05 & Estimate < -0.5 )
	negative$Estimate <- negative$Estimate * (-1)
	negative$invert_P <- negative$invert_P * (-1) 
	i <- 1
	res <- NULL
	while (i >= 0.5) {
		i <- i - 0.01
		x_quartiles <- quantile(negative$Estimate, probs = c(i))
		y_quartiles <- quantile(negative$invert_P, probs = c(i))
		extrem_neg <-	subset(negative[,c('Gene','Estimate','invert_P')], Estimate > x_quartiles & invert_P > y_quartiles)$Gene
		#print(length(extrem_neg))
		if (n < length(extrem_neg)){
			res <- extrem_neg
			break
		}
	}
	return(list(pos= extrem_pos,neg=extrem_neg))
}

expressions2 <- function(results,df){
	# Positive
	positive <- subset(results[order(results$`Pr(>|t|)`, decreasing = FALSE),], `Pr(>|t|)` < 0.05 & Estimate > 0.5 )
	positive$Expression <- rowMeans(assayDataElement(df[positive$Gene, ], elt = "q_norm")) # Add Expression
	positive$Log_Expression <- log10(positive$Expression) # log expression
	positive$invert_fdr <- -log10(positive$FDR) # log False Discovery rate
	positive$invert_fdr_scaled <- scale(positive$invert_fdr) # Scale
	positive$Estimate_scaled <- scale(positive$Estimate) # Scale
	positive$Log_Expression_scaled <- scale(positive$Log_Expression) # Scale
	positive$metric <- sign(positive$invert_fdr_scaled)*(positive$invert_fdr_scaled ^ 2) +sign(positive$Log_Expression_scaled)* (positive$Log_Expression_scaled ^ 2) + sign(positive$Estimate_scaled)*(positive$Estimate_scaled ^ 2) # Define comparison metric
	positive <- positive[order(positive$metric,decreasing = TRUE),] # order by metric
	# Negative
	negative <- subset(results[order(results$`Pr(>|t|)`, decreasing = FALSE),], `Pr(>|t|)` < 0.05 & Estimate < -0.5 )
	negative$Expression <- rowMeans(assayDataElement(df[negative$Gene, ], elt = "q_norm"))
	negative$Log_Expression <- log10(negative$Expression)
	negative$Estimate <- negative$Estimate * (-1)
	negative$invert_fdr <- -log10(negative$FDR) 
	negative$invert_fdr_scaled <- scale(negative$invert_fdr)
	negative$Estimate_scaled <- scale(negative$Estimate)
	negative$Log_Expression_scaled <- scale(negative$Log_Expression)
	negative$metric <- sign(negative$invert_fdr_scaled)*(negative$invert_fdr_scaled ^ 2) +sign(negative$Log_Expression_scaled)* (negative$Log_Expression_scaled ^ 2) + sign(negative$Estimate_scaled)*(negative$Estimate_scaled ^ 2)
	negative <- negative[order(negative$metric,decreasing = TRUE),]
	return(list(pos= positive,neg=negative))
}

enrich <- function(df1,df2,segment){
	df2_column <- paste0(tolower(segment),'_nuclei')
	matches <- match(paste(df2$slide, df2$roi, segment),paste(df1$`slide name`, as.integer(df1$roi),df1$segment))
	for (i in seq_along(matches)) {
		if (!is.na(matches[i])) {
			df1[matches[i],'nuclei'] <- df2[i,df2_column]
			df1[matches[i],'ck_score'] <- df2[i,'ck_score']
			df1[matches[i],'cd45_score'] <- df2[i,'cd45_score']
		}
	}
	return(df1)
}

extreme_genes_edger <- function(results,n){
	results$invert_p <- (-log10(results$FDR)) * sign(results$logFC)
	results$gene <- rownames(results)
	positive <- subset(results[order(results$FDR, decreasing = FALSE),], FDR < 0.05 & logFC > 0.5 )
	pos <- head(positive$gene,n)
	negative <- subset(results[order(results$FDR, decreasing = FALSE),], FDR < 0.05 & logFC < -0.5 )
	neg <- head(negative$gene,n)
	return(list(pos= pos,neg=neg))
}

volcano_edger <- function(df,left,right,title,goi1,goi2,goi3 = NULL){
	df$color <- "ns or fc < 0.5"
	df$color[df$`PValue` < 0.05] <- "p < 0.05"
	df$color[df$FDR < 0.05] <- "FDR < 0.05"
	df$color[df$FDR < 0.001] <- "FDR < 0.001"
	df$color[abs(df$logFC) < 0.5] <- "ns or fc < 0.5"
	df$color <- factor(df$color,levels = c("ns or fc < 0.5", "p < 0.05","FDR < 0.05", "FDR < 0.001"))
	df$gene <- rownames(df)
	df$invert_p <- (-log10(df$`PValue`)) * sign(df$logFC)
	top_g <- c()
	top_g <- c(top_g,df[order(df$invert_p),]$gene[1:15],df[order(-df$invert_p),]$gene[1:15])
	top_g <- unique(top_g)

	p <- ggplot(df,aes(x = logFC, y = -log10(`PValue`),color = color, label = gene)) +
	geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
	geom_hline(yintercept = -log10(0.05), lty = "dashed") +
	geom_point()
	x =  paste("enriched in ",left," <- log2(fc) -> ", right,sep="")


	if (!is.null(goi1)){

		p <- p + geom_text_repel(data = subset(df, gene %in% goi1 ),
					 size = 3, point.padding = 0.15, color = "black",
					 min.segment.length = .1, box.padding = .2, lwd = 2,
					 max.overlaps = 50)
	}
	if (!is.null(goi2)){
		p <- p + geom_text_repel(data = subset(df, gene %in% goi2 ),
					 size = 3, point.padding = 0.15, color = "black",
					 min.segment.length = .1, box.padding = .2, lwd = 2,
					 max.overlaps = 50) 
	}
	if (!is.null(goi3)){
		p <- p + geom_text_repel(data = subset(df, gene %in% goi3 ),
					 size = 3, point.padding = 0.15, color = "red",
					 min.segment.length = .1, box.padding = .2, lwd = 2,
					 max.overlaps = 50) 
	}

	p <- p +
	labs(x = x,
	     y = "significance, -log10(p)",
	     color = "significance" ,
	     title = title) +
scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue",
			      `FDR < 0.05` = "lightblue",
			      `p < 0.05` = "orange2",
			      `ns or fc < 0.5` = "gray"),
		   guide = guide_legend(override.aes = list(size = 3))) +
scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
theme_bw(base_size = 13) +
theme(legend.position = "bottom") 
return(p)
}
perform_pathway_analysis <- function(gene_list, p_cutoff = 0.05, q_cutoff = 0.2, num_categories = 10, plot_title = "pathway enrichment analysis") {
	# convert gene symbols to entrez ids
	entrez_ids <- tryCatch({
		bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
	}, error = function(e) {
		message("error converting gene symbols to entrez ids: ", e$message)
		return(NULL)
	})

	if (is.null(entrez_ids) || nrow(entrez_ids) == 0) {
		stop("no valid entrez ids found. please check your gene list.")
	}

	# perform reactome pathway enrichment analysis
	reactome_enrich <- tryCatch({
		ReactomePA::enrichPathway(gene = entrez_ids$ENTREZID,
			      organism = "human",
			      pvalueCutoff = p_cutoff,
			      qvalueCutoff = q_cutoff)
	}, error = function(e) {
		message("Error in Reactome enrichment: ", e$message)
		return(NULL)
	})

	# Perform GO Biological Process enrichment analysis
	go_enrich <- tryCatch({
		enrichGO(gene = entrez_ids$ENTREZID,
			 OrgDb = org.Hs.eg.db,
			 ont = "BP",
			 pvalueCutoff = p_cutoff,
			 qvalueCutoff = q_cutoff)
	}, error = function(e) {
		message("Error in GO enrichment: ", e$message)
		return(NULL)
	})

	# Plot results if available
	if (!is.null(reactome_enrich) && nrow(reactome_enrich) > 0) {
		reactome_plot <- dotplot(reactome_enrich, showCategory = min(num_categories, nrow(reactome_enrich))) +
		ggtitle(paste("Reactome:", plot_title)) +
		theme(plot.title = element_text(hjust = 0.5))
		print(reactome_plot)
	} else {
		reactome_plot <- NULL
		message("No Reactome enrichment results to plot.")
	}

	if (!is.null(go_enrich) && nrow(go_enrich) > 0) {
		go_plot <- dotplot(go_enrich, showCategory = min(num_categories, nrow(go_enrich))) +
		ggtitle(paste("GO Biological Process:", plot_title)) +
		theme(plot.title = element_text(hjust = 0.5))
		print(go_plot)
	} else {
		go_plot <- NULL
		message("No GO enrichment results to plot.")
	}

	# Return results
#	return(list(reactome = reactome_plot, go = go_plot))
}
