#setwd("~")
#library(tictoc)
#library(Seurat)
library(tidyverse)
library(Matrix)
library(plotly)
library(furrr)
library(jsonlite)
# library(rtracklayer)
library(RVenn)

extract_gene_locs = function(features_list, scrna_df){
  
  feature_loc_list = rownames(scrna_df)[rownames(scrna_df) %in% features_list$GeneID]
  
  extracted_gene_locs = features_list[which(features_list$GeneID %in% feature_loc_list),] %>% 
  select(GeneID, Name) %>%
  unique()

  dup_check = duplicated(extracted_gene_locs$GeneID) %>% 
              sum()
  
  if(dup_check >= 1){

    extracted_gene_locs = extracted_gene_locs[-which(duplicated(extracted_gene_locs$GeneID)),]
    
    return(extracted_gene_locs)
    
    } else {
      
      return(extracted_gene_locs)
    
    }
}
intramutate = function(df, orig_search_list, replacement_list){
  
  for (i in 1:length(replacement_list)){
    
    df = case_when(
      df == orig_search_list[i] ~ replacement_list[i],
      df != orig_search_list[i] ~ df
    )
  }
  
  return(df)
}
dev_stage_cell_pop_stats = function(dev_stage_scrna){

  dev_stage_available_idents = levels(Idents(dev_stage_scrna))

  dev_stage_cell_ident_num_cell = map(.x = dev_stage_available_idents, 
                                      .f = function(x)(length(Cells(subset(dev_stage_scrna, idents = x))))) %>% 
                                  unlist()

  dev_stage_cell_ident_percentage = map(.x = dev_stage_available_idents,
                                        .f = function(x)(length(Cells(subset(dev_stage_scrna, idents = x)))/length(Cells(dev_stage_scrna)))) %>% 
                                    unlist()

  stats_df = tibble(
    dev_stage_available_idents, 
    dev_stage_cell_ident_num_cell, 
    dev_stage_cell_ident_percentage
          )

  return(stats_df)
}
lab_FP = function(scrna_df, feature_df){
  
  plot_list = FeaturePlot(scrna_df, features = feature_df$GeneID)
  
  for (i in 1:nrow(feature_df)){

    plot_list[[i]] = plot_list[[i]] + ggtitle(feature_df$Name[i])

  }

  return(plot_list)

}
lab_DP = function(scrna_df, feature_df, plot_title, by_stage = F, plot_height, interactive, col.min, col.max, dot.min, cols){
  
  
  dp = DotPlot(scrna_df, features = feature_df$GeneID, col.min = col.min, col.max = col.max, dot.min = dot.min, cols = cols) + 
    scale_x_discrete(breaks=c(feature_df$GeneID),labels=c(feature_df$Name)) +
    RotatedAxis() +
    ggtitle(as.character(plot_title)) + 
    coord_flip()
  
    if (interactive == T){
      
      dp = ggplotly(dp, height = plot_height)
      
    }
  

    if (by_stage == T){
      
      dp = DotPlot(scrna_df, features = feature_df$GeneID, group.by = "orig.ident") + 
        scale_x_discrete(breaks=c(feature_df$GeneID),
                         labels=c(feature_df$Name)) +
        RotatedAxis() +
        ggtitle(as.character(plot_title)) + 
        coord_flip()
      
    }
  
  # htmlwidgets::saveWidget(as_widget(dp), paste0(name,"_index.html"))
  print(dp)
  dev.off()

  return(dp)

}
lab_multi_DP = function(scrna_df_list, title_list, lab_DP_feature_df, pdf_title, interactive = F){
  
  DP_list = map2(

    .x = scrna_df_list, 
    .y = title_list, 
    .f = function(x,y) (lab_DP(x, lab_DP_feature_df, y)

      )
    )
  
  if(interactive == F){

    pdf(paste0(pdf_title, ".pdf"), onefile=T, width=12, height=8)
    
    map(DP_list, print) %>% 
    invisible()
    
    dev.off()
    return(DP_list)  
  
  } else {

    dir.create(pdf_title)
    
    DP_list_interactive = map(DP_list, ggplotly) %>% invisible()
    
    map2(
      .x=DP_list_interactive, 
      .y = title_list, 
      .f = function(x,y)(htmlwidgets::saveWidget(as_widget(x), paste0(pdf_title,"/",y,".html"),selfcontained=TRUE))) %>% 
    invisible()
  }
}
find_markers_export = function(devolist, file_names_list){

  map2(.x = devolist,
     .y = file_names_list,
     .f = function(x,y)
      (
        FindAllMarkers(x,only.pos=T) %>% 
        write_csv(paste0("tmp.out/",y,".csv"))

        )
      )

}
cluster_wise_dea = function(scrna_df, avg_expression_df, feature_df, target_cluster, avg_exp_scaled_thresh, stage_name){
  
  q1 = avg_expression_df[avg_expression_df$gene %in% feature_df$GeneID %>% which(),] %>% mutate(GeneID=gene)
  q2 = left_join(q1, feature_df, by = "GeneID")
  q3 = q2[c("GeneID","Name")]
  
  a = lab_DP(scrna_df, unique(q3), "eg", T)

  b = ggplot_build(a)$plot$data
  b = mutate(b, GeneID=rownames(b))

  c = subset(b, avg.exp.scaled>=avg_exp_scaled_thresh)
  c = arrange(c,desc(avg.exp.scaled))
  c = subset(c,id == target_cluster)

  d = c$features.plot %>% unique() %>% droplevels()
  d = tibble(GeneID=d,temp_name=d)

  e = left_join(d,feature_df,by="GeneID")
  e = e[-2]

  # f = unique(e$GeneID)
  # f %in% feature_df$GeneID %>% sum()

  # g = which(feature_df$GeneID %in% f)

  # check = feature_df[g,]
  # plot_df = mutate(check, Name = paste0(GeneID,"-",Name))
  # qq = lab_DP(scrna_df,plot_df, paste0(stage_name,"-",target_cluster),T)

  write_csv(e, paste0("tmp.out/",stage_name,"-",target_cluster,".csv"))

}
get_cells_with_features = function(seurat_obj, features){
  
  list_of_lists_of_cells_with_feature = map(.x = features, .f = function(x)(names(which(GetAssayData(sp48)[which(GetAssayData(sp48) %>% rownames() %in% x),] > 0))))
  return(list_of_lists_of_cells_with_feature)
  
}
cluster_wise_dot_plot = function(scrna_df, dea_list_dir, stage_name, protein_family_name, cluster_name){

  # f = unique(e$GeneID)
  # f %in% feature_df$GeneID %>% sum()
  # g = which(feature_df$GeneID %in% f)
  # check = feature_df[g,]
  # plot_df = mutate(check, Name = paste0(GeneID,"-",Name))
  # qq = lab_DP(scrna_df,plot_df, paste0(stage_name,"-",target_cluster),T)

}
subcluster_suite = function(seurat_object, res, ndims, strat, jackstraw, coord_strat){
  
  tic()
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- FindVariableFeatures(seurat_object, nfeatures=2000)
  seurat_object <- ScaleData(seurat_object)
  
  if (coord_strat == "pca"){
    
    seurat_object <- RunPCA(seurat_object, verbose = FALSE)
    
  }
  
  if (coord_strat == "ica"){
    
    seurat_object <- RunICA(seurat_object, verbose = FALSE, ica.function = "icaimax")
    
  }
  seurat_object <- FindNeighbors(seurat_object, dims = ndims)
  seurat_object <- FindClusters(seurat_object, resolution = res, verbose = FALSE, algorithm = 3)

  if (jackstraw == T){
  seurat_object = JackStraw(seurat_object, num.replicate = 200)
  seurat_object = ScoreJackStraw(seurat_object, dims = 1:20)
  }
  
  if (strat == "umap"){
    seurat_object <- RunUMAP(seurat_object, dims = ndims, n.epochs = 500, n.neighbors = 5L, min.dist = 0.05)
  } else {
    seurat_object <- RunTSNE(seurat_object, dims = ndims)
  }
  
  toc()
  print("****** Job finished. ******")
  
  return(seurat_object)
  
}
get_cell_stats = function(seurat_obj){
  
  cluster_ids = unique(Idents(seurat_obj)) 
  cell_nums = map(.x = unique(Idents(seurat_obj)), .f = function(x)(seurat_obj %>% 
                                                                      subset(idents = x) %>% 
                                                                      GetAssayData() %>% 
                                                                      ncol())
  ) %>% 
    unlist()
  cell_stats_tibble = tibble(cluster_num = cluster_ids, num_cells = cell_nums)
  return(cell_stats_tibble)
}
plot_multiple_res = function(seurat_obj, ndims, res_list){
  tic()
  plan(strategy = multiprocess, workers = 16)
  cluster_list = future_map(.x = res_list, .f = function(x)(subcluster_suite(seurat_obj, x, ndims)))
  dimplot_list = future_map(cluster_list, DimPlot)
  plots = CombinePlots(dimplot_list)
  return(plots)
  toc()
}
plot_multiple_ndims = function(seurat_obj, ndims_list, res){
  tic()
  plan(strategy = multiprocess, workers = 16)
  cluster_list = future_map(.x = ndims_list, .f = function(x)(subcluster_suite(seurat_obj, res, x)))
  dimplot_list = future_map(cluster_list, DimPlot)
  plots = CombinePlots(dimplot_list)
  return(plots)
  toc()
}
plot_smt_specificity = function(seurat_obj, smt_seurat_obj, smt){
  
  p1 = DimPlot(seurat_obj, label = T, order = T)
  p2 = FeaturePlot(seurat_obj, features = smt)
  p3 = DimPlot(smt_seurat_obj, label = T, order = T)
  p4 = FeaturePlot(smt_seurat_obj, features = smt)
  
  p1 + p2 + p3 + p4
}
plot_smt_specificity = function(orig_seurat, smt_seurat, smt){
  
  p1 = DimPlot(orig_seurat)
  p2 = FeaturePlot(orig_seurat, features = smt, order = F)
  p3 = DimPlot(smt_seurat, label = T)
  p4 = FeaturePlot(smt_seurat, features = smt, order = F)
  
  p1 + p2 + p3 + p4
}
extract_subcluster_data = function(stage, file_name, res, ndims, strat, jackstraw){
  
  stage_smts_only = subset(stage, features = c(abc_slc_3_geneid$Name, abc_slc_3_geneid$Name))
  stage_smts_only = subcluster_suite(stage_smts_only, res = res, ndims = ndims, strat = strat, jackstraw = jackstraw)
  stage_smts_only_markers = FindAllMarkers(stage_smts_only, only.pos = F)
  
  t1 = stage_smts_only_markers %>% 
    mutate(pct.diff = pct.1-pct.2) %>%
    subset(avg_log2FC > 1, p_val_adj = 0) %>%
    group_by(cluster) %>% 
    group_map(~ arrange(.x, desc(avg_log2FC)))
  
  exp_1 = map2(.x = t1, .y = c(0:(length(t1)-1)),.f = function(x,y)(add_column(x, cluster = y))) %>% 
    bind_rows()
  
  write_csv(exp_1, file_name)
  return(stage_smts_only)
}
merge_paralog_feature_names_gene_label = function(features_path, features_df){
  
  # features_path = ("~/Downloads/GSE149221_RAW/Sp1/features.tsv.gz")
  features_source = read_tsv(features_path, col_names = F)
  features_df = genes_to_be_merged
  
  list_of_paralog_dfs = features_df %>%
    mutate(source_name = Name) %>% 
    group_by(source_name) %>% 
    group_map(~ (.x))
  
  plan("multisession", workers = 8)
  list_of_paralog_indices = future_map(list_of_paralog_dfs, ~(features_source$X2 %in% .x$GeneID %>% which()))
  first_paralog_index_list = map(.x = list_of_paralog_indices, .f = function(x)(x[1]))
  # merged_paralog_name_list = map(.x = list_of_paralog_dfs, .f = function(x)(unique(x$Name))) %>% unlist()
  # merged_paralog_name_list = merged_paralog_name_list %>%
  #   str_remove_all(" \\[EC([:digit:]*.)*\\]") %>%
  #   str_remove_all("\\/") %>%
  #   str_remove_all(",") %>%
  #   str_replace_all(" ", "_")
  # features_source[first_paralog_index_list %>% unlist(),]
  merged_gene_key = bind_cols(features_source[unlist(first_paralog_index_list),]$X2, paste0("merged_gene_",1:length(features_source[unlist(first_paralog_index_list),]$X2)))
  names(merged_gene_key) = c("GeneID", "Name")
  features_source[unlist(first_paralog_index_list),]$X2 = paste0("merged_gene_",1:nrow(merged_gene_key))
  
  return(features_source)
  
  }
merge_paralog_feature_counts = function(matrix_path, features_path, features_df, key_name, path){

  # matrix_path = "~/sp_transportomics/data_sources/primary/scrna/GSE155427_RAW_modified/sp48gcm_mo/matrix.mtx.gz"
  # features_path = "~/sp_transportomics/data_sources/primary/scrna/GSE155427_RAW_modified/sp48gcm_mo/features.tsv.gz"
  # features_df = genes_to_be_merged
  # key_name = "sp48gcm_mo"
  # path = "~/sp_transportomics/data_sources/primary/scrna/GSE155427_RAW_modified/sp48gcm_mo/"
  # 
  orig_matrix = readMM(matrix_path)
  features_source = read_tsv(features_path, col_names = F)
  list_of_paralog_dfs = features_df %>%
    mutate(source_name = Name) %>% 
    group_by(source_name) %>% 
    group_map(~ (.x))

  list_of_paralog_indices = map(.x = list_of_paralog_dfs, function(x)(features_source$X2 %in% x$GeneID %>% which()))
  # list_of_paralog_indices = list_of_paralog_indices[-c(which((map(list_of_paralog_indices, length) %>% unlist()) == 0))]
  
  new_matrix = Matrix(ncol = ncol(orig_matrix), sparse = T)
  print("Performing matrix addition operations...")

  for (i in 1:length(list_of_paralog_indices)){
    if (length(list_of_paralog_indices[[i]])>1){
      summed_appended_row = Matrix(orig_matrix[list_of_paralog_indices[[i]],] %>% colSums(), nrow = 1, sparse = T)
      new_matrix = rbind(new_matrix, summed_appended_row)
    } else {
      unmodified_appended_row = orig_matrix[list_of_paralog_indices[[i]],]
      new_matrix = rbind(new_matrix, unmodified_appended_row)
    }
    
    if (i%%25 == 0){
      percent_complete = 100*(i/length(list_of_paralog_indices)) %>% round(digits = 4)
      print(paste0(percent_complete, "% of matrix operations for ", key_name, " completed..."))
    }
  }
  print(paste0("100% of matrix operations for ", key_name, " completed."))
  new_matrix = as(new_matrix, "dgCMatrix")
  rownames(new_matrix) = NULL
  new_matrix = new_matrix[-1,]
  
  first_paralog_index_list = map(.x = list_of_paralog_indices, .f = function(x)(x[1])) %>% unlist()
  other_paralog_index_list = map(.x = list_of_paralog_indices, .f = function(x)(x[-c(1)])) %>% unlist()
  print("Copying summed matrix data...")

  orig_matrix[first_paralog_index_list[1:1000],] = new_matrix[1:1000,]
  orig_matrix[first_paralog_index_list[1001:2000],] = new_matrix[1001:2000,]
  orig_matrix[first_paralog_index_list[2001:length(first_paralog_index_list)],] = new_matrix[2001:length(first_paralog_index_list),]

  gene_names_for_key = map(.x = list_of_paralog_dfs, .f = function(x)(x[1,3])) %>% unlist(use.names = F)
  gene_set_names_for_key = map(.x = list_of_paralog_dfs, .f = function(x)(x[1,4])) %>% unlist(use.names = F)
  key = features_source[first_paralog_index_list,] %>% 
    mutate(
      GeneID = paste0("merged-gene-", 
                  as.character(1:length(first_paralog_index_list))),
      Name = gene_names_for_key,
      GeneSet = gene_set_names_for_key
      )
  
  features_source[first_paralog_index_list,]$X2 = key$GeneID
  features_source = features_source[-other_paralog_index_list,]
  write_tsv(features_source, paste0(path, "features.tsv.gz"), col_names = F)
  
  orig_matrix = orig_matrix[-other_paralog_index_list,]
  writeMM(orig_matrix, paste0(path,"matrix.mtx.gz"))
  write_csv(key, paste0(key_name,"_key.csv"))
  
}
get_pos_cells = function(seurat_obj, gene_id, thresh, rev = NULL){
  
  if (rev == T){
    
    cells = which(GetAssayData(seurat_obj, slot = "counts")[gene_id,] < thresh) %>% names()
    
  } else {
    cells = which(GetAssayData(seurat_obj, slot = "counts")[gene_id,] > thresh) %>% names()
    
  }
  
  if (length(cells) > 0){
    
    pos_cells = subset(seurat_obj, cells = cells)
    pos_cells_raw_count_values = GetAssayData(pos_cells, slot = "counts")[gene_id,]
    pos_cells_normalized_cell_type_relative_avg_log2FC_values = GetAssayData(pos_cells)[gene_id,]
    pos_cells_orig_cluster_id = map(.x = cells, .f = function(x)(seurat_obj@meta.data[x,])) %>% bind_rows() %>% select(seurat_clusters)
    
    pos_cells_tibble = tibble(cells = cells, 
                              raw_counts = pos_cells_raw_count_values, 
                              normalized_cell_type_relative_avg_log2FC = pos_cells_normalized_cell_type_relative_avg_log2FC_values, 
                              orig_cluster_id = pos_cells_orig_cluster_id$seurat_clusters, 
                              gene = gene_id, 
                              stage = seurat_obj@meta.data$orig.ident %>% levels(),
                              num_cells = length(cells)
    )
  } else {
    
    pos_cells_tibble = tibble(
      cells = c(NA) ,
      raw_counts = 0, 
      normalized_cell_type_relative_avg_log2FC = 0, 
      orig_cluster_id = c(NA) %>% as.factor(), 
      gene = gene_id, 
      stage = seurat_obj@meta.data$orig.ident %>% levels(),
      num_cells = 0 %>% as.integer()
    )
  }
  return(pos_cells_tibble)
}
get_basic_stats = function(pos_hit_list, tp_cts_seuratobj) {
  
  tp_specific_hit_nums = map(pos_hit_list, length) %>% unlist()
  existing_stages = map(.x = tp_cts_seuratobj, .f = function(x)(x@meta.data$orig.ident %>% unique())) %>% unlist()
  existing_cell_types = map(.x = tp_cts_seuratobj, .f = function(x)(x@active.ident %>% levels())) %>% unlist()
  num_cells_in_cell_type_each_stage = map(.x = tp_cts_seuratobj, .f = function(x)(x@meta.data %>% nrow())) %>% unlist()
  # total_num_cells_in_stage = map(.x = urchin_3_devolist, .f = function(x)(x@meta.data %>% nrow())) %>% unlist()
  
  df = tibble(
    num_genes = tp_specific_hit_nums, 
    stage = existing_stages, 
    cell_type = existing_cell_types, 
    num_cells = num_cells_in_cell_type_each_stage
    # total_cells_in_stage = total_num_cells_in_stage
  )
  
  return(df)
  
}
get_positive_hits = function(stage_specific_devo_expression_matrix, thresh){
  list_of_list_of_genes_with_positive_hits = list()
  
  
  for (i in 1:nrow(stage_specific_devo_expression_matrix)) {
    list_of_list_of_genes_with_positive_hits[[rownames(stage_specific_devo_expression_matrix)[i]]] = which(stage_specific_devo_expression_matrix[i,] > thresh) %>% names()
    # message = paste0(i,": ",rownames(stage_specific_devo_expression_matrix)[i]," complete.")
    
  }
  
  expressed_genes = map(list_of_list_of_genes_with_positive_hits,length) > 0
  list_of_list_of_genes_with_positive_hits = list_of_list_of_genes_with_positive_hits[which(expressed_genes)]
  return(list_of_list_of_genes_with_positive_hits)
  
}
get_basic_stats_2 = function(seurat_obj, gene_list, thresh){
  
  seurat_obj_idents = 
    seurat_obj %>% 
    levels()
  
  seurat_objects_by_ident = 
    map(.x = seurat_obj_idents, 
        .f = function(x)(subset(seurat_obj, idents = x) 
        ) 
    )
  
  seurat_objects_by_ident_gene_list_filtered =
    map(.x = seurat_objects_by_ident, 
        .f = function(x)(GetAssayData(x, slot = "count")[which(GetAssayData(x) %>% rownames() %in% gene_list$GeneID), ]) 
    )
  
  pos_hit_list =
    map(.x = seurat_objects_by_ident_gene_list_filtered, .f = function(x) (get_positive_hits(x, thresh)))
  
  get_basic_stats(pos_hit_list, seurat_objects_by_ident)
  
}
anno_filter = function(annofile, term){
  
  res = annofile[annofile$product %>% str_which(term),]
  res = mutate(res, used_search_term = term)
  return(res)
  
}
parse_kegg_brite = function(kegg_json_file_path){
  awful_kegg_json = tibble(L2_l1 = "L2_l1",
                           L2_l2 = "L2_l2",
                           L3_l1 = "L3_l1",
                           L3_l2 = "L3_l2",
                           L3_l3 = "L3_l3",
                           L4_l1 = "L4_l1",
                           L4_l2 = "L4_l2",
                           L4_l3 = "L4_l3",
                           L4_l4 = "L4_l4",
                           L5_l1 = "L5_l1",
                           L5_l2 = "L5_l2",
                           L5_l3 = "L5_l3",
                           L5_l4 = "L5_l4",
                           L5_l5 = "L5_l5")
  
  kegg_sp_gene_set = read_json(kegg_json_file_path)
  kegg_sp_gene_set_tibble = kegg_sp_gene_set %>% as_tibble()
  kegg_sp_gene_set_tibble_processed = kegg_sp_gene_set_tibble[2]
  
  l1_children_names = map(.x = kegg_sp_gene_set_tibble_processed$children, .f = function(x)(x$name))
  l1_children_lengths = map(.x = kegg_sp_gene_set_tibble_processed$children, .f = function(x)(x$children %>% length()))
  l1_children = map(.x = kegg_sp_gene_set_tibble_processed$children, .f = function(x)(x$children))
  
  for (g in 1:length(l1_children_names)){
    l2_children_names = map(.x = l1_children[[g]], .f = function(x)(x$name))
    l2_children_lengths = map(.x = l1_children[[g]], .f = function(x)(x$children %>% length()))
    
    awful_kegg_json = add_row(awful_kegg_json, 
                              L2_l2 = l2_children_names %>% unlist(),
                              L2_l1 = l1_children_names[[g]]
    )
    
    l2_children = map(.x = l1_children[[g]], .f = function(x)(x$children))
    l2_non_null = which(map(l2_children, is_null) == F)
    
    # print(paste0("g:", g))
    

    for (j in l2_non_null){
      l3_children_names = map(.x = l2_children[[j]], .f = function(x)(x$name))
      l3_children_lengths = map(.x = l2_children[[j]], .f = function(x)(x$children %>% length()))
      
      awful_kegg_json = add_row(awful_kegg_json, 
                                L3_l3 = l3_children_names %>% unlist(),
                                L3_l2 = l2_children_names[[j]],
                                L3_l1 = l1_children_names[[g]]
      )
      
      l3_children = map(.x = l2_children[[j]], .f = function(x)(x$children))
      l3_non_null = which(map(l3_children, is_null) == F)
      
      # print(paste0("j:", j))
      
      for (i in l3_non_null){
        l4_children_names = map(.x = l3_children[[i]], .f = function(x)(x$name))
        awful_kegg_json = add_row(awful_kegg_json, 
                                  L4_l4 = l4_children_names %>% unlist(), 
                                  L4_l3 = l3_children_names[[i]],
                                  L4_l2 = l2_children_names[[j]],
                                  L4_l1 = l1_children_names[[g]]
        )
        
        l4_children = map(.x = l3_children[[i]], .f = function(x)(x$children))
        l4_non_null = which(map(l4_children, is_null) == F)
        
        # print(paste0("i:", i))
        
        for (k in l4_non_null){
          l5_children_names = map(.x = l4_children[[k]], .f = function(x)(x$name))
          awful_kegg_json = add_row(awful_kegg_json,
                                    L5_l5 = l5_children_names %>% unlist(),
                                    L5_l4 = l4_children_names[[k]], 
                                    L5_l3 = l3_children_names[[i]],
                                    L5_l2 = l2_children_names[[j]],
                                    L5_l1 = l1_children_names[[g]]
          )
        
          # print(paste0("k:", k))
        
      }
    }
    }
  }
  
    awful_kegg_json = awful_kegg_json[-1,]
  return(awful_kegg_json)
  
}
in_urchin = function(seurat_obj, gene_id_list){
  
  res = which(gene_id_list %in% rownames(seurat_obj))
  return(res)
  
}
process_sp_kegg_geneset = function(parsed_and_filtered_sp_kegg_geneset, unmerged_feature_df){
  
  res1 = parsed_and_filtered_sp_kegg_geneset[colSums(!is.na(parsed_and_filtered_sp_kegg_geneset)) > 0]
  res2 = res1[[length(res1)]] %>% str_split("; ", simplify = T) %>% as_tibble()
  res3 = res2 %>% mutate(GeneID = paste0("LOC", str_extract(res2[[1]], "[:digit:]*")))
  res4 = tibble(GeneID = res3$GeneID, Name = res3$V2, additional_info = res3$V1)
  res5 = which(res4$GeneID %in% rownames(unmerged_feature_df))
  res6 = res4[res5,]
  
  return(res6)
  
}
make_merged_cell_stats = function(mappable_list_of_seurat_objects, gene_set, gene_set_name, thresh){
  res = map(
    .x = mappable_list_of_seurat_objects,
    .f = function(x)(get_basic_stats_2(x, gene_set, thresh = thresh))
  ) %>%
    bind_rows()
  
  util_rate = res$num_genes/(nrow(gene_set))
  
  res = mutate(
    res, 
    gene_set = gene_set_name, 
    cell_type_grouped_state = "unmerged", 
    gene_util_rate = util_rate)
  
  print(
    paste0(
      "Completed merging statistics for gene set: ", 
      gene_set_name)
  )
  return(res)
}
search_sp4.2_annotations = function(term){
  
  spuranno_v2 = spuranno %>% as_tibble(spuranno)
  spuranno_v3 = spuranno_v2$product
  spuranno_v4 = str_which(spuranno_v3, term)
  spuranno_v5 = spuranno_v2[spuranno_v4,] %>% select(gene, product) %>% unique()
  spuranno_v6 = spuranno_v2[spuranno_v4,] %>% select(gene) %>% unique()
  spuranno_v7 = tibble(GeneID = spuranno_v6$gene, Name = spuranno_v6$gene)
  spuranno_v8 = map(.x = map(.x = spuranno_v7$GeneID, .f = function(x)(subset(spuranno_v5, gene == x))), .f = function(x)(x[1,])) %>% 
    bind_rows() %>%
    mutate(GeneID = gene, Name = product)
  
  return(spuranno_v8)
}
get_positive_hits_per_cell_type = function(seurat_obj, gene_set, thresh, relative){
  
  cell_types = Idents(seurat_obj) %>% levels()
  stage = seurat_obj@meta.data$orig.ident %>% levels()
  num_all_genes_in_gene_set = seurat_obj %>% 
    rownames() %>% 
    length()
  
  # for one stage, create separate seurat_objs by cell type
  t1.b = map(
    .x = cell_types,
    .f = function(x)(
      subset(seurat_obj, idents = x)
    )
  )
  
  if(relative == T){
    assay_data_choice = "scale.data"
  } else{
    assay_data_choice = "counts"
  }
  # for one stage, for each cell type, get expression count matrix data
  t1.c = map(
    .x = t1.b,
    .f = function(x)(
      GetAssayData(
        x,
        slot = assay_data_choice
      )
    )
  )
  
  # for one stage, for each cell type, get cells for genes with any positive hits
  t2 = map(
    .x = t1.c,
    .f = function(x)(
      get_positive_hits(x, thresh = thresh)
    )
  )
  
  # get the number of smts present in each cell type with at least one transcript count
  t3 = map(
    .x = t2,
    .f = function(x)(
      length(x)
    )
  ) %>% 
    unlist()
  
  res = tibble(
    Stage = stage,
    "Cell type" = cell_types,
    "Gene set" = gene_set,
    "Number present" = t3,
    "Total existing" = num_all_genes_in_gene_set
  )
  
  return(res)
}
get_transcript_raw_counts_by_geneset = function(genesets_key, seurat_obj){
  
  stage = 
    seurat_obj@meta.data$orig.ident %>%
    levels()
  
  geneset_names = 
    genesets_key$GeneSet %>% 
    unique() %>% 
    sort()
  
  total_transcript_count = 
    GetAssayData(seurat_obj, slot = "counts") %>% 
    sum()
  
  geneset_key_grouped_by_geneset = 
    genesets_key %>% 
    group_by(GeneSet) %>% 
    group_map(~ (.x),
              .keep = TRUE)
  
  seurat_obj_subsets_by_geneset = map(
    .x = geneset_key_grouped_by_geneset, 
    .f = function(x)(
      subset(
        seurat_obj,
        features = x$GeneID
      )
    )
  )
  
  transcript_counts_by_geneset = map(
    .x = seurat_obj_subsets_by_geneset,
    .f = function(x)(
      GetAssayData(
        x,
        slot = "counts") %>% 
        sum()
    )
  )
  
  res = tibble(
    stage = stage,
    GeneSet = geneset_names,
    geneset_raw_transcript_counts = transcript_counts_by_geneset %>% unlist(),
    total_raw_transcript_counts = total_transcript_count,
    geneset_expression_ratio = geneset_raw_transcript_counts/total_raw_transcript_counts
    
  )
  
  return(res)
}
get_num_pos_cells_multiple_genes = function(seurat_obj, gene_list, thresh, rev = NULL){
  
  res1 = map(
    .x = gene_list, 
    .f = function(x)(
      get_pos_cells(seurat_obj, 
                    x, 
                    thresh = thresh, 
                    rev = rev)
    )
  )  
  
  res2 = map(
    .x = res1,
    .f = function(x)(
      tibble(
        gene = x$gene %>% unique(),
        stage = x$stage %>% unique(),
        num_cells = x$num_cells %>% unique()
      )
    )
  )
  return(res2)
  
}
get_overlap_percentages_one_stage_one_geneset = function(seurat_obj, gene_set){
  
  # seurat_obj = urchin_5_clustered_seurat_objs[[1]]
  # gene_set = urchin_5_smts_final_cleaned_with_mito_smts %>% mutate(GeneSet = "smts")
  
  b = Idents(seurat_obj) %>% levels()
  
  
  # subset a single seurat obj by ident or "cell type"
  c = map(
    .x = b,
    .f = function(x)(
      subset(
        seurat_obj,
        ident = x)
    )
  )
  
  # subset each seurat obj in list of seurat objs by GeneIDs of a single gene set
  d = map(
    .x = c,
    .f = function(x)(
      subset(
        x,
        features = gene_set$GeneID
      )
    )
  )
  
  # get transcript count data from each seurat obj in list of seurat objs that have been subsetted by GeneIDs in a single gene set
  e = map(
    .x = d,
    .f = function(x)(
      GetAssayData(x, slot = "counts")
    )
  )
  
  # two things achieved here: get cell IDs for cells expressing each gene in gene set. 
  # not only will this tell us #1 which and how many cells express a gene, but also #2 how many genes have any expression
  # at all. remember, this is still being done for each cell type.
  f = map(
    .x = e,
    .f = function(x)(
      get_positive_hits(
        x,
        0
      )
    )
  )
  
  # get the names of the genes that have any number of cells expressing each gene for each cell type.
  g = map(
    .x = f,
    .f = function(x)(
      names(x)
    )
  )
  
  # create Venn object- requires a list of vectors to be compared
  g = Venn(g)
  
  # find pairwise differences
  h = discern_pairs(g)
  
  i = map(
    .x = h,
    .f = function(x)(
      length(x)/nrow(gene_set)
    )
  )
  
  j = names(i)
  k = unlist(i, use.names = F)
    
  l = tibble(
    gene_set = gene_set$GeneSet %>% unique(),
    comparison_pair = j,
    difference_percentage = k,
    mean_pairwise_difference = mean(k),
    stage = seurat_obj@meta.data$orig.ident %>% levels()
  )
  return(l)
  
}
get_average_overlap_percentages_one_stage_multiple_genesets = function(seurat_obj, gene_set_list){
  
  res = map(.x = gene_set_list, 
            .f = function(x)(
              get_overlap_percentages_one_stage_one_geneset(seurat_obj, x)
            )
  ) %>% bind_rows()
  
  res = res %>% mutate(Stage = seurat_obj@meta.data$orig.ident %>% levels())
  print(paste0(seurat_obj@meta.data$orig.ident %>% levels(), " done."))
  return(res)
  
}
get_positive_hits_for_geneset_per_cell_type = function(seurat_obj, gene_set, gene_set_name, thresh, relative, from_kegg = NULL, missing_germline_smts = NULL){
  
  # seurat_obj = urchin_5_spmb_cell_type_labelled
  # gene_set = urchin_5_smts_final_cleaned
  # gene_set_name = "smts"
  # thresh = 0
  # relative = F
  # from_kegg = F

  
  if (from_kegg == T){
    geneset_subset = subset(urchin_5_keys[[1]], GeneSet == gene_set_name)
  } else {
    geneset_subset = gene_set
  }
  
  seurat_obj = subset(
    seurat_obj, 
    features = geneset_subset$GeneID
  )
  
  cell_types = Idents(seurat_obj) %>%
    levels()
  
  stage = seurat_obj@meta.data$orig.ident %>%
    levels()
  
  num_all_genes_in_gene_set = seurat_obj %>% 
    rownames() %>% 
    length()
  
  # for one stage, create separate seurat_objs by cell type
  t1.b = map(
    .x = cell_types,
    .f = function(x)(
      subset(seurat_obj, idents = x)
    )
  )
  
  if(relative == T){
    assay_data_choice = "scale.data"
  } else{
    assay_data_choice = "counts"
  }
  # for one stage, for each cell type, get expression count matrix data
  t1.c = map(
    .x = t1.b,
    .f = function(x)(
      GetAssayData(
        x,
        slot = assay_data_choice
      )
    )
  )
  
  # for one stage, for each cell type, get cells for genes with any positive hits
  t2 = map(
    .x = t1.c,
    .f = function(x)(
      get_positive_hits(x, thresh = thresh)
    )
  )
  
  # t2_germline_test = map(t2, names)
  # t2_non_germline = t2_germline_test[-c(8)]
  # t2_germline = t2_germline_test[[8]]
  # 
  # stage_germline_absent_smts = map(.x = t2_non_germline, .f = function(x)(setdiff(x, t2_germline))) %>% 
  #   unlist() %>% 
  #   unique()
  # 
  # stage_germline_absent_smts_tibble = tibble(
  #   germline_missing_smts = stage_germline_absent_smts,
  #   stage = stage
  # )
  # 
  # get the number of smts present in each cell type with at least one transcript count
  t3 = map(
    .x = t2,
    .f = function(x)(
      length(x)
    )
  ) %>% 
    unlist()
  
  res = tibble(
    Stage = stage,
    "Cell type" = cell_types,
    "Gene set" = gene_set_name,
    "Number present" = t3,
    "Total existing" = num_all_genes_in_gene_set
  )
  
  
  
  if (missing_germline_smts == T){
    return(stage_germline_absent_smts_tibble)
  } else {
    return(res)
  }
  
}
search_purp = function(term){
  
  spuranno_res_1 = 
    spuranno[spuranno$product %>% str_which(term),] %>% 
    as_tibble() %>% 
    select(gene, product) %>% 
    transmute(GeneID = gene, Name = product) %>%
    unique()
  
  unique_gene_names = spuranno_res_1 %>% 
    group_by(GeneID) %>% 
    group_map(~.x[1,], .keep = T) %>% 
    bind_rows()
  
  spuranno_res_2 = tibble(
    GeneID = spuranno_res_1$GeneID %>% unique(),
    Name = spuranno_res_1$GeneID %>% unique()
  ) %>% 
    right_join(unique_gene_names, by = "GeneID") %>% 
    transmute(
      GeneID = GeneID,
      Name = Name.y
      )

    return(spuranno_res_2)
}
label_marked_cells = function(seurat_obj, marker_list, cell_type_labels, thresh){
  
  marked_cells_df = map(
    .x = marker_list,
    .f = function(x)(
      get_pos_cells(
        seurat_obj, 
        gene_id = x, 
        thresh = thresh, 
        rev = F)
    )
  ) %>% 
    bind_rows() 
  
  marked_cells = marked_cells_df$cells %>% unique()
  print(marked_cells_df)
  Idents(seurat_obj, cells = marked_cells) = cell_type_labels
  
  return(seurat_obj)
}
get_expression_data_from_cells_with_gradients_of_marker_expression = function(seurat_obj, gene){
  
  gene_pos_cells = get_pos_cells(seurat_obj, gene_id = gene, thresh = 0, rev = F)
  gene_transcript_count_quantiles = quantile(gene_pos_cells$raw_counts, c(0.95, 0.90, 0.85, 0.80, 0.75))
  gene_transcript_count_quantiles_tibble = tibble(
    
    quantile = names(gene_transcript_count_quantiles),
    transcript_count = gene_transcript_count_quantiles
    
  )
  
  mapped_gene_transcript_count_quantiles_tibble = 
    map(
      .x = names(gene_transcript_count_quantiles), 
      .f = function(x)(
        subset(
          gene_transcript_count_quantiles_tibble,
          quantile == x)
      )
    )
  
  gene_pos_cells_95_100_quantile = gene_pos_cells %>% subset(raw_counts >= mapped_gene_transcript_count_quantiles_tibble[[1]]$transcript_count)
  gene_pos_cells_90_95_quantile = gene_pos_cells %>% subset(raw_counts >= mapped_gene_transcript_count_quantiles_tibble[[2]]$transcript_count & raw_counts < mapped_gene_transcript_count_quantiles_tibble[[1]]$transcript_count) 
  gene_pos_cells_85_90_quantile = gene_pos_cells %>% subset(raw_counts >= mapped_gene_transcript_count_quantiles_tibble[[3]]$transcript_count & raw_counts < mapped_gene_transcript_count_quantiles_tibble[[2]]$transcript_count)
  gene_pos_cells_80_85_quantile = gene_pos_cells %>% subset(raw_counts >= mapped_gene_transcript_count_quantiles_tibble[[4]]$transcript_count & raw_counts < mapped_gene_transcript_count_quantiles_tibble[[3]]$transcript_count)
  gene_pos_cells_75_80_quantile = gene_pos_cells %>% subset(raw_counts >= 1 & raw_counts < mapped_gene_transcript_count_quantiles_tibble[[4]]$transcript_count) 
  
  gene_pos_cells_quantile_list = list(
    
    gene_pos_cells_95_100_quantile,
    gene_pos_cells_90_95_quantile,
    gene_pos_cells_85_90_quantile,
    gene_pos_cells_80_85_quantile,
    gene_pos_cells_75_80_quantile
    
  )
  
  gene_pos_cells_quantile_names_list = list(
    "_95_100_quantile",
    "_90_95_quantile",
    "_85_90_quantile",
    "_80_85_quantile",
    "_75_80_quantile"
  )
  
  gene_pos_cells_nrows = map(gene_pos_cells_quantile_list, nrow)
  
  for (i in 1:length(gene_pos_cells_quantile_list)){
    
    if(gene_pos_cells_nrows[[i]] > 0){
      Idents(seurat_obj, cells = gene_pos_cells_quantile_list[[i]]$cells) = paste0(gene, gene_pos_cells_quantile_names_list[[i]])
    }
  }
  
  return(seurat_obj)
  
}
find_new_cluster_cell_idents_in_old_cluster_cell_idents = function(new_seurat_obj, old_seurat_obj){
  
  new_seurat_obj_idents = levels(Idents(new_seurat_obj))
  new_seurat_obj_idents_subsetted_list = map(.x = new_seurat_obj_idents,
                                             .f = function(x)(subset(new_seurat_obj, 
                                                                     idents = x)
                                             )
  )
  cells_in_new_seurat_obj_for_each_ident_in_old_seurat_obj = map(.x = new_seurat_obj_idents_subsetted_list, 
                                                                 .f = function(x)(subset(old_seurat_obj,
                                                                                         cells = Cells(x)
                                                                 )
                                                                 )
  )
  
  cell_stats_of_cells_in_new_seurat_obj_for_each_ident_in_old_seurat_obj = map(.x = cells_in_new_seurat_obj_for_each_ident_in_old_seurat_obj, 
                                                                               .f = function(x)(
                                                                                 get_cell_stats(x)
                                                                               )
  )
  
  cells_stats_of_old_seurat_obj = get_cell_stats(old_seurat_obj)
  cells_stats_of_old_seurat_obj$cluster_num = levels(cells_stats_of_old_seurat_obj$cluster_num)
  
  cell_stats_of_cells_in_new_seurat_obj_for_each_ident_in_old_seurat_obj_with_new_cluster_id_origin_column = map2(.x = cell_stats_of_cells_in_new_seurat_obj_for_each_ident_in_old_seurat_obj,
                                                                                                                  .y = new_seurat_obj_idents,
                                                                                                                  .f = function(x,y)(
                                                                                                                    mutate(x,
                                                                                                                           cluster_origin = as.character(y)
                                                                                                                    )
                                                                                                                  )
  )
  
  for (j in 1:length(cell_stats_of_cells_in_new_seurat_obj_for_each_ident_in_old_seurat_obj_with_new_cluster_id_origin_column)){
    missing_cell_identities = which(cells_stats_of_old_seurat_obj$cluster_num %in% cell_stats_of_cells_in_new_seurat_obj_for_each_ident_in_old_seurat_obj_with_new_cluster_id_origin_column[[j]]$cluster_num == F)
    
    if (length(missing_cell_identities >= 1)){
      for (i in 1:length(missing_cell_identities)){
        cell_stats_of_cells_in_new_seurat_obj_for_each_ident_in_old_seurat_obj_with_new_cluster_id_origin_column[[j]] =
          add_row(cell_stats_of_cells_in_new_seurat_obj_for_each_ident_in_old_seurat_obj_with_new_cluster_id_origin_column[[j]],
                  cluster_num = cells_stats_of_old_seurat_obj$cluster_num[missing_cell_identities[i]],
                  num_cells = 0,
                  cluster_origin = cell_stats_of_cells_in_new_seurat_obj_for_each_ident_in_old_seurat_obj_with_new_cluster_id_origin_column[[j]]$cluster_origin %>% unique())
        
      }
      
    } else {print("mop")}
    
    
    
  }
  res = map(.x = cell_stats_of_cells_in_new_seurat_obj_for_each_ident_in_old_seurat_obj_with_new_cluster_id_origin_column,
            .f = function(x)(arrange(x, desc(cluster_num))))
  # res = bind_rows(cell_stats_of_cells_in_new_seurat_obj_for_each_ident_in_old_seurat_obj_with_new_cluster_id_origin_column)
  # res$cluster_origin = as.factor(res$cluster_origin)
  return(res)
  
}

