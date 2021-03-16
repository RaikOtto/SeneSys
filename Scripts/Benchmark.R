run_benchmark = function(
    dataset_query,
    dataset_training,
    algorithm,
    type,
    transcriptome_file,
    visualization_file,
    path_benchmark_files,
    confidence_threshold,
    high_threshold,
    low_threshold,
    path_benchmark_files_dec_res
){
    
    ### prep
    
    algorithm = str_to_lower(algorithm)
    transcriptome_data = read.table(transcriptome_file, sep ="\t",header = T, row.names = 1, stringsAsFactors = F)
    colnames(transcriptome_data) = str_replace_all(colnames(transcriptome_data),pattern="^X","")
    meta_data = meta_info[colnames(transcriptome_data),]
    if(  !nrow(meta_data) == ncol(transcriptome_data))
        stop("Inconclusive meta data and transcriptome file dimensions")
    row_names = as.character(rownames(transcriptome_data))
    col_names = as.character(colnames(transcriptome_data))
    transcriptome_data = matrix(as.double(as.character(unlist(transcriptome_data))),ncol = length(col_names))
    rownames(transcriptome_data) = row_names
    colnames(transcriptome_data) = col_names
    
    visualization_data = read.table(visualization_file, sep ="\t",header = T, row.names = 1, stringsAsFactors = F)
    colnames(visualization_data) = str_replace_all(colnames(visualization_data),pattern="^X","")
    
    name_query_data = tail(as.character(unlist(str_split(transcriptome_file,pattern ="/"))),1)
    name_query_data = str_replace_all(name_query_data,pattern =".tsv","")
    
    name_training_data = tail(as.character(unlist(str_split(dataset_training,pattern ="_"))),1)
    
    name_algorithm_datatype = tail(as.character(unlist(str_split(path_benchmark_files,pattern ="/"))),1)
    name_algorithm_datatype = str_replace_all(name_algorithm_datatype,pattern ="(.tsv)|(Benchmark_results.)","")
    name_algorithm = head(as.character(unlist(str_split(name_algorithm_datatype,pattern ="\\."))),1)
    name_datatype = tail(as.character(unlist(str_split(name_algorithm_datatype,pattern ="\\."))),1)
    
    identifier = tail( as.character(unlist(str_split(dataset_training, pattern = "_"))), 1)
    ki_index = which(rownames(transcriptome_data) == "MKI67")
    
    ###
    
    model_path = paste0(
        c(
            "~/Deko_Projekt/Data/Bench_data/Models/",
            dataset_query,
            type,
            identifier,
            algorithm,
            "RDS"
        ),
        collapse = "."
    )
    model_path = str_replace_all(model_path, pattern = "/\\.","/")
    
    ###
    if (
        file.exists(str_replace(path_benchmark_files_dec_res, pattern = ".dec_res.tsv",".dec_res.RDS")))
    {
        #benchmark_results_t = read.table(path_benchmark_files,sep="\t",stringsAsFactors = F,header = T)
        deconvolution_results = readRDS(
            str_replace(path_benchmark_files_dec_res, pattern = ".dec_res.tsv",".dec_res.RDS"))

    } else {
    
        deconvolution_results = Deconvolve_transcriptome(
            transcriptome_data = transcriptome_data,
            deconvolution_algorithm = algorithm,
            models = dataset_training,
            nr_permutations = 1000,
            output_file = ""
        )
        
        if ( type == "ductal")
            deconvolution_results = deconvolution_results[
                grep(deconvolution_results$model, pattern = "Ductal", value = F) ,]
        if ( type == "hisc")
            deconvolution_results = deconvolution_results[
                grep(deconvolution_results$model, pattern = "Hisc", value = F) ,]
        
        if (algorithm == "bseqsc")
            deconvolution_results$P_value = as.double(deconvolution_results$P_value)
        
        col_names = colnames(deconvolution_results)
        width = ncol(deconvolution_results)
        height = nrow(deconvolution_results) 
        res_table = matrix(unlist(deconvolution_results),ncol = width, nrow = height)
        colnames(res_table) = col_names
        res_table= as.data.frame(res_table)

        saveRDS(
            deconvolution_results,
            str_replace(path_benchmark_files_dec_res, pattern = ".dec_res.tsv",".dec_res.RDS")
        )
        
        write.table(
            deconvolution_results,
            path_benchmark_files_dec_res,
            sep="\t",
            quote =F, row.names = F
        )
    }
    
    #return(deconvolution_results)
    
    #if (type == "hisc")
    #    deconvolution_results = deconvolution_results[grep(deconvolution_results$model, pattern = "Hisc", ignore.case = T),]
    #if (type == "ductal")
    #    deconvolution_results = deconvolution_results[grep(deconvolution_results$model, pattern = "ductal",ignore.case = T),]
    
    meta_data = meta_info[ rownames(deconvolution_results),]
    #meta_data = meta_info[ deconvolution_results$Sample,]
    #rownames(deconvolution_results) = deconvolution_results$Sample

    if (sum( meta_data[rownames(deconvolution_results),"Grading"] != "") > 0){
        
        meta_data = meta_data[meta_data$Grading!="",]
        deconvolution_results = deconvolution_results[rownames(meta_data),]
        deconvolution_results$Grading = meta_data$Grading
        visualization_data  = visualization_data[,rownames(meta_data)]
        transcriptome_data  = transcriptome_data[,rownames(meta_data)]
        
    }
    #write.table(deconvolution_results,"~/Deco/Results/Cell_fraction_predictions/RepSet_Cibersort_Baron.tsv",sep ="\t", row.names =T , quote=F)
    
    if( length(ki_index) != 0 ){
        
        deconvolution_results[,"MKI67"] = rep(0,nrow(deconvolution_results))
        deconvolution_results[,"MKI67"] = log(as.double(transcriptome_data[ki_index[1],])+1)
        
    } else {
        
        deconvolution_results[,"MKI67"] = as.double(meta_data$KI67)
    }
    
    deconvolution_results$Strength_de_differentiation = 0
    deconvolution_results$Confidence_score_dif = 0
    deconvolution_results$Subtype = ""

    vis_mat = create_visualization_matrix(
        visualization_data = visualization_data,
        deconvolution_results = deconvolution_results,
        confidence_threshold = confidence_threshold,
        high_threshold = high_threshold,
        low_threshold = low_threshold
    )
    
    ### results parsing
    
    graphics_path_heatmap = paste("~/Deko_Projekt/Results/Images",algorithm, sep ="/")
    graphics_path_heatmap = paste(graphics_path_heatmap,"Heatmap", sep ="/")
    if (! dir.exists(graphics_path_heatmap))
        dir.create(graphics_path_heatmap)
    graphics_path_heatmap = paste(graphics_path_heatmap,name_datatype, sep ="/")
    if (! dir.exists(graphics_path_heatmap))
        dir.create(graphics_path_heatmap)
    graphics_path_heatmap = paste(graphics_path_heatmap,name_query_data, sep ="/")
    if (! dir.exists(graphics_path_heatmap))
        dir.create(graphics_path_heatmap)
    graphics_path_heatmap = paste(graphics_path_heatmap,paste0(name_training_data,".pdf"),sep = "/")
    
    ### survival curve
    
    if (sum( meta_data$OS_Tissue != "") != 0){

        vis_mat = vis_mat[rownames(deconvolution_results),]
        vis_mat$OS_Tissue = as.double(str_replace_all(meta_data$OS_Tissue, pattern = ",", "\\."))
        vis_mat$OS_Tissue[is.na(vis_mat$OS_Tissue)] = 1
        vis_mat$Grading = meta_data$Grading
        vis_mat$Zensur = meta_data$Zensur

        graphics_path_survival = paste("~/Deko_Projekt/Results/Images",algorithm, sep ="/")
        graphics_path_survival = paste(graphics_path_survival,"Survival", sep ="/")
        if (! dir.exists(graphics_path_survival))
            dir.create(graphics_path_survival)
        graphics_path_survival = paste(graphics_path_survival,name_datatype, sep ="/")
        if (! dir.exists(graphics_path_survival))
            dir.create(graphics_path_survival)
        graphics_path_survival = paste(graphics_path_survival,name_query_data, sep ="/")
        if (! dir.exists(graphics_path_survival))
            dir.create(graphics_path_survival)
        
        #if ( "hisc" %in% colnames(deconvolution_results)){
        if ( sum( deconvolution_results$hisc ) > 0 ){
            
            hisc = deconvolution_results$hisc
            agg = aggregate(deconvolution_results$hisc, FUN = mean, by = list(meta_info[rownames(deconvolution_results),"Grading"]))
            
            hisc[hisc <= agg$x[3] ] = "low"
            hisc[hisc != "low"] = "high"
            
            data = vis_mat[,c("OS_Tissue","Zensur")]
            data = cbind(data, hisc)
            data = data[!is.na(data$Zensur),]

            graphics_path_survival_hisc = paste(graphics_path_survival,paste0(paste(name_training_data,"hisc",sep="_"),".pdf"),sep = "/")
            
            fit = survival::survfit( survival::Surv( as.double(data$OS_Tissue), data$Zensur ) ~ data$hisc, data = data)
            surv_hisc = survminer::surv_pvalue(fit, data = data)$pval

            pdf(graphics_path_survival_hisc,onefile = FALSE)#,width="1024px",height="768px")
                print(survminer::ggsurvplot(fit, data = data, risk.table = T, pval = T, censor.size = 10))
            dev.off()
            #surv_ductal = 1
        
        }

        graphics_path_survival_ductal = paste(graphics_path_survival,paste0(paste(name_training_data,"ductal",sep="_"),".svg"),sep = "/")
        
        #deconvolution_results$
        
        ductal_values = deconvolution_results$ductal
        agg = aggregate(ductal_values, FUN = mean, by = list(meta_info[rownames(deconvolution_results),"Grading"]))
        
        ductal = rep("high",length(ductal_values))
        #ductal[ductal_values <= mean( ductal_values )]   = "low"
        #ductal[ductal <= median( ductal )]   = "low"
        ductal[ductal_values <= tail( agg$x, 1 )]   = "low"
        #ductal[ductal <= ( ( tail( agg$x, 1 ) + tail( agg$x, 2 )[1] ) / 2 ) ]    = "low"
        
        data = vis_mat[,c("OS_Tissue","Zensur")]
        data = cbind(data, ductal)
        data = data[!is.na(data$Zensur),]
        
        fit = survival::survfit( survival::Surv( as.double(data$OS_Tissue), data$Zensur ) ~ data$ductal)
        surv_ductal = survminer::surv_pvalue(fit, data = data)$pval
        surv_ductal
        
        svg(filename = graphics_path_survival_ductal, width = 10, height = 10)
        print(survminer::ggsurvplot(fit, data = data, risk.table = F, pval = T, censor.size = 10))
        dev.off()

        #surv_hisc = 1

        
        # mki67
        
        ki67_values = deconvolution_results$MKI67
        agg = aggregate(ki67_values, FUN = mean, by = list(meta_info[rownames(deconvolution_results),"Grading"]))
        
        ki67 = rep("high",length(ki67_values))
        #ki67[ki67_values <= tail( agg$x, 1 )]   = "low"
        ki67[ki67_values <= mean( ki67_values)]   = "low"
        
        data = vis_mat[,c("OS_Tissue","Zensur")]
        data = cbind(data, ki67)
        data = data[!is.na(data$Zensur),]

        surv_mki67 = survminer::surv_pvalue(survival::survfit( survival::Surv( as.double(data$OS_Tissue) ) ~ data$ki67), data = data, method = "survdiff")$pval
        if ( sum(ki67 == "high") < 5) surv_mki67 = 1

        graphics_path_survival_mki67 = paste(graphics_path_survival,paste0(paste(name_training_data,"mki67",sep="_"),".pdf"),sep = "/")
        
        fit = survival::survfit( survival::Surv( as.double(data$OS_Tissue), data$Zensur ) ~ data$ki67)
        surv_mki67 = survminer::surv_pvalue(fit, data = data)$pval
        
        pdf(graphics_path_survival_mki67,onefile = FALSE)#,width="1024px",height="768px")
            print(survminer::ggsurvplot(fit, data = data, risk.table = F, pval = T, censor.size = 10))
        dev.off()
        

    } else {
        surv_mki67 = surv_ductal = surv_hisc = 1
    }
    
    ###
    
    graphics_path_pca = paste("~/Deko_Projekt/Results/Images",algorithm, sep ="/")
    graphics_path_pca = paste(graphics_path_pca,"PCA", sep ="/")
    if (! dir.exists(graphics_path_pca))
        dir.create(graphics_path_pca)
    graphics_path_pca = paste(graphics_path_pca,name_datatype, sep ="/")
    if (! dir.exists(graphics_path_pca))
        dir.create(graphics_path_pca)
    graphics_path_pca = paste(graphics_path_pca,name_query_data, sep ="/")
    if (! dir.exists(graphics_path_pca))
        dir.create(graphics_path_pca)
    graphics_path_pca = paste(graphics_path_pca,paste0(name_training_data,".pdf"),sep = "/")
    
    ### MKI-67 to ratio correlation
    
    if( length(ki_index) != 0 ){
        
        deconvolution_results$MKI67 = transcriptome_data[ki_index,rownames(deconvolution_results)]
        
        graphics_path_mki67 = paste("~/Deko_Projekt/Results/Images",algorithm, sep ="/")
        graphics_path_mki67 = paste(graphics_path_mki67,"MKI67", sep ="/")
        if (! dir.exists(graphics_path_mki67))
            dir.create(graphics_path_mki67)
        graphics_path_mki67 = paste(graphics_path_mki67,name_datatype, sep ="/")
        if (! dir.exists(graphics_path_mki67))
            dir.create(graphics_path_mki67)
        graphics_path_mki67 = paste(graphics_path_mki67,name_query_data, sep ="/")
        if (! dir.exists(graphics_path_mki67))
            dir.create(graphics_path_mki67)
        graphics_path_mki67 = paste(graphics_path_mki67,paste0(name_training_data,".pdf"),sep = "/")
        
        vis_mat = vis_mat[rownames(deconvolution_results),]
        
        if (type == "ductal"){
            
            scale_mat = data.frame(
                "MKI67" = deconvolution_results$MKI67,
                "Ductal" = deconvolution_results$ductal
            )
            lm.model <- lm(scale_mat$MKI67 ~ scale_mat$Ductal) # Fit linear model
            summary(lm.model)
            correlation = round(cor(scale_mat$MKI67, scale_mat$Ductal),2)
            cor.test(scale_mat$MKI67, scale_mat$Ductal)
            cor_ductal_MKI_67 = correlation
            cor_ductal_MKI_67_p_value = cor.test(scale_mat$MKI67, scale_mat$Ductal)$p.value
            
            cor_hisc_MKI_67 = 1
            cor_hisc_MKI_67_p_value = 1
            
        } else {
            
            scale_mat = data.frame(
                "MKI67" = deconvolution_results$MKI67,
                "HISC" = deconvolution_results$hisc
            )
            lm.model <- lm(scale_mat$MKI67 ~ scale_mat$HISC) # Fit linear model
            summary(lm.model)
            correlation = round(cor(scale_mat$MKI67, scale_mat$HISC),2)
            cor.test(scale_mat$MKI67, scale_mat$HISC)
            cor_hisc_MKI_67 = correlation
            cor_hisc_MKI_67_p_value = cor.test(scale_mat$MKI67, scale_mat$HISC)$p.value
            
            cor_ductal_MKI_67 = 1
            cor_ductal_MKI_67_p_value = 1
        }
    }
    
    ### stats
    
    if ( length(vis_mat$grading) > 0){ ### case grading available
        
        off_set = rnorm(nrow(deconvolution_results),mean=0.001,sd=0.001)
        
        # numerical grading
        grading_numeric = vis_mat$grading
        grading_numeric = as.integer(str_replace_all(grading_numeric,pattern ="G",""))
        
        # MKI-67 vs. grading numeric
        
        cor_MKI67_grading = cor(deconvolution_results$MKI67 + off_set,grading_numeric)
        cor_MKI67_grading_p_value = cor.test(deconvolution_results$MKI67,grading_numeric)$p.value
        
        # ductal vs. grading numeric
        
        cor_ductal_grading = cor(deconvolution_results$ductal + off_set,grading_numeric)
        cor_ductal_grading_p_value = cor.test(deconvolution_results$ductal  + off_set,grading_numeric)$p.value
        
        # hisc vs. grading numeric
        
        if ( type == "hisc"){
            
            off_set = rnorm(length(deconvolution_results$hisc),mean=0.001,sd=0.001)
            cor_hisc_grading = cor(deconvolution_results$hisc + off_set,grading_numeric)
            cor_hisc_grading_p_value = cor.test(deconvolution_results$hisc + off_set,grading_numeric)$p.value
        } else {cor_hisc_grading = cor_hisc_grading_p_value = 1.0}
        

    } else {
        
        cor_MKI67_grading = cor_ductal_grading = cor_hisc_grading = cor_ratio_grading = 0
        cor_MKI67_grading_p_value = cor_ductal_grading_p_value = cor_hisc_grading_p_value = cor_ratio_grading_p_value = 1
    }
    
    if( length(ki_index) != 0 ){
        
        off_set = rnorm( nrow(deconvolution_results),mean=0.001)
        # hisc vs. MKI-67
        
        if ( type == "hisc"){
            
            cor_hisc_MKI_67 = cor(deconvolution_results$hisc + off_set,deconvolution_results$MKI67)
            cor_hisc_MKI_67_p_value = cor.test(deconvolution_results$hisc + off_set,deconvolution_results$MKI67)$p.value
            
            cor_ductal_MKI_67 = 0
            cor_ductal_MKI_67_p_value = 1
            
        } else {
            
            cor_ductal_MKI_67 = cor(deconvolution_results$ductal + off_set,deconvolution_results$MKI67)
            cor_ductal_MKI_67_p_value = cor.test(deconvolution_results$ductal + off_set,deconvolution_results$MKI67)$p.value
            
            cor_hisc_MKI_67 = 0
            cor_hisc_MKI_67_p_value = 1
        }
    }
    
    ### anova
    
    # anova mki-67 versus grading
    
    if (
        (length(vis_mat$Grading) > 0) &
        (length(ki_index) != 0) &
        (length(unique(vis_mat$Grading)) > 1)
    ) {
        
        # anova 1
        
        anova_1 = aov(deconvolution_results$MKI67 ~ as.factor(as.character(vis_mat$Grading)) )
        anova_1_p_value = TukeyHSD(anova_1)$`as.factor(as.character(vis_mat$Grading))`
        G1_G2_index = which(rownames(anova_1_p_value) == "G2-G1")
        G1_G3_index = which(rownames(anova_1_p_value) == "G3-G1")
        G2_G3_index = which(rownames(anova_1_p_value) == "G3-G2")
        if(length(G1_G2_index) != 0){
            mki_67_g1_g2 = anova_1_p_value[G1_G2_index,4]
        } else {mki_67_g1_g2 = 1}
        if(length(G1_G3_index) != 0){
            mki_67_g1_g3 = anova_1_p_value[G1_G3_index,4]
        } else {mki_67_g1_g2 = 1}
        if(length(G2_G3_index) != 0){
            mki_67_g2_g3 = anova_1_p_value[G2_G3_index,4]
        } else {mki_67_g1_g2 = 1}

        #data_mat = reshape2::melt(deconvolution_results)
        data_mat = deconvolution_results[,c("Grading","MKI67")]
        data_mat$Sample = rownames(deconvolution_results)
        data_mat$Sample = factor(data_mat$Sample, levels = data_mat$Sample[order(data_mat$MKI67)] )
        data_mat$MKI67 = data_mat$MKI67
        
        color_vec = data_mat$Grading
        color_vec[color_vec == "G1"] = "green"
        color_vec[color_vec == "G2"] = "yellow"
        color_vec[color_vec == "G3"] = "red"
        color_vec = color_vec[order(data_mat$MKI67)]
        g = ggplot(data_mat,aes( x = Grading, y = MKI67, fill = Grading )) + geom_boxplot( )
        
        mki67_grading_path = paste("~/Deko_Projekt/Results/Images",algorithm, sep ="/")
        mki67_grading_path = paste(mki67_grading_path,"Grading", sep ="/")
        if (! dir.exists(mki67_grading_path)) dir.create(mki67_grading_path)
        mki67_grading_path = paste(mki67_grading_path,name_datatype, sep ="/")
        if (! dir.exists(mki67_grading_path)) dir.create(mki67_grading_path)
        mki67_grading_path = paste(mki67_grading_path,name_query_data, sep ="/")
        if (! dir.exists(mki67_grading_path)) dir.create(mki67_grading_path)
        mki67_grading_path = paste(mki67_grading_path,paste0(paste(name_training_data,"mki67",sep="_"),".pdf"),sep = "/")
        pdf(mki67_grading_path,onefile = FALSE)#,width="1024px",height="768px")
            plot(g)
        dev.off()
        
    } else {
        
        mki_67_g1_g2 = mki_67_g1_g3 = mki_67_g2_g3 = 1
    }

    # anova ductal versus grading
    
    if (
        (length(vis_mat$Grading) > 0) & 
        (length(unique(vis_mat$Grading)) > 1)
    ) {
        
        # anova 2 ductal
        
        off_set = rnorm(nrow(deconvolution_results),mean = 0.0001, sd = 0.0001)
        
        anova_2 = aov(deconvolution_results$ductal + off_set ~ as.factor(as.character(vis_mat$Grading)) )
        anova_2_p_value = TukeyHSD(anova_2)$`as.factor(as.character(vis_mat$Grading))`
        G1_G2_index = which(rownames(anova_2_p_value) == "G2-G1")
        G1_G3_index = which(rownames(anova_2_p_value) == "G3-G1")
        G2_G3_index = which(rownames(anova_2_p_value) == "G3-G2")
        if(length(G1_G2_index) != 0){
            ductal_g1_g2 = anova_2_p_value[G1_G2_index,4]
        } else {ductal_g1_g2 = 1}
        if(length(G1_G3_index) != 0){
            ductal_g1_g3 = anova_2_p_value[G1_G3_index,4]
        } else {ductal_g1_g3 = 1}
        if(length(G2_G3_index) != 0){
            ductal_g2_g3 = anova_2_p_value[G2_G3_index,4]
        } else {ductal_g2_g3 = 1}
        
        data_mat = reshape2::melt(deconvolution_results )
        data_mat = deconvolution_results[,c("Grading","ductal")]
        data_mat$Sample = rownames(data_mat)
        data_mat$Sample = factor(data_mat$Sample, levels = data_mat$Sample[order(data_mat$ductal)] )
        data_mat$ductal = data_mat$ductal
        
        color_vec = data_mat$Grading
        color_vec[color_vec == "G1"] = "green"
        color_vec[color_vec == "G2"] = "yellow"
        color_vec[color_vec == "G3"] = "red"
        color_vec = color_vec[order(data_mat$ductal)]
        g=ggplot(data_mat,aes( x = Grading, y = ductal, fill = Grading )) + geom_boxplot( )
        
        ductal_grading_path = paste("~/Deko_Projekt/Results/Images",algorithm, sep ="/")
        ductal_grading_path = paste(ductal_grading_path,"Grading", sep ="/")
        if (! dir.exists(ductal_grading_path)) dir.create(ductal_grading_path)
        ductal_grading_path = paste(ductal_grading_path,name_datatype, sep ="/")
        if (! dir.exists(ductal_grading_path)) dir.create(ductal_grading_path)
        ductal_grading_path = paste(ductal_grading_path,name_query_data, sep ="/")
        if (! dir.exists(ductal_grading_path)) dir.create(ductal_grading_path)
        ductal_grading_path = paste(ductal_grading_path,paste0(paste(name_training_data,"ductal",sep="_"),".pdf"),sep = "/")
        pdf(ductal_grading_path,onefile = FALSE)#,width="1024px",height="768px")
            plot(g)
        dev.off()
    } else { ductal_g1_g2 = ductal_g1_g3 = ductal_g2_g3 = 1 }
    
    # anova hisc versus grading
    
    if (
        (length(deconvolution_results$hisc) > 0) &
        (length(vis_mat$Grading) > 0) &
        (length(unique(vis_mat$Grading)) > 1)
    ){
        
        # anova 3
        
        anova_3 = aov(deconvolution_results$hisc ~ as.factor(as.character(vis_mat$Grading)) )
        anova_3_p_value = TukeyHSD(anova_3)$`as.factor(as.character(vis_mat$Grading))`
        G1_G2_index = which(rownames(anova_3_p_value) == "G2-G1")
        G1_G3_index = which(rownames(anova_3_p_value) == "G3-G1")
        G2_G3_index = which(rownames(anova_3_p_value) == "G3-G2")
        if(length(G1_G2_index) != 0){
            hisc_g1_g2 = anova_3_p_value[G1_G2_index,4]
        } else {hisc_g1_g2 = 1}
        if(length(G1_G3_index) != 0){
            hisc_g1_g3 = anova_3_p_value[G1_G3_index,4]
        } else {hisc_g1_g3 = 1}
        if(length(G2_G3_index) != 0){
            hisc_g2_g3 = anova_3_p_value[G2_G3_index,4]
        } else {hisc_g2_g3 = 1}

        data_mat = reshape2::melt(deconvolution_results )
        data_mat = deconvolution_results[,c("Grading","hisc")]
        data_mat$Sample = rownames(data_mat)
        data_mat$Sample = factor(data_mat$Sample, levels = data_mat$Sample[order(data_mat$hisc)] )
        data_mat$hisc = data_mat$hisc
        
        color_vec = data_mat$Grading
        color_vec[color_vec == "G1"] = "green"
        color_vec[color_vec == "G2"] = "yellow"
        color_vec[color_vec == "G3"] = "red"
        color_vec = color_vec[order(data_mat$hisc)]
        g=ggplot(data_mat,aes( x = Grading, y = hisc, fill = Grading )) + geom_boxplot( )
        
        hisc_grading_path = paste("~/Deko_Projekt/Results/Images",algorithm, sep ="/")
        hisc_grading_path = paste(hisc_grading_path,"Grading", sep ="/")
        if (! dir.exists(hisc_grading_path)) dir.create(hisc_grading_path)
        hisc_grading_path = paste(hisc_grading_path,name_datatype, sep ="/")
        if (! dir.exists(hisc_grading_path)) dir.create(hisc_grading_path)
        hisc_grading_path = paste(hisc_grading_path,name_query_data, sep ="/")
        if (! dir.exists(hisc_grading_path)) dir.create(hisc_grading_path)
        hisc_grading_path = paste(hisc_grading_path,paste0(paste(name_training_data,"hisc",sep="_"),".pdf"),sep = "/")
        pdf(hisc_grading_path,onefile = FALSE)#,width="1024px",height="768px")
            plot(g)
        dev.off()
        
    } else {
        
        hisc_g1_g2 = hisc_g1_g3 = hisc_g2_g3 = 1.0
    }
    
    ### output
    
    results_vec = c(
        dataset_query,
        dataset_training,
        cor_ductal_grading,
        cor_ductal_grading_p_value,
        cor_hisc_grading,
        cor_hisc_grading_p_value,
        cor_ductal_MKI_67,
        cor_ductal_MKI_67_p_value,
        cor_hisc_MKI_67,
        cor_hisc_MKI_67_p_value,
        ductal_g1_g2,
        ductal_g1_g3,
        ductal_g2_g3,
        hisc_g1_g2,
        hisc_g1_g3,
        hisc_g2_g3,
        surv_mki67,
        surv_ductal,
        surv_hisc,
        type
    )
    
    if (! file.exists(path_benchmark_files)){
        
        benchmark_results_t = matrix(results_vec,ncol=length(results_vec))
        colnames(benchmark_results_t) = c(
            "dataset_query",
            "dataset_training_label",
            "cor_ductal_grading",
            "cor_ductal_grading_p_value",
            "cor_hisc_grading",
            "cor_hisc_grading_p_value",
            "cor_ductal_MKI_67",
            "cor_ductal_MKI_67_p_value",
            "cor_hisc_MKI_67",
            "cor_hisc_MKI_67_p_value",
            "ductal_g1_g2",
            "ductal_g1_g3",
            "ductal_g2_ g3",
            "hisc_g1_g2",
            "hisc_g1_g3",
            "hisc_g2_g3",
            "surv_mki67",
            "surv_ductal",
            "surv_hisc",
            "type"
        )
    } else {
        
        benchmark_results_t = read.table(path_benchmark_files,sep="\t",stringsAsFactors = F,header = T)
        benchmark_results_t = rbind(benchmark_results_t,results_vec)    
    }

    # output 
    
    write.table(benchmark_results_t, path_benchmark_files,sep="\t",quote=F,row.names= F, col.names = T)
    return(benchmark_results_t)
}
    