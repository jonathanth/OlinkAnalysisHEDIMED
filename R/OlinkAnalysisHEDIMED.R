
#   __    __  .___________. __   __       __  .___________. __   _______      _______.
#  |  |  |  | |           ||  | |  |     |  | |           ||  | |   ____|    /       |
#  |  |  |  | `---|  |----`|  | |  |     |  | `---|  |----`|  | |  |__      |   (----`
#  |  |  |  |     |  |     |  | |  |     |  |     |  |     |  | |   __|      \   \
#  |  `--'  |     |  |     |  | |  `----.|  |     |  |     |  | |  |____ .----)   |
#   \______/      |__|     |__| |_______||__|     |__|     |__| |_______||_______/
#

show_outcome_distribution <- function(olink_data, outcome, annotation = NULL){
  outdat <- olink_data %>% group_by(sampleid) %>% slice(1) %>% ungroup %>% count({{outcome}})
  if(!is.null(annotation)){
    outdat <- dplyr::bind_cols(outdat, annotation = annotation)
  }
}

make_pca <- function(npx_data, getWide = FALSE){
  wide <- npx_data %>% dplyr::select(assay, npx, sampleid) %>%
    tidyr::pivot_wider(names_from = assay, values_from = npx) %>%
    dplyr::mutate(across(-sampleid, ~ifelse(is.na(.x), mean(.x, na.rm = T), .x)))
  if(getWide)
    return(wide)
  prcomp(wide %>% dplyr::select(-sampleid))
}

get_calibrated_residuals <- function(formula, data){
  model <- lm(formula, data)
  resids <- augment(model, newdata = data) %>% pull(.resid)
  intercept <- tidy(model) %>% slice(1) %>% pull(estimate)
  intercept + resids
}
#   _______   __   _______  _______     _______ ___   ___ .______
#  |       \ |  | |   ____||   ____|   |   ____|\  \ /  / |   _  \
#  |  .--.  ||  | |  |__   |  |__      |  |__    \  V  /  |  |_)  |
#  |  |  |  ||  | |   __|  |   __|     |   __|    >   <   |   ___/
#  |  '--'  ||  | |  |     |  |        |  |____  /  .  \  |  |
#  |_______/ |__| |__|     |__|        |_______|/__/ \__\ | _|
#

olink_lm <- function(olink_data, outcome, covariates = "1", annotation = ""){
  require(broom)
  res <- olink_data %>%
    filter(!is.na(olinkid)) %>%
    group_by(olinkid, uniprot, assay) %>%
    do(
      lm(paste0("npx ~ ", outcome, " + ", paste0(covariates, collapse = " + ")), data = .) %>% tidy(conf.int = T) %>%
        filter(grepl(outcome, term))
    ) %>%
    ungroup %>%
    dplyr::mutate(outcome = outcome, annotation = annotation, q.value = p.adjust(p.value, "fdr")) %>%
    dplyr::select(annotation, outcome, everything()) %>%
    arrange(p.value)
  res
}

olink_lmer <- function(olink_data, outcome, covariates = "1", annotation = ""){
  require(broom)
  require(lme4)
  require(lmerTest)
  require(broom.mixed)
  res <- olink_data %>%
    filter(!is.na(olinkid)) %>%
    group_by(olinkid, uniprot, assay) %>%
    do(
      lmer(paste0("npx ~ ", outcome, " + ", paste0(covariates, collapse = " + ")), data = .) %>% tidy(conf.int = T) %>%
        filter(grepl(outcome, term))
    ) %>%
    ungroup %>%
    dplyr::mutate(outcome = outcome, annotation = annotation, q.value = p.adjust(p.value, "fdr")) %>%
    dplyr::select(annotation, outcome, everything()) %>%
    arrange(p.value)
  res
}

#  ____    ____  ______    __        ______      ___      .__   __.   ______
#  \   \  /   / /  __  \  |  |      /      |    /   \     |  \ |  |  /  __  \
#   \   \/   / |  |  |  | |  |     |  ,----'   /  ^  \    |   \|  | |  |  |  |
#    \      /  |  |  |  | |  |     |  |       /  /_\  \   |  . `  | |  |  |  |
#     \    /   |  `--'  | |  `----.|  `----. /  _____  \  |  |\   | |  `--'  |
#      \__/     \______/  |_______| \______|/__/     \__\ |__| \__|  \______/
#

olink_volcano_jt <- function(olink_res){

  reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    scales::trans_new(paste0("reverselog-", format(base)), trans, inv,
                      scales::log_breaks(base = base),
                      domain = c(1e-100, Inf))
  }

  plottitle <- olink_res %>%
    slice(1) %>%
    dplyr::mutate(title = paste0(annotation, ": ", outcome)) %>%
    pull(title)
  olink_res %>%
    dplyr::mutate(siglevel = case_when(q.value < 0.05 ~ "FDR",
                                p.value < 0.05 ~ "Nominal",
                                TRUE ~ "N.S.") %>%
             factor(levels = c("N.S.", "Nominal", "FDR"))) %>%
    ggplot(aes(estimate, p.value)) +
    geom_point(aes(color = siglevel), size = 2) +
    ggrepel::geom_text_repel(aes(label = assay), data = . %>% filter(p.value < 0.05)) +
    geom_hline(yintercept = 0.05) +
    xlab("Estimate (NPX units)") +
    scale_y_continuous(trans = reverselog_trans(10), name = "P-value") +
    scale_color_manual(name = "Significance level", values = c("grey20", "blue", "green")) +
    ggtitle(plottitle) +
    theme_bw() +
    theme(legend.position = c(.02,.02), legend.justification = c(0,0),
          legend.box.background = element_rect(color = "black", fill = "white", linewidth = 1))
}

#    ______   ______   .______      .______           _______.
#   /      | /  __  \  |   _  \     |   _  \         /       |
#  |  ,----'|  |  |  | |  |_)  |    |  |_)  |       |   (----`
#  |  |     |  |  |  | |      /     |      /         \   \
#  |  `----.|  `--'  | |  |\  \----.|  |\  \----..----)   |
#   \______| \______/  | _| `._____|| _| `._____||_______/
#

olink_lm_corrs <- function(reslist, method = "spearman"){
  require(broom)
  resnames <- sapply(reslist, pull, annotation) %>% unlist %>% unique %>% as.character
  names(reslist) <- resnames
  corrlist <- list()
  corrlist_clean <- list()
  for(i in resnames){
    for(j in resnames){
      testdf <- full_join(reslist[[i]] %>% dplyr::select(annotation, assay, estimate),
                          reslist[[j]] %>% dplyr::select(annotation, assay, estimate),
                          by = "assay")
      corrlist[[length(corrlist) + 1]] <- cor.test(testdf$estimate.x, testdf$estimate.y, method = method)
      corrlist_clean [[length(corrlist) + 1]] <- data.frame(x = i, y = j, cor.test(testdf$estimate.x, testdf$estimate.y, method = method) %>% tidy)
    }
  }
  # corrlist
  bind_rows(corrlist_clean)
}


#     _____ ____  _____  _____    _    _ ______       _______ __  __          _____
#    / ____/ __ \|  __ \|  __ \  | |  | |  ____|   /\|__   __|  \/  |   /\   |  __ \
#   | |   | |  | | |__) | |__) | | |__| | |__     /  \  | |  | \  / |  /  \  | |__) |
#   | |   | |  | |  _  /|  _  /  |  __  |  __|   / /\ \ | |  | |\/| | / /\ \ |  ___/
#   | |___| |__| | | \ \| | \ \  | |  | | |____ / ____ \| |  | |  | |/ ____ \| |
#    \_____\____/|_|  \_\_|  \_\ |_|  |_|______/_/    \_\_|  |_|  |_/_/    \_\_|
#
#
olink_corrs_heatmap <- function(corrdata, xtextangle = 30){
  require(copiome)
  my_dist <- corrdata %>%
    dplyr::mutate(dist = 1 - abs(estimate)) %>%
    ungroup %>%
    tidyr::pivot_wider(id_cols = x, names_from = y, values_from = dist) %>%
    tibble::column_to_rownames("x") %>%
    as("matrix") %>%
    as.dist
  my_hc <- hclust(my_dist)
  corrdata %>%
    dplyr::mutate(x = factor(x, levels = my_hc$labels[my_hc$order]),
           y = factor(y, levels = my_hc$labels[my_hc$order])) %>%
    # dplyr::mutate(across(x:y, factor)) %>%
    ggplot(aes(x, y)) +
    geom_tile(aes(fill = estimate)) +
    geom_text(aes(label = paste0(roundex(estimate, 2), ifelse(p.value < 0.05, "*", ""))), data = . %>% filter(as.numeric(x) < as.numeric(y))) +
    scale_fill_gradientn(name = corrdata$method[1], colors = heatmap_palette, limits = c(-1, 1)) +
    theme_bw() +
    xlab(NULL) +
    ylab(NULL) +
    theme(axis.text.x = element_text(angle = xtextangle, hjust = 1)) +
    theme(legend.position = "bottom") +
    coord_equal()
}


#    ______ _ _         _____  _       _____
#   |  ____(_) |       |  __ \| |     / ____|
#   | |__   _| |_   ___| |__) | |    | (___
#   |  __| | | __| / __|  ___/| |     \___ \
#   | |    | | |_  \__ \ |    | |____ ____) |
#   |_|    |_|\__| |___/_|    |______|_____/
#
#

olink_fit_sPLS <- function(olink_data, outcome, seed = 123, fixX = c(), annotation = "sPLS analysis"){

  require(tidyverse)
  require(caret)
  require(mixOmics)
  require(mixOmicsCaret)
  require(doMC)
  require(ggforce)

  registerDoMC(cores = 9)

  repCV10 <- caret::trainControl(method = "repeatedcv",
                                 number = 10,
                                 repeats = 5,
                                 returnResamp = "all",
                                 savePredictions = "all",
                                 allowParallel = T,
                                 verboseIter = F)

  cytokine_list <- c("IL8", "VEGFA", "CD8A", "MCP-3", "GDNF", "CDCP1", "CD244",
                     "IL7", "OPG", "LAP TGF-beta-1", "uPA", "IL6", "IL-17C", "MCP-1",
                     "IL-17A", "CXCL11", "AXIN1", "TRAIL", "IL-20RA", "CXCL9", "CST5",
                     "IL-2RB", "IL-1 alpha", "OSM", "IL2", "CXCL1", "TSLP", "CCL4",
                     "CD6", "SCF", "IL18", "SLAMF1", "TGF-alpha", "MCP-4", "CCL11",
                     "TNFSF14", "FGF-23", "IL-10RA", "FGF-5", "MMP-1", "LIF-R", "FGF-21",
                     "CCL19", "IL-15RA", "IL-10RB", "IL-22 RA1", "IL-18R1", "PD-L1",
                     "Beta-NGF", "CXCL5", "TRANCE", "HGF", "IL-12B", "IL-24", "IL13",
                     "ARTN", "MMP-10", "IL10", "TNF", "CCL23", "CD5", "CCL3", "Flt3L",
                     "CXCL6", "CXCL10", "4E-BP1", "IL-20", "SIRT2", "CCL28", "DNER",
                     "EN-RAGE", "CD40", "IL33", "IFN-gamma", "FGF-19", "IL4", "LIF",
                     "NRTN", "MCP-2", "CASP-8", "CCL25", "CX3CL1", "TNFRSF9", "NT-3",
                     "TWEAK", "CCL20", "ST1A1", "STAMBP", "IL5", "ADA", "TNFB", "CSF-1")

  myformula <- paste0("as.numeric(", outcome, ") ~ ", paste(paste0("`", cytokine_list, "`"), collapse = " + ")) %>% as.formula()

  plsdat <- olink_data %>%
    filter(!is.na(.data[[outcome]]) & !is.na(sampleid) & !is.na(npx)) %>%
    tidyr::pivot_wider(id_cols = c("sampleid", outcome), names_from = assay, values_from = npx) %>%
    filter(complete.cases(.))

  set.seed(seed)
  pls_model_caret_train <- caret::train(myformula,
                                        data = plsdat,
                                        method = get_mixOmics_spls(),
                                        preProc = c("center", "scale"),
                                        tuneGrid = expand.grid(ncomp = 1,
                                                               keepX = c(1, 2, 4, 8, 16, 25, 40, 50, 60, 92),
                                                               keepY = 1),
                                        trControl = repCV10,
                                        fixX = fixX)

  bestTune <- pls_model_caret_train$pred %>%
    separate(Resample, c("Fold", "Rep")) %>%
    group_by(ncomp, keepX, Rep) %>%
    summarize(auc = pROC::auc(predictor = pred, obs, direction = "<", levels = c(0, 1)) %>% as.numeric) %>%
    group_by(ncomp, keepX) %>%
    filter(auc == median(auc)) %>%
    ungroup %>%
    filter(auc == max(auc)) %>%
    slice(1)

  CV_AUC <- bestTune$auc

  pls_model_caret_train$bestTune <- data.frame(ncomp = bestTune$ncomp, keepX = bestTune$keepX, keepY = 1)

  bestPreds <- pls_model_caret_train$pred %>%
    separate(Resample, c("Fold", "Rep")) %>%
    filter(ncomp == bestTune$ncomp & keepX == bestTune$keepX & Rep == bestTune$Rep)

  score <- bestPreds %>% arrange(rowIndex) %>%
    dplyr::bind_cols(plsdat) %>%
    dplyr::select(sampleid, pred, obs, one_of(outcome))

  aucplot <- pls_model_caret_train$pred %>%
    separate(Resample, c("Fold", "Rep")) %>%
    group_by(ncomp, keepX, Rep) %>%
    summarize(auc = pROC::auc(predictor = pred, obs, direction = "<", levels = c(0, 1)) %>% as.numeric) %>%
    ggplot(aes(x = factor(keepX), y = auc)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = factor(ncomp)), alpha = 0.2, show.legend = FALSE) +
    geom_mark_rect(aes(filter = ncomp == bestTune$ncomp & keepX == bestTune$keepX), linetype = "dashed") +
    geom_hline(yintercept = 0.5) +
    facet_grid(. ~ ncomp) +
    scale_color_brewer(palette = "Set1", name = NULL) +
    theme_bw() + theme(strip.background = element_blank()) +
    xlab("Number of variables") +
    ggtitle(annotation)

  scoreplot <- bestPreds %>%
    ggplot(aes(factor(obs), pred)) +
    geom_boxplot() +
    ggpubr::stat_compare_means(label.x.npc = 0.5, hjust = 0.5) +
    ylab("Score") +
    xlab(outcome) +
    theme_bw()

  loadings <- get_loadings(pls_model_caret_train, "CV", remove_empty = F)

  finalModel <- mixOmics::spls(plsdat %>% dplyr::select(-sampleid, -one_of(outcome)),
                               plsdat %>% pull(outcome) %>% as.numeric,
                               keepX = rep(bestTune$keepX, bestTune$ncomp),
                               ncomp = bestTune$ncomp,
                               scale = T)

  mypred <- predict(finalModel,
                    plsdat %>% dplyr::select(-sampleid, -one_of(outcome)))$predict

  train_AUC <- pROC::auc(predictor = mypred[,1,bestTune$ncomp],
                         plsdat %>% pull(outcome) %>% as.numeric,
                         direction = "<", levels = c(0, 1)) %>% as.numeric

  glm_model <- glm(obs ~ scale(pred), data = bestPreds, family = binomial)

  glm_result <- broom::tidy(glm_model, exp = T, conf.int = T)

  list(annotation = annotation,
       CV_AUC = CV_AUC,
       train_AUC = train_AUC,
       bestTune = bestTune,
       aucplot = aucplot,
       scoreplot = scoreplot,
       bestPreds = bestPreds,
       loadings = loadings,
       finalModel = finalModel,
       pls_model_caret_train = pls_model_caret_train,
       glm_model = glm_model,
       glm_result = glm_result,
       olink_data = olink_data,
       plsdat = plsdat,
       outcome = outcome,
       annotation = annotation,
       seed = seed)
}

make_cross_pls_list <- function(pls_objects = list()){
  n_models <- length(pls_objects)
  cross_results_list <- list()
  model_number <- 0
  for(i in 1:n_models){  # model i
    for(j in 1:n_models){   # dataset and predictions j
      model_number <- model_number + 1
      message(i, ", ", j)
      prediction_obj <- predict(pls_objects[[i]]$finalModel,
                                pls_objects[[j]]$plsdat %>%
                                  dplyr::select(one_of(colnames(pls_objects[[i]]$finalModel$input.X))))
      preds <- dplyr::bind_cols(pls_objects[[j]]$plsdat,
                         cross_preds = prediction_obj$predict[,1,pls_objects[[i]]$bestTune$ncomp] %>% scale %>% as.numeric)
      glm_model <- glm(as.formula(paste0(pls_objects[[j]]$outcome, " ~ cross_preds")), data = preds, family = binomial)
      glm_result <- glm_model %>% tidy(exp = T, conf.int = T)
      crossAUC <- pROC::auc(predictor = preds$cross_preds,
                            preds[,pls_objects[[j]]$outcome] %>% unlist %>% as.numeric,
                            direction = "<", levels = c(0, 1)) %>% as.numeric
      OR <- glm_result %>% slice(-1) %>% pull(estimate)
      p.value <- glm_result %>% slice(-1) %>% pull(p.value)
      cross_results_list[[model_number]] <- list(
        model_annotation = pls_objects[[i]]$annotation,
        prediction_annotation = pls_objects[[j]]$annotation,
        preds = preds,
        glm_model = glm_model,
        glm_result = glm_result,
        crossAUC = crossAUC,
        OR = OR,
        p.value = p.value
      )
    }
  }
  cross_results_list
}

make_cross_pls_list_plotdat <- function(cross_pls_list = list()){
  data.frame(model_annotation = sapply(cross_pls_list, function(x) x[["model_annotation"]]) %>% fct_inorder,
             prediction_annotation = sapply(cross_pls_list, function(x) x[["prediction_annotation"]]) %>% fct_inorder,
             crossAUC = sapply(cross_pls_list, function(x) x[["crossAUC"]]),
             OR = sapply(cross_pls_list, function(x) x[["OR"]]),
             p.value = sapply(cross_pls_list, function(x) x[["p.value"]]),
             conf.low = sapply(cross_pls_list, function(x) x[["glm_result"]][2,"conf.low"] %>% unlist %>% unname),
             conf.high = sapply(cross_pls_list, function(x) x[["glm_result"]][2,"conf.high"] %>% unlist %>% unname)

  )
}

plot_results_of_pls_models <- function(pls_model_list = list()){
  pls_model_list %>% lapply(function(x) x$glm_result) %>%
    setNames(pls_model_list %>% lapply(function(x) x$annotation)) %>%
    bind_rows(.id = "Annotation") %>%
    group_by(Annotation) %>%
    slice(-1) %>%
    ggplot(aes(estimate, Annotation)) +
    geom_vline(xintercept = 1) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
    geom_point(size = 2) +
    scale_x_log10(name = "Odds Ratio per SD in fingerprint score") +
    theme_bw() +
    ggtitle("Summary of cross-validated fingerprint scores")
}

plot_cross_auc_heatmap <- function(cross_results = data.frame()){
  # Make cross-AUC heatmap
  cross_results %>%
    make_cross_pls_list_plotdat %>%
    ggplot(aes(model_annotation, prediction_annotation)) +
    geom_tile(aes(fill = crossAUC)) +
    geom_label(aes(label = paste0("AUC=", roundex(crossAUC, 2), "\np=", roundex(p.value, 3))), data = . %>% filter(crossAUC > 0.5 & p.value < 0.05)) +
    scale_fill_gradient2(midpoint = 0.5) +
    coord_equal() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

}
