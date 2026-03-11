
library(terra)
library(sf)
library(tidyverse)
library(exactextractr)
library(xlsx)
library(openxlsx)
library(xlsx)
library(signal)
library(zoo)
library(stringr)
library(lubridate)
library(gslnls)
# library(tidymodels)
library(minpack.lm)
library(skimr)
library(minpack.lm)
library(dplyr)
library(tidyr)


######## 
########  UTILS GEOESPTIAL PRODUCTS
######## 


#### MASK

calculate_cc_pixel <- function(r, g, b, p1=0.95, p2=0.95, p3=20, p4=NULL, p5=NULL) {
  cond1 <- (r / g < p1) & (b / g < p2) & ((2 * g - r - b) > p3)
  cond2 <- FALSE  # <--- CHANGED FROM 0 TO FALSE
  if (!is.null(p4) && !is.null(p5)) {
    cond2 <- ((r + g + b) < p4) & ((g - b) > p5)
  }
  cc_map <- (cond1 | cond2) 
  cc_map[is.na(cc_map)] <- 0
  return(cc_map)
}


calculate_canopy_cover <- function(cc_mask, polygons, id_col = "id") {
  stats <- exactextractr::exact_extract(cc_mask, polygons, function(values, coverage_fraction) {
    sum(values == 1, na.rm = TRUE) / length(values)
  })
  
  results <- data.frame(
    COD = polygons[[id_col]],
    CanopyCover_Pct = stats * 100  
  )
  
  return(results)
}


extract_cv <- function(dsm, dtm, polygons, id_col = "id", cc_mask = NULL) {
  
  chm <- dsm - dtm
  chm[chm < 0.05] <- 0
  chm <- terra::mask(chm, cc_mask, maskvalues = 0, updatevalue = 0)
  pixel_area <- terra::res(chm)[1] * terra::res(chm)[2]
  vol_raster <- chm * pixel_area
  vol_sums <- exactextractr::exact_extract(vol_raster, polygons, 'sum')
  results <- data.frame(
    COD = polygons[[id_col]],
    Volume_m3 = vol_sums
  )
  return(results)
}

extract_exg_mean <- function(red, green, blue, cc_mask, polygons, id_col = "id") {
  
  total <- red + green + blue
  r_norm <- red / total
  g_norm <- green / total
  b_norm <- blue / total
  
  exg <- 2 * g_norm - r_norm - b_norm
  
  exg_masked <- terra::mask(exg, cc_mask, maskvalues = 0)
  
  mean_vals <- exactextractr::exact_extract(exg_masked, polygons, 'mean')
  
  results <- data.frame(
    COD = polygons[[id_col]],
    ExG = mean_vals
  )
  
  return(results)
}

extract_ndvi_mean <- function(red, nir, cc_mask, polygons, id_col = "id") {
  
  ndvi <- (nir - red) / (nir + red)
  
  ndvi_masked <- terra::mask(ndvi, cc_mask, maskvalues = 0)
  
  mean_vals <- exactextractr::exact_extract(ndvi_masked, polygons, 'mean')
  
  results <- data.frame(
    COD = polygons[[id_col]],
    NDVI = mean_vals
  )
  
  return(results)
}


########### WRAPPER

extract_cc_pct <- function(cc_mask, polygons, id_col = "id") {
  cc_vals <- exactextractr::exact_extract(cc_mask, polygons, function(values, coverage_fraction) {
    sum(values == 1, na.rm = TRUE) / length(values)
  }, progress = FALSE)
  
  results <- data.frame(
    COD = polygons[[id_col]],
    CC_pct = cc_vals * 100
  )
  return(results)
}


###### ZONAL STATS

extract_htp_zonal <- function(multi_path   = "PROCESING/MULTI/", 
                              dsm_path     = "PROCESING/DSM/", 
                              border_path  = "PROCESING/BORDERS/",
                              dtm_file     = "PROCESING/DSM/dtm_dap_62_r.tif",
                              save_rasters = FALSE,
                              output_path  = "PROCESING/STACKS/",
                              group_var    = "COD",
                              daps) {
  
  message("--- Starting Full Remote Sensing Pipeline ---")
  
  dtm_ref <- terra::rast(dtm_file)
  
  if (save_rasters) {
    products <- c("CC", "CV", "ExG", "NDVI")
    for(prod in products) {
      dir.create(file.path(output_path, prod), showWarnings = FALSE, recursive = TRUE)
    }
  }
  
  all_results <- purrr::map_dfr(daps, function(dap) {
    
    message(paste("--- Processing DAP:", dap, "---"))
    
    # Strict uppercase "DAP_" pattern search
    m_file <- list.files(multi_path,  pattern = paste0("DAP_", dap), full.names = TRUE)[1]
    s_file <- list.files(dsm_path,    pattern = paste0("DAP_", dap), full.names = TRUE)[1]
    b_file <- list.files(border_path, pattern = paste0("DAP_", dap), full.names = TRUE)[1]
    
    if (any(is.na(c(m_file, s_file, b_file)))) {
      warning(paste("Missing files for DAP", dap, "- Check if filenames use 'DAP_' (Uppercase)"))
      return(NULL)
    }
    
    message(paste("   > Using Multi:", basename(m_file)))
    message(paste("   > Using DSM:  ", basename(s_file)))
    message(paste("   > Using Border:", basename(b_file)))
    
    r_list   <- terra::rast(m_file)
    dsm      <- terra::rast(s_file)
    polygons <- sf::st_read(b_file, quiet = TRUE)
    
    r <- r_list[["Red"]]; g <- r_list[["Green"]]; b <- r_list[["Blue"]]; nir <- r_list[["NIR"]]
    
    cc_mask  <- calculate_cc_pixel(r * 255, g * 255, b * 255, p1 = 1.05, p2 = 1.00, p3 = 5)
    
    chm <- dsm - dtm_ref
    chm[chm < 0.05] <- 0
    cv_raster <- terra::mask(chm, cc_mask, maskvalues = 0, updatevalue = NA)
    
    total <- r + g + b
    exg   <- 2 * (g/total) - (r/total) - (b/total)
    exg_masked <- terra::mask(exg, cc_mask, maskvalues = 0, updatevalue = NA)
    
    ndvi  <- (nir - r) / (nir + r)
    ndvi_masked <- terra::mask(ndvi, cc_mask, maskvalues = 0, updatevalue = NA)
    
    if (save_rasters) {
      terra::writeRaster(cc_mask,    file.path(output_path, "CC",   paste0("CC_DAP_", dap, ".tif")),   overwrite = TRUE)
      terra::writeRaster(cv_raster,  file.path(output_path, "CV",   paste0("CV_DAP_", dap, ".tif")),   overwrite = TRUE)
      terra::writeRaster(exg_masked, file.path(output_path, "ExG",  paste0("ExG_DAP_", dap, ".tif")),  overwrite = TRUE)
      terra::writeRaster(ndvi_masked,file.path(output_path, "NDVI", paste0("NDVI_DAP_", dap, ".tif")), overwrite = TRUE)
    }
    
    cv_res   <- extract_cv(dsm, dtm_ref, polygons, group_var, cc_mask)
    exg_res  <- extract_exg_mean(r, g, b, cc_mask, polygons, group_var)
    ndvi_res <- extract_ndvi_mean(r, nir, cc_mask, polygons, group_var)
    cc_res   <- extract_cc_pct(cc_mask, polygons, group_var)
    
    dap_combined <- data.frame(
      id    = polygons$id,
      COD   = polygons$COD,
      BLOQ  = if("BLOQ" %in% colnames(polygons)) polygons$BLOQ else NA,
      CV    = cv_res$Volume_m3,
      ExG   = exg_res$ExG,
      NDVI  = ndvi_res$NDVI,
      CC    = cc_res$CC_pct,
      DAP   = as.numeric(dap)
    )
    
    return(dap_combined)
  })
  
  message("--- Pipeline Complete ---")
  return(all_results)
}

##### PLOT RASTERS


plot_htp <- function(output_path = "PROCESING/STACKS_07_03_26V2/",
                     multi_path  = "PROCESING/MULTI/",
                     target_cod  = "CQC-183",
                     target_bloq = "B1",
                     border_file = "PROCESING/BORDERS/BORDERS_DAP_62_bloq.gpkg",
                     daps        = c("62", "86", "93", "121", "128"),
                     plot_path = NULL) {
  
  library(terra)
  library(sf)
  library(ggplot2)
  library(tidyterra)
  library(dplyr)
  library(patchwork)
  
  boundaries <- st_read(border_file, quiet = TRUE)
  plot_poly  <- boundaries %>% dplyr::filter(COD == target_cod, BLOQ == target_bloq)
  
  plot_list <- list()
  row_names <- c("rgb", "cc", "cv", "ndvi", "exg")
  
  for (dap in daps) {
    rgb_path <- list.files(multi_path, pattern = paste0("DAP_", dap), full.names = TRUE)[1]
    rgb_raw  <- rast(rgb_path)[[c("Red", "Green", "Blue")]]
    rgb_crop <- crop(rgb_raw, ext(plot_poly)) %>% mask(plot_poly)
    rgb_bright <- stretch(rgb_crop, minv=0, maxv=255, minq=0.02, maxq=0.98)
    
    ndvi_r <- rast(file.path(output_path, "NDVI", paste0("NDVI_DAP_", dap, ".tif"))) %>% crop(ext(plot_poly)) %>% mask(plot_poly)
    cc_r   <- rast(file.path(output_path, "CC",   paste0("CC_DAP_", dap, ".tif")))   %>% crop(ext(plot_poly)) %>% mask(plot_poly)
    cv_r   <- rast(file.path(output_path, "CV",   paste0("CV_DAP_", dap, ".tif")))   %>% crop(ext(plot_poly)) %>% mask(plot_poly)
    exg_r  <- rast(file.path(output_path, "ExG",  paste0("ExG_DAP_", dap, ".tif")))  %>% crop(ext(plot_poly)) %>% mask(plot_poly)
    
    p_rgb  <- ggplot() + geom_spatraster_rgb(data = rgb_bright) + theme_void() + labs(title = paste("dap", dap))
    p_cc   <- ggplot() + geom_spatraster(data = cc_r)   + scale_fill_gradient(low="white", high="darkgreen", na.value=NA) + theme_void() + theme(legend.position="none")
    p_cv   <- ggplot() + geom_spatraster(data = cv_r)   + scale_fill_viridis_c(option = "magma", na.value=NA) + theme_void() + theme(legend.position="none")
    p_ndvi <- ggplot() + geom_spatraster(data = ndvi_r) + scale_fill_gradientn(colors = rev(terrain.colors(10)), na.value=NA) + theme_void() + theme(legend.position="none")
    p_exg  <- ggplot() + geom_spatraster(data = exg_r)  + scale_fill_gradient2(low="brown", mid="white", high="forestgreen", na.value=NA) + theme_void() + theme(legend.position="none")
    
    if (dap == daps[1]) {
      label_style <- theme(axis.title.y = element_text(angle=90, size=8, face="bold", vjust=1),
                           display.axis.titles = TRUE)
      
      
      p_rgb  <- p_rgb  + ylab("rgb")  + theme(axis.title.y = element_text(angle=90, size=14, face="bold", vjust=1))
      p_cc   <- p_cc   + ylab("cc")   + theme(axis.title.y = element_text(angle=90, size=14, face="bold", vjust=1))
      p_cv   <- p_cv   + ylab("cv")   + theme(axis.title.y = element_text(angle=90, size=14, face="bold", vjust=1))
      p_ndvi <- p_ndvi + ylab("ndvi") + theme(axis.title.y = element_text(angle=90, size=14, face="bold", vjust=1))
      p_exg  <- p_exg  + ylab("exg")  + theme(axis.title.y = element_text(angle=90, size=14, face="bold", vjust=1))
    }
    
    plot_list[[paste0(dap, "_1")]] <- p_rgb
    plot_list[[paste0(dap, "_2")]] <- p_cc
    plot_list[[paste0(dap, "_3")]] <- p_cv
    plot_list[[paste0(dap, "_4")]] <- p_ndvi
    plot_list[[paste0(dap, "_5")]] <- p_exg
  }
  
  final_plot <- wrap_plots(plot_list, ncol = length(daps), byrow = FALSE) + 
    plot_annotation(title = paste(target_cod, "-", target_bloq)) & 
    theme(plot.title = element_text(hjust = 0.5, size = 14)) #face = "bold"
  
  ggsave(file.path(plot_path, paste0("geoProducst_", target_bloq , ".png")), 
         plot = final_plot, width = 10, height = 12, limitsize = FALSE)
  
  return(final_plot)
}





######################
###################### UTILS FITING
######################

rmse <- function(observed, predicted) {
  sqrt(mean((observed - predicted)^2, na.rm = TRUE))
}

r2 <- function(observed, predicted) {
  1 - sum((observed - predicted)^2, na.rm = TRUE) / 
    sum((observed - mean(observed, na.rm = TRUE))^2, na.rm = TRUE)
}

extract_features_spline <- function(data, response, time_var) {
  data <- as.data.frame(data)
  
  data <- data[!is.na(data[[response]]) & 
                 !is.na(data[[time_var]]) & 
                 data[[response]] > 0, ]
  
  if (nrow(data) < 4) return(NULL)
  data <- data[order(data[[time_var]]), ]
  if (length(unique(data[[time_var]])) < 4) return(NULL)
  
  fit <- tryCatch({
    smooth.spline(x = data[[time_var]], y = data[[response]])
  }, error = function(e) return(NULL))
  
  if (is.null(fit)) return(NULL)
  
  obs_y <- data[[response]]
  pred_y_at_obs <- predict(fit, data[[time_var]])$y
  
  fit_rmse <- rmse(obs_y, pred_y_at_obs)
  fit_r2 <- r2(obs_y, pred_y_at_obs)
  
  t_seq <- seq(min(data[[time_var]]), max(data[[time_var]]), length.out = 1000)
  pred <- predict(fit, t_seq)
  y <- pred$y
  d1 <- predict(fit, t_seq, deriv = 1)$y
  dt <- t_seq[2] - t_seq[1]
  
  idx_max_y <- which.max(y)
  idx_max_d1 <- which.max(d1)
  
  max_val <- y[idx_max_y]
  dap_at_max <- t_seq[idx_max_y]
  max_growth_rate <- d1[idx_max_d1]
  dap_at_max_growth <- t_seq[idx_max_d1]
  
  half_max <- max_val / 2
  is_above_half <- y >= half_max
  duration_half <- if(any(is_above_half)) diff(range(t_seq[is_above_half])) else 0
  
  inc_mask <- t_seq <= dap_at_max
  dec_mask <- t_seq > dap_at_max
  
  metrics <- switch(
    response,
    "CC" = data.frame(F1=max_val, F2=max_growth_rate, F3=dap_at_max_growth, F4=duration_half),
    "CV" = data.frame(F5=max_val, F6=max_growth_rate, F7=dap_at_max_growth, F8=duration_half),
    "ExG" = data.frame(
      F9=max_val, F10=dap_at_max, 
      F11=max(d1[inc_mask], na.rm=TRUE), F12=y[idx_max_y], 
      F13=sum(inc_mask)*dt, F14=sum(y[inc_mask])*dt,
      F15=min(d1[dec_mask], na.rm=TRUE), F16=y[idx_max_y],
      F17=sum(dec_mask)*dt, F18=sum(y[dec_mask])*dt
    ),
    "NDVI" = data.frame(
      F19=max_val, F20=dap_at_max, 
      F21=max(d1[inc_mask], na.rm=TRUE), F22=min(d1[dec_mask], na.rm=TRUE)
    )
  )
  
  metrics[[paste0("R2_", response)]] <- fit_r2
  metrics[[paste0("RMSE_", response)]] <- fit_rmse
  
  return(list(fit = fit, t_seq = t_seq, y = y, metrics = metrics))
}

extract_htp_pheno <- function(df, 
                              group_var    = "COD", 
                              bloq_var     = "BLOQ", 
                              time_var     = "DAP", 
                              plot         = FALSE, 
                              plot_path    = "PROCESING/PLOTS/",
                              genotypes    = NULL) {
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  if (!is.null(bloq_var) && bloq_var %in% names(df)) {
    df <- df %>% dplyr::mutate(PLOT_ID = paste(.data[[group_var]], .data[[bloq_var]], sep = "_"))
    fit_var <- "PLOT_ID"
  } else {
    fit_var <- group_var
    bloq_var <- NULL 
  }
  
  vars_to_process <- intersect(c("CC", "CV", "ExG", "NDVI"), names(df))
  unique_fits <- unique(df[[fit_var]])
  
  all_index_results <- list()
  if (plot) dir.create(plot_path, showWarnings = FALSE, recursive = TRUE)
  
  for (v in vars_to_process) {
    group_results <- list()
    line_list <- list() 
    
    for (g in unique_fits) {
      sub_df <- df[df[[fit_var]] == g, c(fit_var, time_var, v)]
      res <- extract_features_spline(sub_df, v, time_var)
      
      if (!is.null(res)) {
        row_data <- res$metrics
        row_data[[fit_var]] <- g
        group_results[[g]] <- row_data
        
        if (plot) {
          orig_meta <- df %>% 
            dplyr::filter(!!dplyr::sym(fit_var) == g) %>% 
            dplyr::select(dplyr::any_of(c(group_var, bloq_var))) %>% 
            dplyr::distinct()
          
          line_list[[g]] <- data.frame(
            FIT_ID = g,
            time_seq = res$t_seq,
            pred_y = res$y,
            R2 = round(res$metrics[[paste0("R2_", v)]], 3),
            RMSE = round(res$metrics[[paste0("RMSE_", v)]], 4)
          ) %>% dplyr::bind_cols(orig_meta)
        }
      }
    }
    
    all_index_results[[v]] <- do.call(rbind, group_results)
    
    if (plot && length(line_list) > 0) {
      plot_data <- do.call(rbind, line_list)
      
      if (!is.null(genotypes)) {
        df_p <- df %>% dplyr::filter(!!dplyr::sym(group_var) %in% genotypes)
        ln_p <- plot_data %>% dplyr::filter(!!dplyr::sym(group_var) %in% genotypes)
        suffix <- "_subset"
      } else {
        df_p <- df
        ln_p <- plot_data
        suffix <- "_full"
      }
      
      if (nrow(ln_p) > 0) {
        label_data <- ln_p %>% 
          dplyr::group_by(across(dplyr::any_of(c(group_var, bloq_var)))) %>% 
          dplyr::summarize(R2 = dplyr::first(R2), RMSE = dplyr::first(RMSE), .groups = 'drop') %>%
          dplyr::mutate(label = paste0("R2: ", R2, "\nRMSE: ", RMSE))
        
        p <- ggplot() +
          geom_point(data = df_p, aes(x = .data[[time_var]], y = .data[[v]]), alpha = 0.5) +
          geom_line(data = ln_p, aes(x = time_seq, y = pred_y), color = "blue") +
          geom_text(data = label_data, aes(x = -Inf, y = Inf, label = label), 
                    hjust = -0.1, vjust = 1.1, size = 2.5, inherit.aes = FALSE) +
          labs(title = paste("Fit:", v), x = time_var, y = v)
        
        if (!is.null(bloq_var)) {
          p <- p + facet_grid(rows = vars(!!dplyr::sym(group_var)), cols = vars(!!dplyr::sym(bloq_var)), scales = "free_y")
        } else {
          p <- p + facet_wrap(as.formula(paste("~", group_var)), ncol = 5, scales = "free_y")
        }
        
        h_val <- length(unique(ln_p[[group_var]])) * 2 + 1
        ggsave(file.path(plot_path, paste0("Grid_Fit_", v, suffix, ".png")), 
               plot = p, width = 10, height = h_val, limitsize = FALSE)
      }
    }
  }
  
  master_df <- Reduce(function(x, y) merge(x, y, by = fit_var, all = TRUE), all_index_results)
  
  if (fit_var == "PLOT_ID") {
    master_df <- master_df %>% tidyr::separate(PLOT_ID, into = c(group_var, "BLOQ_TEMP"), sep = "_", remove = FALSE)
    colnames(master_df)[colnames(master_df) == "BLOQ_TEMP"] <- bloq_var
  }
  
  return(list(
    features = master_df %>% dplyr::select(dplyr::any_of(c(group_var, bloq_var)), starts_with("F")),
    quality  = master_df %>% dplyr::select(dplyr::any_of(c(group_var, bloq_var)), contains("R2_"), contains("RMSE_"))
  ))
}

# BOXPLOTS

boxplot <- function(df_features, 
                    htp_feature = "F5", 
                    genotypes = NULL, 
                    save_plot = FALSE,
                    plot_path = "PROCESING/PLOTS/BOXPLOTS/") {
  
  library(ggplot2)
  library(dplyr)
  
  df_plot <- if(!is.null(genotypes)) {
    df_features[df_features$COD %in% genotypes, ]
  } else {
    df_features
  }
  
  p <- ggplot(df_plot, aes(x = COD, y = .data[[htp_feature]])) +
    geom_boxplot(outlier.shape = NA, width = 0.6, fill = "white", color = "black", linewidth = 0.4) +
    geom_jitter(width = 0.2, size = 0.8, alpha = 0.6) +
    labs(x = "genotype", y = paste( htp_feature)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7),
      axis.title = element_text(size = 8),
      panel.grid.major = element_line(color = "white"),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "#EBEBEB", color = NA),
      plot.background = element_blank()
    )
  
  if (save_plot) {
    dir.create(plot_path, showWarnings = FALSE, recursive = TRUE)
    ggsave(file.path(plot_path, paste0("Small_Boxplot_", htp_feature, ".png")), 
           plot = p, 
           width = 4,      # Fixed Width in inches
           height = 2.5,   # Fixed Height in inches
           units = "in", 
           dpi = 300)
  }
  
  return(p)
}



######################
###################### UTILS REGRESION MODELING 
######################


##### CORRELATION

htp_correlations <- function(df, htp_features, trait = "PMSEM") {
  library(dplyr)
  library(purrr)
  combined_df <- df %>%
    inner_join(htp_features, by = c("COD", "BLOQ"))
  
  feat_names <- htp_features %>% 
    dplyr::select(starts_with("F")) %>% 
    names()
  
  cor_results <- map_dfr(feat_names, function(f) {
    
    if (!f %in% names(combined_df) || !is.numeric(combined_df[[f]])) return(NULL)
    test <- cor.test(combined_df[[f]], combined_df[[trait]], use = "pairwise.complete.obs")
    
    data.frame(
      Variable = f,
      Correlation = round(test$estimate, 3),
      P_Value = round(test$p.value, 4),
      Significance = case_when(
        test$p.value < 0.001 ~ "***",
        test$p.value < 0.01  ~ "**",
        test$p.value < 0.05  ~ "*",
        test$p.value < 0.1   ~ ".",
        TRUE                 ~ ""
      )
    )
  })
  
  return(cor_results %>% arrange(desc(abs(Correlation))))
}


#### REGRESION

ht_regression <- function(df, htp_df, trait = "PMSEM", na_threshold = 0.3) {
  library(dplyr)
  library(tidyr)
  library(purrr)
  
  htp_clean_cols <- htp_df %>%
    dplyr::select(where(~ mean(is.na(.)) <= na_threshold))
  
  combined_df <- dplyr::inner_join(df, htp_clean_cols, by = c("COD", "BLOQ"))
  
  feat_names <- htp_clean_cols %>% 
    dplyr::select(dplyr::starts_with("F"), where(is.numeric)) %>% 
    dplyr::select(-dplyr::any_of(c("COD", "BLOQ"))) %>% 
    names()
  
  data_clean <- combined_df %>%
    dplyr::select(dplyr::all_of(trait), dplyr::all_of(feat_names)) %>%
    dplyr::filter(dplyr::if_all(dplyr::everything(), ~ is.finite(.)))
  
  data_clean <- data_clean %>%
    dplyr::select(dplyr::all_of(trait), where(~ sd(., na.rm = TRUE) > 0))
  
  feat_names <- setdiff(names(data_clean), trait)
  
  if (nrow(data_clean) <= (length(feat_names) + 1)) {
    stop(paste("Overfitting risk: Only", nrow(data_clean), "rows for", length(feat_names), "features."))
  }
  
  data_clean[[trait]] <- log10(data_clean[[trait]] + 1e-6)
  
  data_scaled <- data_clean %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(feat_names), ~ as.numeric(scale(.))))
  
  full_f <- as.formula(paste(trait, "~ ."))
  null_f <- as.formula(paste(trait, "~ 1"))
  
  m_full <- lm(full_f, data = data_scaled)
  
  m_fwd <- tryCatch(
    step(lm(null_f, data = data_scaled), scope = list(lower = null_f, upper = m_full), direction = "forward", trace = 0),
    error = function(e) return(NULL)
  )
  
  m_back <- tryCatch(
    step(m_full, direction = "backward", trace = 0),
    error = function(e) return(NULL)
  )
  
  extract_stats <- function(m, name) {
    if(is.null(m)) return(data.frame(Approach = name, R2 = NA, Adj_R2 = NA, RSE = NA, Num_Features = NA))
    s <- summary(m)
    data.frame(
      Approach = name,
      R2 = round(s$r.squared, 4),
      Adj_R2 = round(s$adj.r.squared, 4),
      RSE = round(s$sigma, 4),
      Num_Features = length(coef(m)) - 1,
      stringsAsFactors = FALSE
    )
  }
  
  summary_tab <- dplyr::bind_rows(
    extract_stats(m_full, "Full Model"),
    extract_stats(m_fwd, "Forward Selection"),
    extract_stats(m_back, "Backward Elimination")
  )
  
  return(list(
    summary = summary_tab,
    models = list("Full" = m_full, "Forward" = m_fwd, "Backward" = m_back),
    data = data_scaled,
    trait = trait
  ))
}

### PLOT REGRESION MODEL


plot_comparison_grid <- function(analysis_output, plot_path) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  target <- analysis_output$trait
  mods <- analysis_output$models
  df <- analysis_output$data
  
  valid_names <- names(mods)[sapply(mods, function(x) !is.null(x))]
  
  all_plots_df <- lapply(valid_names, function(nm) {
    m <- mods[[nm]]
    s <- summary(m)
    data.frame(
      Obs = df[[target]], 
      Fit = predict(m),
      ModelName = nm,
      Stats = paste0("R² = ", round(s$r.squared, 3), "\nRSE = ", round(s$sigma, 3))
    )
  }) %>% bind_rows()
  
  final_plot <- ggplot(all_plots_df, aes(x = Fit, y = Obs)) +
    geom_point(size = 1.5) +
    geom_smooth(method = "lm", color = "blue", se = FALSE, linewidth = 0.8) +
    geom_label(aes(label = Stats), x = Inf, y = -Inf, 
               hjust = 1.05, vjust = -0.2, size = 3.5, 
               fill = "white", label.size = 0.2, label.r = unit(0, "lines")) +
    facet_wrap(~ModelName, scales = "free") +
    labs(x = "predicted.values", y = "observed.values") +
    theme_grey() +
    theme(
      aspect.ratio = 1,
      strip.text = element_text(size = 10),
      panel.grid.major = element_line(color = "white"),
      panel.grid.minor = element_line(color = "white"),
      axis.title = element_text(size = 12)
    )
  
  if (!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)
  
  ggsave(file.path(plot_path, "regresion_plots.png"), 
         plot = final_plot, width = 10, height = 5, dpi = 300)
  
  return(final_plot)
}

########## 
########## ANOVA
########## 

run_glmer_anova <- function(df, htp_df, response = "PMSEM", treatment = "COD", block = "BLOQ") {
  library(lme4)
  library(broom.mixed)
  library(MuMIn)
  library(DHARMa)
  library(emmeans)
  library(multcomp)
  library(dplyr)
  
  combined_df <- df %>% inner_join(htp_df, by = c("COD", "BLOQ"))
  
  clean_df <- combined_df %>%
    filter(is.finite(!!sym(response)), 
           !is.na(!!sym(treatment)), 
           !is.na(!!sym(block)))
  
  if (nrow(clean_df) < 5) {
    warning(paste("Not enough finite data for response:", response))
    return(NULL)
  }
  
  formula_str <- paste(response, "~", treatment, "+ (1|", block, ") - 1")
  
  m1 <- tryCatch({
    lmer(as.formula(formula_str), data = clean_df)
  }, error = function(e) {
    message(paste("Model failed for", response, ":", e$message))
    return(NULL)
  })
  
  if (is.null(m1)) return(NULL)
  
  r2 <- r.squaredGLMM(m1)[1,2]
  sim_resid <- simulateResiduals(m1)
  dispersion <- testDispersion(sim_resid, plot = FALSE)
  anova_res <- car::Anova(m1, type = "III")
  
  em <- emmeans(m1, as.formula(paste("pairwise ~", treatment)), adjust = "tukey")
  res_cld <- cld(em$emmeans, Letters = letters)
  
  list(
    "Fixed_Effects" = tidy(m1, effects = "fixed"),
    "Random_Effects" = tidy(m1, effects = "ran_pars"),
    "Fit" = glance(m1) %>% mutate(r2 = r2, dispersion = as.numeric(dispersion$statistic)),
    "Anova" = tidy(anova_res),
    "Means" = tidy(res_cld)
  )
}


get_model_info <- function(model, var) {
  r2_val <- MuMIn::r.squaredGLMM(model)[1, 2]
  sim_resid <- DHARMa::simulateResiduals(model)
  disp_test <- DHARMa::testDispersion(sim_resid, plot = FALSE)
  
  model_stats <- broom.mixed::glance(model) %>%
    dplyr::mutate(
      variable = var,
      r2 = r2_val,
      dispersion = as.numeric(disp_test$statistic)
    )
  
  list(
    "Fixed_Effects" = broom.mixed::tidy(model, effects = "fixed"),
    "Random_Effects" = broom.mixed::tidy(model, effects = "ran_pars"),
    "Fit" = model_stats
  )
}

####
#### PLOT OBS VS PRED
####


# get_predictions <- function(model, data, source){
#   
#   fitted_values <- fitted(model)
#   
#   plot_data <- data.frame(
#     Observed = data$PGRA,
#     Fitted = fitted_values,
#     VAR = data$VAR  # Add the VAR column
#   ) %>% dplyr::mutate(input_data = source)
#   
#   return(plot_data)
#   
# }
# 
# observed_vs_fitted_plot <- function(data) {
#   plot <- ggplot(data, aes(x = Observed, y = Fitted)) +
#     geom_point(color = "blue", fill = "blue", size = 3, shape = 21, alpha = 0.3) + # Blue points
#     geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black", size = 1) + # 1:1 line
#     geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "red", size = 0.5) + # Regression line
#     
#     stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., sep = "~~~")),
#                  formula = y ~ x, parse = TRUE, 
#                  color = "black", size = 4, 
#                  label.y = 0.1, label.x = 0.5) +
#     
#     facet_wrap(~ input_data, scales = "fixed") + # Custom facet labels
#     theme_bw() +
#     theme(
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_line(colour = "gray", size = 0.5),
#       panel.border = element_rect(colour = "black", fill = NA),
#       strip.background = element_blank(),
#       strip.text = element_text(size = 15),
#       plot.title = element_text(hjust = 0.5, size = 20),
#       axis.text = element_text(size = 12.5),
#       axis.text.x = element_text(angle = 45, hjust = 1),
#       axis.title = element_text(size = 15),
#       legend.text = element_text(size = 15),
#       legend.title = element_text(size = 15)
#     ) +
#     labs(
#       title = " ",
#       x = expression(Observed ~  weight ~ (g ~ plant^{-1})), 
#       y = expression(Predicted ~  weight ~ (g ~ plant^{-1})) 
#     )
#   return(plot)
# }




