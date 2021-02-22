sessionInfo()

# Initialise system ############################################################
################################################################################.

# Cleaning
rm(list = ls()); graphics.off()

# max.print
options(max.print = 500)

# Libraries
library(parallel)
library(tidyverse)
library(rstan)
library(bipartite)
library(wrswoR)

# fail in robustness function (bipartite version 2.14)
if (packageVersion("bipartite") == "2.14"){
  robustness <- function (object) {
    if (!is(object, "bipartite")) { # "!" was missing here
      stop("This function cannot be meaningfully applied to objects of this class.") }
    N <- colSums(object)
    if (attr(object, "exterminated") == "lower") {
      y <- -object[, 3]
    }
    else {
      y <- -object[, 2]
    }
    y <- (sum(y) - cumsum(y))/sum(y)
    x <- (object[, "no"]/max(object[, "no"]))
    ext.curve <- splinefun(x, y)
    ext.area <- integrate(ext.curve, 0, 1)
    return(as.numeric(ext.area[[1]]))
  }
}


print(paste("#cores =", detectCores()))

rstan_options(auto_write = TRUE)
options(mc.cores = round(detectCores()/2, 0))

# get run ID and correction
perm_i <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
corr_i <- commandArgs(trailingOnly = TRUE)[2]

# get modularity and robustness arguments (if present)
mod_i <- T
out_full <- F
data_source <- "swnet_data.Rdata"
round_max <- 999
if (length(commandArgs(trailingOnly = TRUE)) > 2) mod_i <- as.logical(commandArgs(trailingOnly = TRUE)[3])
if (length(commandArgs(trailingOnly = TRUE)) > 3) rob_i <- as.logical(commandArgs(trailingOnly = TRUE)[4])
if (length(commandArgs(trailingOnly = TRUE)) > 4) out_full <- as.logical(commandArgs(trailingOnly = TRUE)[5])
if (length(commandArgs(trailingOnly = TRUE)) > 5) data_source <- as.character(commandArgs(trailingOnly = TRUE)[6])
if (length(commandArgs(trailingOnly = TRUE)) > 6) round_max <- as.numeric(commandArgs(trailingOnly = TRUE)[7])

# set working directory
setwd("XXX")

# change extinction function in bipartite (to prevent a bug)
source("f.extinction.R")
environment(extinction) <- asNamespace('bipartite')
assignInNamespace("extinction", extinction, ns = "bipartite")

# set parameters
responses <- c("nested", "qmod", "rob.A")

# load data
load(data_source)

# print settings
cat(paste("perm_i =", perm_i, 
          "\ncorr_i =", corr_i,
          "\nmod_i =", mod_i,
          "\nround_max =", round_max,
          "\ntimestamp:", Sys.time()))

# Define functions #############################################################
################################################################################.

f.null.models.swnet <- function (plots, 
                                 int_global, mod) {
  
  

  
  f.null.model <- function(pl_i){
    
    int_matrix_plot <- int_global %>% 
      filter(PlotID == pl_i) %>% 
      select(-PlotID) %>% 
      spread(Herb_Name_std, interaction, fill = 0) %>% 
      column_to_rownames("Host_Name_std") %>% 
      as.data.frame
    

    if (corr_i == "pres.abs") int_matrix_plot[int_matrix_plot > 0] <- 1
    
    
    n_LL <- nrow(int_matrix_plot)
    n_UL <- ncol(int_matrix_plot)
    n_int <- length(int_matrix_plot[int_matrix_plot > 0])
    
    cont <- TRUE
    nm.solution <- FALSE
    n_rep <- 1
    
    
    probs <- int_global %>%
      filter(grepl(substr(pl_i, 1, 3), PlotID)) %>% # only consider interactions of the region
      group_by(PlotID, Herb_Name_std) %>%
      summarise(spec_abund = sum(interaction)) %>%
      ungroup() %>%
      group_by(PlotID) %>%
      mutate(spec_abund_std = spec_abund / max(spec_abund)) %>%
      ungroup() %>%
      group_by(Herb_Name_std) %>%
      summarise(abundsum = sum(spec_abund_std)) %>%
      ungroup()
    
    while (cont) { # if a drawing round does not yield a compatible null matrix
      
      int_global_sample <- int_global %>%
        filter(grepl(substr(pl_i, 1, 3), PlotID)) %>% # only consider interactions of the region
        mutate(ID = paste(Host_Name_std, Herb_Name_std)) %>%
        group_by(Herb_Name_std) %>% 
        mutate(n_int = n()) %>% 
        ungroup() %>%
        left_join(probs, by = "Herb_Name_std") %>%
        mutate(prob = abundsum / n_int) 
      

      # random sample (probability according to overall abundance of the herbivore)
      if (corr_i == "pres.abs"){
        probs <- int_global_swnet_raw %>%
          filter(grepl(substr(pl_i, 1, 3), PlotID)) %>% # only consider interactions of the region
          group_by(PlotID, Herb_Name_std) %>%
          summarise(spec_abund = sum(interaction)) %>%
          ungroup() %>%
          group_by(PlotID) %>%
          mutate(spec_abund_std = spec_abund / max(spec_abund)) %>%
          ungroup() %>%
          group_by(Herb_Name_std) %>%
          summarise(abundsum = sum(spec_abund_std)) %>%
          ungroup()
        
        int_global_raw <- int_global_swnet_raw %>%
          filter(grepl(substr(pl_i, 1, 3), PlotID)) %>% # only consider interactions of the region
          mutate(ID = paste(Host_Name_std, Herb_Name_std)) %>%
          group_by(Herb_Name_std) %>% 
          mutate(n_int = n()) %>% 
          ungroup() %>% 
          left_join(probs, by = "Herb_Name_std") %>%
          mutate(prob = abundsum / n_int) 
        
        int_global_sample <- 
          int_global_sample[sample_int_rank(nrow(int_global_sample), 
                                   nrow(int_global_sample),
                                   prob = int_global_raw$prob), ] %>% 
          filter(!duplicated(ID)) # no double interactions allowed
        
      } else {
        int_global_sample <- 
          int_global_sample[sample_int_rank(nrow(int_global_sample), 
                                   nrow(int_global_sample), 
                                   prob = int_global_sample$prob), ] %>% 
          filter(!duplicated(ID)) # no double interactions allowed
      }

      
      int_matrix_plot_null <- int_global_sample[c(1:min(n_LL, n_UL)), ]
      int_global_sample <- int_global_sample[-c(1:min(n_LL, n_UL)), ]
      
      repeat {
        int_add <- int_global_sample[1, ]
        if(length(unique(c(int_matrix_plot_null$Host_Name_std, 
                           int_add$Host_Name_std))) <= n_LL &
           length(unique(c(int_matrix_plot_null$Herb_Name_std, 
                           int_add$Herb_Name_std))) <= n_UL) {
          int_matrix_plot_null <- rbind(int_matrix_plot_null, int_add)
        } 
        
        int_global_sample <- int_global_sample[-1, ]
        
        # if plants are "filled", remove all the other ones from the sample data frame
        if(length(unique(int_matrix_plot_null$Host_Name_std)) == n_LL){
          int_global_sample <- int_global_sample %>%
            filter(Host_Name_std %in% int_matrix_plot_null$Host_Name_std)
        }
        
        # if herbivores are "filled", remove all the other ones from the sample data frame
        if(length(unique(int_matrix_plot_null$Herb_Name_std)) == n_UL){
          int_global_sample <- int_global_sample %>%
            filter(Herb_Name_std %in% int_matrix_plot_null$Herb_Name_std)
        }
        
        
        # same #interactions break
        if(nrow(int_matrix_plot_null) == n_int | 
           nrow(int_global_sample) == 0) {
          break
        }

      }
      
      # same #interactions break
      if(nrow(int_matrix_plot_null) == n_int &
         length(unique(int_matrix_plot_null$Host_Name_std)) == n_LL &
         length(unique(int_matrix_plot_null$Herb_Name_std)) == n_UL) {
        cont <- FALSE
        nm.solution <- TRUE
        # break
      } else {
        n_rep <- n_rep + 1
        if (n_rep > round_max){
          cont <- FALSE
          print(paste("No solution found for", pl_i))
        } else {
          if (n_rep %% 10 == 0){
            print(paste(pl_i, "goes to round", n_rep, ";",
                        nrow(int_matrix_plot_null), "of", n_int, "interactions /",
                        length(unique(int_matrix_plot_null$Host_Name_std)), "of", n_LL, "hosts /",
                        length(unique(int_matrix_plot_null$Herb_Name_std)), "of", n_UL, "herbivores"))
          }
        }
      }
      
      
    }
    
    if (nm.solution){
      int_matrix_plot_null <- int_matrix_plot_null %>%
        select(Host_Name_std, Herb_Name_std, interaction) %>%
        spread(Herb_Name_std, interaction, fill = 0) %>%
        column_to_rownames("Host_Name_std")
      
      if (corr_i == "pres.abs"){ # for presence-absence networks (weighted equals zero)
        nested.target <- nested(int_matrix_plot_null, method = "NODF2")
      } else {
        nested.target <- nested(int_matrix_plot_null, method = "weighted NODF")
      }

      scnd.ext <- second.extinct(int_matrix_plot_null, participant = "lower",
                                 method = "abund")
      robA.target <- robustness(scnd.ext)
      
    if (mod){
        modules.target <- computeModules(int_matrix_plot_null)
        qmod <- modules.target@likelihood
    }
    } else {
      nested.target <- NA
      robA.target <- NA
      qmod <- NA
    }
    
    
    out <- data.frame(PlotID = pl_i,
                      nested = nested.target,
                      rob.A = robA.target,
                      stringsAsFactors = F)
    

    if(mod) out <- cbind(out, qmod)
    
    rownames(out) <- NULL
    out
  }
  
  data.nullmodel <- mclapply(plots, f.null.model)
  data.nullmodel <- do.call(rbind, data.nullmodel)
  
  
  # output ---------------------------------------------------------------------.
  data.nullmodel
  


}


# hierarchical models
stan_lui_mod <- stan_model(file = "scr_stan/Stan_LUI_model.stan")
stan_luisplit_mod <- stan_model(file = "scr_stan/Stan_LUIsplit_model.stan")

hier_mod_swnet_null <- function(data, datatype){
  
  data <- data %>% 
    filter(!is.na(nested)) %>% # Exclude plots for which no null model solution could be obtained
    left_join(lui.glob, by = "PlotID")
  
  out <- list()
  for (res in responses[responses %in% names(data)]){
    
    data$response <- (data[, res] - mean(nw.metrics.swnet[[datatype]][, res]))/sd(nw.metrics.swnet[[datatype]][, res])
    
    data.comb <- list(J = length(unique(data$Exploratory)),
                      N = nrow(data),
                      Exploratory = as.integer(as.factor(data$Exploratory)),
                      LUI_z = as.vector(data$LUI_z),
                      y = as.vector(data$response))
    
    fit <- sampling(object = stan_lui_mod, data = data.comb,
                    chains = 4, warmup = 500, iter = 2000, control = list(adapt_delta = 0.99))
    if (out_full) {
      out[[res]] <- list(lui_comb = extract(fit)[c("a", "b", "mu_a", "sigma_a", "sigma_b", "sigma_y")])
    } else {
      out[[res]] <- list(lui_comb = extract(fit)[c("a", "b")])
    }
    
    data.split <- list(J = length(unique(data$Exploratory)),
                       N = nrow(data),
                       Exploratory = as.integer(as.factor(data$Exploratory)),
                       M_z = as.vector(data$M_z),
                       F_z = as.vector(data$F_z),
                       G_z = as.vector(data$G_z),
                       y = as.vector(data$response))
    
    
    fit <- sampling(object = stan_luisplit_mod, data = data.split,
                    chains = 4, warmup = 500, iter = 2000, control = list(adapt_delta = 0.99))
    
    if (out_full) {
      out[[res]] <- c(out[[res]], list(lui_split = extract(fit)[c("a", "b", "mu_a", "sigma_a", "sigma_b_m", "sigma_b_f", "sigma_b_g", "sigma_y")]))
    } else {
      out[[res]] <- c(out[[res]], list(lui_split = extract(fit)[c("a", "b")]))
    }
    
    
    # cleaning (from https://github.com/stan-dev/rstan/issues/448) -------------.
    dso_filename <-  stan_lui_mod@dso@dso_filename
    loaded_dlls  <-  getLoadedDLLs()
    if (dso_filename %in% names(loaded_dlls)) {
      message("Unloading DLL for model dso ", dso_filename)
      model.dll  <-  loaded_dlls[[dso_filename]][['path']]
      dyn.unload(model.dll)
    } 
    
    dso_filename <-  stan_luisplit_mod@dso@dso_filename
    loaded_dlls  <-  getLoadedDLLs()
    if (dso_filename %in% names(loaded_dlls)) {
      message("Unloading DLL for model dso ", dso_filename)
      model.dll  <-  loaded_dlls[[dso_filename]][['path']]
      dyn.unload(model.dll)
    } 
    
    loaded_dlls <-  getLoadedDLLs()
    loaded_dlls  <-  loaded_dlls[str_detect(names(loaded_dlls), '^file')]
    if (length(loaded_dlls) > 10) {
      for (dll in head(loaded_dlls, -10)) {
        dyn.unload(dll[['path']])
      }
    }
    
  }
  
  out
}

# Network metrics ##############################################################
################################################################################.

assign("int_global_target", eval(parse(text  = paste0("int_global_swnet_", gsub("\\.", "_", corr_i)))))

nw.metrics_target <- f.null.models.swnet(plots = plots, 
                                         int_global = int_global_target, 
                                         mod = mod_i)

# Analyses  ####################################################################
################################################################################.


out <- hier_mod_swnet_null(nw.metrics_target, corr_i)

out %>%
  map(cbind.data.frame)  %>%
  enframe("response") %>% 
  unnest() %>% 
  write.table(paste0("out_data/swnet_", gsub("\\.", "_", corr_i), "/data_",
                     formatC(perm_i, width = 5, format = "d", flag = "0"), ".txt"))



