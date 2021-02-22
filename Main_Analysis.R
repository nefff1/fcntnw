# INITIAL PREPARATIONS ----------------------------------- #####################
################################################################################.

# read packages
library(cowplot)
library(bipartite)
library(magick)
library(nlme)
library(parallel)
library(tidyverse); theme_set(theme_classic())
library(furrr)
library(rstan)
library(bayestestR)
library(lmerTest)

# fix functions
select <- dplyr::select
rename <- dplyr::rename

# set global options for Stan
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
options(stringsAsFactors = F)

# read Stan models
stan_lui_mod <- stan_model(file = "Stan_models/Stan_LUI_model.stan")
stan_luisplit_mod <- stan_model(file = "Stan_models/Stan_LUIsplit_model.stan")
stan_lui_mod_poiss <- stan_model(file = "Stan_models/Stan_LUI_model_poisson.stan")
stan_luisplit_mod_poiss <- stan_model(file = "Stan_models/Stan_LUIsplit_model_poisson.stan")

stan_formi_mod <- stan_model(file = "Stan_models/Stan_ForMI_model.stan")
stan_formisplit_mod <- stan_model(file = "Stan_models/Stan_ForMIsplit_model.stan")
stan_formi_mod_poiss <- stan_model(file = "Stan_models/Stan_ForMI_model_poisson.stan")
stan_formisplit_mod_poiss <- stan_model(file = "Stan_models/Stan_ForMIsplit_model_poisson.stan")

stan_type_mod <- stan_model(file = "Stan_models/Stan_type_model.stan")

# set main response variables (network metrics)
responses <- c("nested", "qmod", "rob.A")

# stage weights for network construction
stage_weights <-  c(juvenil = 1, adult = 1)

# set order of host taxonomic ranks:

host.rank.levels <- c("Polyphyletic", "Kingdom",
                      "Phylum", "Subdivision", "Class", 
                      "Cladus1", "Order", "Family", 
                      "Tribus", "Genus", "Species",
                      "Subspecies")

# data used to match polyphagous herbivores and their hosts:

pl.groups.def <- data.frame(Name = c("tree", "forb", "grass", "shrub", "fern", "clubmoss", 
                                     "epiphyte", "conifer", "deciduous", "perennial",
                                     "dwarf shrub", "deciduous tree", "deciduous shrub",
                                     "fruit tree", "flower"),
                            Level = c(1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 5)) %>% 
  mutate(Level = paste0("Host_Group", Level))


# READ DATA --------------------------------------------- ######################
################################################################################.

# interaction database ---------------------------------------------------------.
# available from the BEXIS repository (https://www.bexis.uni-jena.de), ID 26926

int.db <- read.delim("Data/26926.txt", header = T, sep = "\t")

int.db <- int.db %>%
  mutate(Host_Rank_red = paste0(toupper(substr(Host_Rank, 1, 1)),
                                substr(Host_Rank, 2, nchar(Host_Rank))),
         Host_Rank_red = ifelse(Host_Rank_red %in% c("Variety", "Form"),
                                "Subspecies", Host_Rank_red),
         Host_Rank_red = ifelse(Host_Rank_red %in% c("Aggregate", "Section"),
                                "Species", Host_Rank_red),
         Host_Rank_red = ifelse(Host_Rank_red == "Polyphagous", NA, Host_Rank_red)) %>%
  select(Herb_Name_std, Herb_Stage, Host_Name_std, Host_Genus_Species, Host_Rank_red) %>%
  rename(stage = Herb_Stage)

# make taxonomic rank an ordered factor
int.db <- int.db %>% 
  mutate(Host_Rank_red = ordered(Host_Rank_red, levels = host.rank.levels))

# collated data ----------------------------------------------------------------.
# available from the BEXIS repository (https://www.bexis.uni-jena.de), ID 30902

data_coll <- read.delim("Data/30902.txt", header = T, sep = "\t")

# extract arthropod data -------------------------------------------------------.
 
herb_agg <- data_coll %>% 
  filter(Kingdom == "Animalia") %>% 
  select(PlotID, Name_std, NumberAdults) %>% 
  rename(Herb_Name_std = "Name_std")

# define grassland subset (sweep netting):

swnet_agg <- herb_agg %>% 
  filter(grepl("G", PlotID))

# define forest subset (window traps):

wintr_agg <- herb_agg %>% 
  filter(grepl("W", PlotID))

rm(herb_agg)

# extract a herbivore species list including feeding guild and body size -------.

herb_splist <- data_coll %>% 
  filter(Kingdom == "Animalia") %>% 
  select(Name_std, Genus, Family, Suborder, Order, 
         Feeding_guild_short, Body_Size) %>% 
  distinct() %>% 
  rename_at(vars(-c(Feeding_guild_short, Body_Size)),
            ~paste0("Herb_", .))

# extract plant cover and biomass data -----------------------------------------.

plant_agg <- data_coll %>% 
  filter(Kingdom == "Plantae") %>% 
  select(PlotID, Name_std, Rank_red, Kingdom, Phylum, Subdivision,
         Class, Cladus1, Order, Family, Tribus, Genus, Species, Subspecies,
         Group1, Group2, Group3, Group4, Group5, Cover, biomass) %>% 
  rename_at(vars(-c(PlotID, Cover, biomass)),
            ~paste0("Host_", .))

# make taxonomic rank an ordered factor
plant_agg <- plant_agg %>% 
  mutate(Host_Rank_red = ordered(Host_Rank_red, levels = host.rank.levels))

# define grassland subset:

plant_agg_g <- plant_agg %>% 
  filter(grepl("G", PlotID))


# define forest subset:

plant_agg_f <- plant_agg %>% 
  filter(grepl("W", PlotID))

rm(plant_agg)

# extract land use data --------------------------------------------------------.

landuse <- data_coll %>% 
  select(PlotID, Exploratory,
         LUI_z, M_z, F_z, G_z, 
         ForMI_z, Inonat_z, Idwcut_z, Iharv_z) %>% 
  distinct()

# define grassland subset (LUI):

lui.glob <- landuse %>% 
  filter(grepl("G", PlotID)) %>% 
  select(-c(ForMI_z, Inonat_z, Idwcut_z, Iharv_z))

# define forest subset (ForMI):

formi <- landuse %>% 
  filter(grepl("W", PlotID)) %>% 
  select(-c(LUI_z, M_z, F_z, G_z))

rm(landuse)
rm(data_coll)

# Insect length / weight allometric relation -----------------------------------.
# data from Sohlström et al. (2018), Dryad Digital Repository, DOI 10.5061/dryad.vk24fr1
# and from Sohlström et al. (2018), Ecology and Evolution, DOI 10.1002/ece3.4702

regs_mass_length <- read.table("Data/regs_mass_length.txt", header = T) 

# Columns: *Herb_Name_std* (standardized taxon name), *a* (intercept),
# *b* (slope)

# ANALYSES ---------------------------------------------- ######################
################################################################################.

#  functions ###################################################################
################################################################################.

# unique function not returning NAs --------------------------------------------.

unique_noNA <- function(x) {
  x <- unique(x)
  x <- x[!is.na(x)]
  x
}

# extract allometric regression coefficient for herbivores ---------------------.

f_lw <- function(herb_i, data){
  genus_i <- unique(herb_splist$Herb_Genus[herb_splist$Herb_Name_std == herb_i])
  family_i <- unique(herb_splist$Herb_Family[herb_splist$Herb_Name_std == herb_i])
  suborder_i <- unique(herb_splist$Herb_Suborder[herb_splist$Herb_Name_std == herb_i])
  order_i <- unique(herb_splist$Herb_Order[herb_splist$Herb_Name_std == herb_i])
  
  if (genus_i %in% data$Herb_Name_std) {
    out <- data.frame(Herb_Name_std = herb_i,
                      b = data$b[data$Herb_Name_std == genus_i],
                      a = data$a[data$Herb_Name_std == genus_i],
                      stringsAsFactors = F)
    out
  } else if (family_i %in% data$Herb_Name_std) {
    out <- data.frame(Herb_Name_std = herb_i,
                      b = data$b[data$Herb_Name_std == family_i],
                      a = data$a[data$Herb_Name_std == family_i],
                      stringsAsFactors = F)
    out
  } else if (suborder_i %in% data$Herb_Name_std) {
    out <- data.frame(Herb_Name_std = herb_i,
                      b = data$b[data$Herb_Name_std == suborder_i],
                      a = data$a[data$Herb_Name_std == suborder_i],
                      stringsAsFactors = F)
    out
  } else if (order_i %in% data$Herb_Name_std) {
    out <- data.frame(Herb_Name_std = herb_i,
                      b = data$b[data$Herb_Name_std == order_i],
                      a = data$a[data$Herb_Name_std == order_i],
                      stringsAsFactors = F)
    out
  }
}

# construct global interaction data base ---------------------------------------.


construct_int_matrix <- function(herb_i, plant_agg, tax_weights, poly_weight, poly_all = T, poly_sub = T){
  int_target_i <- int.db %>% 
    filter(Herb_Name_std == herb_i,
           !is.na(Host_Genus_Species))
  
  f.stage <- function(int_target){
    if (nrow(int_target) > 0){
      plants_herb_i <- data.frame()
      
      
      # 1. Hosts based on taxonomy
      for (plant_i in unique_noNA(int_target$Host_Name_std)){
        
        host_rank_target <- unique_noNA(int_target$Host_Rank_red[int_target$Host_Name_std == plant_i])
        host_column_target <- paste0("Host_", host_rank_target)
        
        plants_herb_i <- plant_agg %>% 
          filter(!! sym(host_column_target) == plant_i) %>% 
          group_by(Host_Name_std) %>% 
          summarise() %>% 
          ungroup() %>% 
          mutate(weight = tax_weights$weight[tax_weights$group == host_rank_target]) %>% 
          bind_rows(plants_herb_i, .)
        
        # also include Genus reports of this plant species (particularly important for Taraxacum!!)
        if (host_rank_target > "Genus"){
          genus_target <- unique(int_target$Host_Genus[int_target$Host_Name_std == 
                                                       plant_i])
          host_rank_target <- "Genus"
          
          plants_herb_i <- plant_agg %>% 
            filter(Host_Rank_red == "Genus",
                   Host_Genus == genus_target) %>% 
            group_by(Host_Name_std) %>% 
            summarise() %>% 
            ungroup() %>% 
            mutate(weight = tax_weights$weight[tax_weights$group == host_rank_target]) %>% 
            bind_rows(plants_herb_i, .)
        }
      }
      
      # 2. Hosts of polyphagous species
      poly <- int_target$Host_Genus_Species[grepl("polyphagous", 
                                                  int_target$Host_Genus_Species)]
      
      if (any(poly == "polyphagous" & poly_all)) {
        plants_herb_i <- data.frame(Host_Name_std = sort(unique(plant_agg$Host_Name_std)),
                                    weight = poly_weight,
                                    stringsAsFactors = F) %>% 
          bind_rows(plants_herb_i, .)
      } else if (poly_sub) {
        poly <- poly[poly != "polyphagous"]
        poly_groups <- data.frame(Name = gsub("polyphagous - ", "", poly),
                                  stringsAsFactors = F)
        poly_groups <- poly_groups %>% 
          left_join(pl.groups.def, by = "Name")
        
        for (poly_i in poly_groups$Name){
          host_column_target <- unique(poly_groups$Level[poly_groups$Name == poly_i])
          
          plants_herb_i <- plant_agg %>% 
            filter(!! sym(host_column_target) == poly_i) %>% 
            group_by(Host_Name_std) %>% 
            summarise() %>% 
            ungroup() %>% 
            mutate(weight = poly_weight) %>% 
            bind_rows(plants_herb_i, .)
        }
      }
      
      # exclude duplicates (due to e.g. both lower and higher order plants in int.db, e.g. family and species)
      plants_herb_i <- plants_herb_i %>% 
        arrange(-weight) %>% # so that lower weights are excluded
        filter(!duplicated(Host_Name_std))
      
      
      
      if(nrow(plants_herb_i) > 0){
        plants_herb_i <- plants_herb_i %>% 
          spread(Host_Name_std, weight)
        
        out <- plants_herb_i %>% 
          bind_cols(Herb_Name_std = herb_i, .) %>% 
          as.data.frame()
        
      } else {
        out <- data.frame(Herb_Name_std = herb_i, stringsAsFactors = F)
      }
      
    } else {
      out <- data.frame(Herb_Name_std = herb_i, stringsAsFactors = F)
    }
    
    out
  }

  out_juv <- int_target_i %>%
    filter(stage %in% c("juvenil", "both") |
             is.na(stage)) %>%
    f.stage()

  out_ad <- int_target_i %>%
    filter(stage %in% c("adult", "both") |
             is.na(stage)) %>%
    f.stage()

  out <- bind_rows(out_ad, out_juv) %>% 
      mutate_all(~ ifelse(is.na(.), 0, .)) %>% 
      add_column(stage = c("adult", "juvenil")) %>% 
      select(Herb_Name_std, stage, everything())
  
  out
}


# global collection of all interactions ----------------------------------------.

f.int.global <- function(pl_i, int_matrix, herb_ct, plant_ct, stage_weights) {
  
  if (any(colnames(int_matrix) !=   names(unlist(plant_ct[pl_i, ])))){
    warning("plant_ct and int_matrix not covering the same plants")
    stop()
  }

  int_plant_target <- unlist(plant_ct[pl_i, ]) * t(int_matrix)
  int_plant_target <- int_plant_target[ , colSums(int_plant_target) > 0]

  herb_ct_juv <- herb_ct[pl_i, ] * stage_weights["juvenil"] / sum(stage_weights)
  names(herb_ct_juv) <- paste(names(herb_ct_juv), "juvenil", sep = "_")
  
  herb_ct_ad <- herb_ct[pl_i, ] * stage_weights["adult"] / sum(stage_weights)
  names(herb_ct_ad) <- paste(names(herb_ct_ad), "adult", sep = "_")
  
  herb_ct_target <- cbind(herb_ct_juv, herb_ct_ad)
  herb_ct_target <- herb_ct_target[, colnames(int_plant_target)]
  
  if (any(names(herb_ct_target) != colnames(int_plant_target))){
    warning("herb_ct_target and int_matrix not covering the same insects")
    stop()
  }
  
  int_matrix_plot <- int_plant_target %*% diag(1 / colSums(int_plant_target)) %*% diag(herb_ct_target)
  colnames(int_matrix_plot) <- colnames(int_plant_target)
  int_matrix_plot <- as.data.frame(int_matrix_plot)
  int_matrix_plot <- int_matrix_plot[rowSums(int_matrix_plot) > 0, colSums(int_matrix_plot) > 0]
  
  # summarise rows per species:
  
  sp.list <- unique(substr(names(int_matrix_plot), 1, 
                           regexpr("_", names(int_matrix_plot)) - 1))
  
  # create data frame for condensed data (per species)
  int_matrix_plot_cond <- data.frame(row.names = rownames(int_matrix_plot))
  
  for (sp_i in sp.list){
    
    # if there are entries for both stages
    if (sum(grepl(sp_i, names(int_matrix_plot))) == 2){
      add <- data.frame(rowSums(int_matrix_plot[, grepl(sp_i, 
                                                        names(int_matrix_plot))]))
      names(add) <- sp_i
      
      int_matrix_plot_cond <- cbind(int_matrix_plot_cond,
                                    add)
    } else if (sum(grepl(sp_i, names(int_matrix_plot))) == 1){
      
      add <- int_matrix_plot[, grepl(sp_i, names(int_matrix_plot)), drop = F]
      
      
      # check whether both stages are actually herbivorous
      stage.n <- int.db %>% 
        filter(Herb_Name_std == sp_i) %>% 
        select(stage) %>% 
        unique() %>% nrow
      
      # if both would be herbivorous, remove stage-weighting
      if (stage.n > 1) {
        stage_present <- substr(names(add), regexpr("_", names(add)) + 1,
                                nchar(names(add)))
        
        add <- add / stage_weights[stage_present] * sum(stage_weights)
        
      }
      
      names(add) <- sp_i
      
      int_matrix_plot_cond <- cbind(int_matrix_plot_cond,
                                    add)
      
    } else {
      warning("Something wrong")
    }
  }
  
  int_matrix_plot <- as.data.frame(int_matrix_plot_cond)
  
  # exclude very small interactions when all plants are included in very small abundances:
  if(all(plant_ct > 0)){
    int_matrix_plot[int_matrix_plot < 10^-10] <- 0 # empirically found threshold to define the ditch between two clearly separated groups 
  }
  
  # exclude empty rows / columns. this will affect some metrics! (e.g. nestedness)
  int_matrix_plot <- int_matrix_plot[rowSums(int_matrix_plot) > 0, colSums(int_matrix_plot) > 0]
  

  int_matrix_plot %>% 
    rownames_to_column("Host_Name_std") %>% 
    gather(Herb_Name_std, interaction, -Host_Name_std) %>% 
    filter(interaction > 0) %>% 
    data.frame(PlotID = pl_i, ., stringsAsFactors = F)
}


# network plotting and metrics -------------------------------------------------.

f.nw.metrics <- function(pl_i, int_global,
                         mod = T, 
                         colset_herb = F, colset_host = F) {

  int_matrix_plot <- int_global %>% 
    filter(PlotID == pl_i) %>% 
    select(-PlotID) %>% 
    spread(Herb_Name_std, interaction, fill = 0) %>% 
    column_to_rownames("Host_Name_std") %>% 
    as.data.frame()
  
  # exclude empty rows / columns. this will affect some metrics! (e.g. nestedness) (should already be the case, just to be sure)
  int_matrix_plot <- int_matrix_plot[rowSums(int_matrix_plot) > 0, colSums(int_matrix_plot) > 0]

  
  if (length(unique(unlist(int_matrix_plot))) == 2){ # for presence-absence networks (weighted equals zero)
    nested.target <- nested(int_matrix_plot, method = "NODF2")
  } else {
    nested.target <- nested(int_matrix_plot, method = "weighted NODF")
  }
  
  scnd.ext <- second.extinct(int_matrix_plot, participant = "lower",
                             method = "abund")
  robA.target <- robustness(scnd.ext)
  
  if (mod){
    modules.target <- computeModules(int_matrix_plot)
    qmod <- modules.target@likelihood
  }
  
  out <- data.frame(PlotID = pl_i,
                    nested = nested.target,
                    rob.A = robA.target,
                    stringsAsFactors = F)
  
  if(mod) out <- cbind(out, qmod)
  
  rownames(out) <- NULL
  
  out
  
}


#  GRASSLANDS / SWEEP NETS ------------------------------- #####################
################################################################################.

#  prepare data and interaction matrices #######################################
################################################################################.

# herbivores -------------------------------------------------------------------.

herb_list_swnet <- swnet_agg %>% 
  select(Herb_Name_std) %>% 
  arrange(Herb_Name_std) %>% 
  unique()

swnet_ct <- swnet_agg %>% 
  select(PlotID, Herb_Name_std, NumberAdults) %>% 
  spread(Herb_Name_std, NumberAdults, fill = 0) %>% 
  as.data.frame() %>% 
  column_to_rownames("PlotID")

# weight omnivore and mostly carnivore/fungivore species less:

feeding_weights <- herb_list_swnet %>% 
  left_join(herb_splist, by = "Herb_Name_std") %>% 
  mutate(feeding_weight = ifelse(Feeding_guild_short == "h", 1,
                                 ifelse(Feeding_guild_short == "o",
                                        0.5, 0.25))) %>% 
  select(Herb_Name_std, feeding_weight)
swnet_ct_fw <- as.data.frame(t(feeding_weights$feeding_weight * t(swnet_ct)))

# weight with metabolic rate:
lw_swnet <- herb_list_swnet$Herb_Name_std %>% 
  future_map(f_lw, regs_mass_length) %>% 
  do.call(rbind, .)

lw_swnet <- herb_list_swnet %>% 
  left_join(lw_swnet, by = "Herb_Name_std") %>% 
  left_join(select(herb_splist, Herb_Name_std, Body_Size), by = "Herb_Name_std") 


lw_swnet <- lw_swnet %>% 
  mutate(weight = 10 ^ a * (Body_Size ^ b),
         metabolic_rel = weight ^ 0.75)

swnet_ct_mr <- 
  swnet_agg %>% 
  filter(Herb_Name_std %in% herb_list_swnet$Herb_Name_std) %>%
  left_join(lw_swnet, by = "Herb_Name_std") %>%
  mutate(metabolic_entities = NumberAdults * metabolic_rel) %>% 
  select(PlotID, Herb_Name_std, metabolic_entities) %>% 
  spread(Herb_Name_std, metabolic_entities, fill = 0) %>% 
  as.data.frame() %>% 
  column_to_rownames("PlotID")


# combine feeding weighting and metabolic rate:
swnet_ct_fw_mr <- as.data.frame(t(feeding_weights$feeding_weight * t(swnet_ct_mr)))


# plants -----------------------------------------------------------------------.

plant_ct_g <- plant_agg_g %>% 
  filter(PlotID %in% rownames(swnet_ct)) %>% 
  select(PlotID, Host_Name_std, Cover) %>% 
  spread(Host_Name_std, Cover, fill = 0) %>% 
  as.data.frame() %>% 
  column_to_rownames("PlotID")

# biomass instead of cover:
plant_ct_g_bm <- plant_agg_g %>% 
  filter(PlotID %in% rownames(swnet_ct)) %>% 
  select(PlotID, Host_Name_std, biomass) %>% 
  spread(Host_Name_std, biomass, fill = 0) %>% 
  as.data.frame() %>% 
  column_to_rownames("PlotID")

# all plants present in very low numbers
plant_ct_g_nozero <- plant_ct_g
plant_ct_g_nozero[plant_ct_g_nozero == 0] <- 10^-25

# global interaction matrix ----------------------------------------------------.

# weighting all taxonomic level information the same
tax_weights_nw <- data.frame(group = ordered(host.rank.levels, 
                                             levels = host.rank.levels))
tax_weights_nw$weight <- 1
poly_weight_nw <- 1


int_mat_swnet <- mclapply(herb_list_swnet$Herb_Name_std, 
                          construct_int_matrix, plant_agg_g, 
                          tax_weights_nw, poly_weight_nw)
int_mat_swnet <- do.call(bind_rows, int_mat_swnet)
int_mat_swnet[is.na(int_mat_swnet)] <- 0
int_mat_swnet <- int_mat_swnet[ , sort(names(int_mat_swnet))]
int_mat_swnet <- int_mat_swnet %>% 
  mutate(ID = paste(Herb_Name_std, stage, sep = "_")) %>% 
  column_to_rownames("ID") %>% 
  select(-c(Herb_Name_std, stage))


# weighting  different taxonomic level information differently
tax_weights_pw <- data.frame(group = ordered(host.rank.levels, 
                                             levels = host.rank.levels))
tax_weights_pw <- tax_weights_pw %>% 
  mutate(weight = ifelse(group < "Family", 0.25,
                         ifelse(group < "Genus", 0.5, 1)))
poly_weight_pw <- 0.25

int_mat_swnet_pw <- mclapply(herb_list_swnet$Herb_Name_std, 
                             construct_int_matrix, plant_agg_g, tax_weights_pw, poly_weight_pw)
int_mat_swnet_pw <- do.call(bind_rows, int_mat_swnet_pw)
int_mat_swnet_pw[is.na(int_mat_swnet_pw)] <- 0
int_mat_swnet_pw <- int_mat_swnet_pw[ , sort(names(int_mat_swnet_pw))]
int_mat_swnet_pw <- int_mat_swnet_pw %>% 
  mutate(ID = paste(Herb_Name_std, stage, sep = "_")) %>% 
  column_to_rownames("ID") %>% 
  select(-c(Herb_Name_std, stage))


# no polyphagous interactions
tax_weights_no_poly <- data.frame(group = ordered(host.rank.levels, 
                                                  levels = host.rank.levels))
tax_weights_no_poly <- tax_weights_no_poly %>% 
  mutate(weight = ifelse(group < "Family", 0 ,1))
int_mat_swnet_no_poly <- mclapply(herb_list_swnet$Herb_Name_std, 
                                  construct_int_matrix, plant_agg_g, 
                                  tax_weights_no_poly, poly_weight = 0,
                                  poly_all = F, poly_sub = F)
int_mat_swnet_no_poly <- do.call(bind_rows, int_mat_swnet_no_poly)
# add null columns:
nullcols  <-  unique(plant_agg_g$Host_Name_std[!plant_agg_g$Host_Name_std %in% names(int_mat_swnet_no_poly)])
int_mat_swnet_no_poly[, nullcols] <- NA; rm(nullcols)
int_mat_swnet_no_poly[is.na(int_mat_swnet_no_poly)] <- 0
int_mat_swnet_no_poly <- int_mat_swnet_no_poly[ , sort(names(int_mat_swnet_no_poly))]
int_mat_swnet_no_poly <- int_mat_swnet_no_poly %>% 
  mutate(ID = paste(Herb_Name_std, stage, sep = "_")) %>% 
  column_to_rownames("ID") %>% 
  select(-c(Herb_Name_std, stage))

# calculate global interaction matrices ########################################
################################################################################.

int_global_swnet_raw <- rownames(swnet_ct) %>% 
  map(f.int.global, int_mat_swnet, swnet_ct, plant_ct_g, stage_weights) %>% 
  do.call(rbind, .)

int_global_swnet_fw <- rownames(swnet_ct_fw) %>% 
  map(f.int.global, int_mat_swnet, swnet_ct_fw, plant_ct_g, stage_weights) %>% 
  do.call(rbind, .)

int_global_swnet_pw <- rownames(swnet_ct) %>% 
  map(f.int.global, int_mat_swnet_pw, swnet_ct, plant_ct_g, stage_weights) %>% 
  do.call(rbind, .)

int_global_swnet_mr <- rownames(swnet_ct_mr) %>% 
  map(f.int.global, int_mat_swnet, swnet_ct_mr, plant_ct_g, stage_weights) %>% 
  do.call(rbind, .)

int_global_swnet_bm <- rownames(swnet_ct) %>% 
  map(f.int.global, int_mat_swnet, swnet_ct, plant_ct_g_bm, stage_weights) %>% 
  do.call(rbind, .)

int_global_swnet_fw_pw <- rownames(swnet_ct_fw) %>% 
  map(f.int.global, int_mat_swnet_pw, swnet_ct_fw, plant_ct_g, stage_weights) %>% 
  do.call(rbind, .)

int_global_swnet_fw_mr <- rownames(swnet_ct_fw_mr) %>% 
  map(f.int.global, int_mat_swnet, swnet_ct_fw_mr, plant_ct_g, stage_weights) %>% 
  do.call(rbind, .)

int_global_swnet_fw_bm <- rownames(swnet_ct_fw) %>% 
  map(f.int.global, int_mat_swnet, swnet_ct_fw, plant_ct_g_bm, stage_weights) %>% 
  do.call(rbind, .)

int_global_swnet_pw_mr <- rownames(swnet_ct_mr) %>% 
  map(f.int.global, int_mat_swnet_pw, swnet_ct_mr, plant_ct_g, stage_weights) %>% 
  do.call(rbind, .)

int_global_swnet_pw_bm <- rownames(swnet_ct) %>% 
  map(f.int.global, int_mat_swnet_pw, swnet_ct, plant_ct_g_bm, stage_weights) %>% 
  do.call(rbind, .)

int_global_swnet_mr_bm <- rownames(swnet_ct_mr) %>% 
  map(f.int.global, int_mat_swnet, swnet_ct_mr, plant_ct_g_bm, stage_weights) %>% 
  do.call(rbind, .)

int_global_swnet_fw_pw_mr <- rownames(swnet_ct_fw_mr) %>% 
  map(f.int.global, int_mat_swnet_pw, swnet_ct_fw_mr, plant_ct_g, stage_weights) %>% 
  do.call(rbind, .)

int_global_swnet_fw_pw_bm <- rownames(swnet_ct_fw) %>% 
  map(f.int.global, int_mat_swnet_pw, swnet_ct_fw, plant_ct_g_bm, stage_weights) %>% 
  do.call(rbind, .)

int_global_swnet_fw_mr_bm <- rownames(swnet_ct_fw_mr) %>% 
  map(f.int.global, int_mat_swnet, swnet_ct_fw_mr, plant_ct_g_bm, stage_weights) %>% 
  do.call(rbind, .)

int_global_swnet_pw_mr_bm <- rownames(swnet_ct_mr) %>% 
  map(f.int.global, int_mat_swnet_pw, swnet_ct_mr, plant_ct_g_bm, stage_weights) %>% 
  do.call(rbind, .)

int_global_swnet_fw_pw_mr_bm <- rownames(swnet_ct_fw_mr) %>% 
  map(f.int.global, int_mat_swnet_pw, swnet_ct_fw_mr, plant_ct_g_bm, stage_weights) %>% 
  do.call(rbind, .)

int_global_swnet_no_poly <- rownames(swnet_ct) %>% 
  map(f.int.global, int_mat_swnet_no_poly, swnet_ct, plant_ct_g, stage_weights) %>% 
  do.call(rbind, .)

int_global_swnet_all_plants <- rownames(swnet_ct) %>% 
  map(f.int.global, int_mat_swnet, swnet_ct, plant_ct_g_nozero, stage_weights) %>% 
  do.call(rbind, .)

int_global_swnet_pres_abs <- int_global_swnet_raw %>% 
  mutate(interaction = 1)

# calculate  metrics ###########################################################
################################################################################.

nw.metrics.swnet <- list()
mod <- T

# "raw" metrics ----------------------------------------------------------------.
tmp <- mclapply(rownames(swnet_ct), f.nw.metrics,
                int_global_swnet_raw,
                mod = mod)
nw.metrics.swnet[["raw"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# correction: feeding type -----------------------------------------------------.
tmp <- mclapply(rownames(swnet_ct_fw), f.nw.metrics,
                int_global_swnet_fw,
                mod = mod)
nw.metrics.swnet[["fw"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# correction: taxonomic information --------------------------------------------.
tmp<- mclapply(rownames(swnet_ct), f.nw.metrics,
               int_global_swnet_pw,
               mod = mod)
nw.metrics.swnet[["pw"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# correction: metabolic rate ---------------------------------------------------.
tmp<- mclapply(rownames(swnet_ct_mr), f.nw.metrics,
               int_global_swnet_mr,
               mod = mod)
nw.metrics.swnet[["mr"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# correction: plant biomass ----------------------------------------------------.
tmp<- mclapply(rownames(swnet_ct), f.nw.metrics,
               int_global_swnet_bm,
               mod = mod)
nw.metrics.swnet[["bm"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# correction: feeding type, taxonomic information ------------------------------.
tmp <- mclapply(rownames(swnet_ct_fw), f.nw.metrics,
                int_global_swnet_fw_pw,
                mod = mod)
nw.metrics.swnet[["fw.pw"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# correction: feeding type, metabolic rate -------------------------------------.
tmp <- mclapply(rownames(swnet_ct_fw_mr), f.nw.metrics,
                int_global_swnet_fw_mr,
                mod = mod)
nw.metrics.swnet[["fw.mr"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# correction: feeding type, plant biomass --------------------------------------.
tmp <- mclapply(rownames(swnet_ct_fw), f.nw.metrics,
                int_global_swnet_fw_bm,
                mod = mod)
nw.metrics.swnet[["fw.bm"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# correction: taxonomic information, metabolic rate ----------------------------.
tmp <- mclapply(rownames(swnet_ct_mr), f.nw.metrics,
                int_global_swnet_pw_mr,
                mod = mod)
nw.metrics.swnet[["pw.mr"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# correction: taxonomic information, plant biomass -----------------------------.
tmp <- mclapply(rownames(swnet_ct), f.nw.metrics,
                int_global_swnet_pw_bm,
                mod = mod)
nw.metrics.swnet[["pw.bm"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# correction: metabolic rate, plant biomass ------------------------------------.
tmp <- mclapply(rownames(swnet_ct_mr), f.nw.metrics,
                int_global_swnet_mr_bm,
                mod = mod)
nw.metrics.swnet[["mr.bm"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# correction: feeding type, taxonomic information, metabolic rate --------------.
tmp <- mclapply(rownames(swnet_ct_fw_mr), f.nw.metrics,
                int_global_swnet_fw_pw_mr,
                mod = mod)
nw.metrics.swnet[["fw.pw.mr"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# correction: feeding type, taxonomic information, plant biomass ---------------.
tmp <- mclapply(rownames(swnet_ct_fw), f.nw.metrics,
                int_global_swnet_fw_pw_bm,
                mod = mod)
nw.metrics.swnet[["fw.pw.bm"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# correction: feeding type, metabolic rate, plant biomass ----------------------.
tmp <- mclapply(rownames(swnet_ct_fw_mr), f.nw.metrics,
                int_global_swnet_fw_mr_bm,
                mod = mod)
nw.metrics.swnet[["fw.mr.bm"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# correction: taxonomic information, metabolic rate, plant biomass -------------.
tmp <- mclapply(rownames(swnet_ct_mr), f.nw.metrics,
                int_global_swnet_pw_mr_bm,
                mod = mod)
nw.metrics.swnet[["pw.mr.bm"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# correction: feeding type, taxonomic information, metabolic rate, plant biomass
tmp <- mclapply(rownames(swnet_ct_fw_mr), f.nw.metrics,
                int_global_swnet_fw_pw_mr_bm,
                mod = mod)
nw.metrics.swnet[["fw.pw.mr.bm"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# no polyphagous species -------------------------------------------------------.
tmp <- mclapply(rownames(swnet_ct), f.nw.metrics,
                int_global_swnet_no_poly,
                mod = mod)
nw.metrics.swnet[["no.poly"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# all plants present -----------------------------------------------------------.
tmp <- mclapply(rownames(swnet_ct), f.nw.metrics,
                int_global_swnet_all_plants,
                mod = mod)
nw.metrics.swnet[["all.plants"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# PRESENCE - ABSENCE -----------------------------------------------------------.
tmp <- mclapply(rownames(swnet_ct), f.nw.metrics,
                int_global_swnet_pres_abs,
                mod = mod)
nw.metrics.swnet[["pres.abs"]] <- do.call(rbind.data.frame, tmp); rm(tmp)


# hierarchical models ##########################################################
################################################################################.

# define funvtion
hier_mod_swnet <- function(data){

  data <- data %>%
    left_join(lui.glob, by = "PlotID")
  
  out <- list()
  for (res in responses[responses %in% names(data)]){
    
    data$response <- (data[, res] - mean(data[, res]))/sd(data[, res])
    
    data.comb <- list(J = length(unique(data$Exploratory)),
                      N = nrow(data),
                      Exploratory = as.integer(as.factor(data$Exploratory)),
                      LUI_z = as.vector(data$LUI_z),
                      y = as.vector(data$response))
    
    
    fit <- sampling(object = stan_lui_mod, data = data.comb,
                    chains = 4, warmup = 500, iter = 5500, control = list(adapt_delta = 0.99))
    out_tmp <- extract(fit)[c("a", "b", "mu_a", "sigma_a", "sigma_b", "sigma_y")]
    pred_tmp <- extract(fit)$y_hat
    out_tmp[["resid"]] <- matrix(rep(data.comb$y, nrow(pred_tmp)),
                                 nrow = nrow(pred_tmp),
                                 byrow = T) - pred_tmp
    
    out[[res]] <- list(lui_comb = out_tmp)
    

    data.split <- list(J = length(unique(data$Exploratory)),
                       N = nrow(data),
                       Exploratory = as.integer(as.factor(data$Exploratory)),
                       M_z = as.vector(data$M_z),
                       F_z = as.vector(data$F_z),
                       G_z = as.vector(data$G_z),
                       y = as.vector(data$response))
    
    
    fit <- sampling(object = stan_luisplit_mod, data = data.split,
                    chains = 4, warmup = 500, iter = 5500, control = list(adapt_delta = 0.99))
    out_tmp <- extract(fit)[c("a", "b", "mu_a", "sigma_a", "sigma_b_m", 
                              "sigma_b_f", "sigma_b_g", "sigma_y")]
    pred_tmp <- extract(fit)$y_hat

    out[[res]] <- c(out[[res]], list(lui_split = out_tmp))
  }
  
  out
}

# apply function
stats_b_swnet <- mclapply(nw.metrics.swnet, hier_mod_swnet)
names(stats_b_swnet) <- names(nw.metrics.swnet)

# null models: on external cluster  ############################################
################################################################################.

# data export ------------------------------------------------------------------.

plots <- rownames(swnet_ct)

save(
     int_global_swnet_raw, 
     int_global_swnet_fw, int_global_swnet_pw,
     int_global_swnet_mr, int_global_swnet_bm, int_global_swnet_fw_pw,
     int_global_swnet_fw_mr, int_global_swnet_fw_bm, int_global_swnet_pw_mr,
     int_global_swnet_pw_bm, int_global_swnet_mr_bm, int_global_swnet_fw_pw_mr,
     int_global_swnet_fw_pw_bm, int_global_swnet_fw_mr_bm,
     int_global_swnet_pw_mr_bm, int_global_swnet_fw_pw_mr_bm,
     int_global_swnet_no_poly, int_global_swnet_all_plants,
     int_global_swnet_pres_abs,
     plots, 
     lui.glob, nw.metrics.swnet,
     file = "swnet_data.Rdata")


# data import ------------------------------------------------------------------.
# after cluster analysis

fit_swnet_null <- list()

loop.files <- list.files("dir_nullmodels",
                         full.names = T)
loop.files <- loop.files[grepl("swnet", loop.files)]
loop.names <- gsub("swnet_", "", list.files("dir_nullmodels",
                         full.names = F))
loop.names <- loop.names[!grepl("wintr", loop.names) & !grepl("system", loop.names)]

for (j in 1:length(loop.files)){
  tmp.files <- list.files(loop.files[j], full.names = T)
  if(length(tmp.files) > 0){
    tmp.names <- gsub(".txt", "", list.files(loop.files[j], full.names = F))
    
    # set values for the colClasses argument in read.csv. Improves speed (maybe? not too much)
    col.test <- read.table(tmp.files[1])
    cols <- ifelse(names(col.test) %in% c("response", "null.model", 
                                          "lui_comb.a.1", "lui_comb.a.2", "lui_comb.a.3",
                                          "lui_comb.b", 
                                          "lui_split.a.1", "lui_split.a.2", "lui_split.a.3",
                                          "lui_split.b.1","lui_split.b.2", 
                                          "lui_split.b.3"),
                   ifelse(names(col.test) %in% c("response", "null.model"), "factor", "numeric"), 
                   "NULL")
    names(cols) <- names(col.test)
    
    for (i in 1:length(tmp.files)){
      fit_swnet_null[[gsub("_", ".", loop.names[j])]][[tmp.names[i]]] <- read.table(tmp.files[i], colClasses = cols)
      
    }
  }
}


# plotting main figure #########################################################
################################################################################.

y_min_thresh <- 0.025

n.nullshades <- 100

for (datatype.def in c("raw", "no.poly", "pres.abs", "all.plants")){
  
  posterior_plots_swnet_def <- list()
  posterior_plots_swnet_int <- list()
  y_axes_swnet <- list()
  range_x_swnet <- list()
  
  for (res in responses){
    if (res %in% levels(fit_swnet_null[[datatype.def]]$data_00001$response)){
      
      # gather all data ----------------------------------------------------------.
      
      comb_nullmodel <- list()
      density_obs <- list()
      
      
      density_obs[["LUI"]] <- data.frame(density(stats_b_swnet[[datatype.def]][[res]][["lui_comb"]]$b)) %>% 
        filter(y > y_min_thresh)
      
      for (i in 1:length(fit_swnet_null[[datatype.def]])){
        
        comb_nullmodel[["LUI"]] <- c(comb_nullmodel[["LUI"]],
                                     filter(fit_swnet_null[[datatype.def]][[i]],
                                            response == res)$lui_comb.b)
      }
      
      mfg_names <- c("Mowing", "Fertilisation", "Grazing")
      for (mfg in 1:3){
        
        density_obs[[mfg_names[mfg]]] <- data.frame(density(stats_b_swnet[[datatype.def]][[res]][["lui_split"]]$b[, mfg])) %>% 
          filter(y > y_min_thresh)
        
        
        
        for (i in 1:length(fit_swnet_null[[datatype.def]])){
          
          comb_nullmodel[[mfg_names[mfg]]] <- c(comb_nullmodel[[mfg_names[mfg]]],
                                                filter(fit_swnet_null[[datatype.def]][[i]],
                                                       response == res)[ , paste0("lui_split.b.", mfg)])
          
        }
      }
      
      max_y <- max(c(sapply(comb_nullmodel, function(x) max(density(x)$y)),
                     sapply(density_obs, function(x) max(x$y))))
      
      
      # LUI comb -----------------------------------------------------------------.
      p <- ggplot(data = density_obs[["LUI"]], 
                  aes(x, y)) +
        xlab(res)
      
      
      density_null1 <- data.frame(density(comb_nullmodel[["LUI"]])) %>% 
        filter(y > y_min_thresh) 
      
      
      for (i in 1:n.nullshades){
        d_target <- data.frame(density(filter(fit_swnet_null[[datatype.def]][[i]],
                                              response == res)$lui_comb.b)) %>% 
          filter(y > y_min_thresh) %>%
          mutate(y = y / (.75 * max_y),
                 y = y + 1)
        
        p <-
          p +
          geom_ribbon(data = d_target, 
                      aes(ymin = 1, ymax = y),
                      fill = "#7570b3", alpha = .1) 
        
      }
      
      hdis_core <- rbind(nullmodel = hdi(comb_nullmodel[["LUI"]], ci = 1 - sqrt(0.05)), # argument: chance to be outside of both core areas (no overlap) is alpha * alpha
                         observed = hdi(stats_b_swnet[[datatype.def]][[res]][["lui_comb"]]$b, ci = 1 - sqrt(0.05))) %>% 
        rownames_to_column("data") %>% 
        mutate(data = factor(data, levels = .$data))
      
      hdis_outer <- rbind(nullmodel = hdi(comb_nullmodel[["LUI"]], ci = .95), 
                          observed = hdi(stats_b_swnet[[datatype.def]][[res]][["lui_comb"]]$b, ci = .95)) %>% 
        rownames_to_column("data") %>% 
        mutate(data = factor(data, levels = .$data))
      
      est <- rbind(nullmodel = map_estimate(comb_nullmodel[["LUI"]]), 
                   observed = map_estimate(stats_b_swnet[[datatype.def]][[res]][["lui_comb"]]$b)) %>% 
        rownames_to_column("data") %>% 
        mutate(data = factor(data, levels = .$data))
      
      p <-
        p +
        geom_vline(xintercept = 0, lty = "dashed") +
        geom_line(data = density_null1 %>% 
                    mutate(y = y / (.75 * max_y),
                           y = y + 1),
                  col = "white", size = 1.25) +
        geom_line(data = density_null1 %>% 
                    mutate(y = y / (.75 * max_y),
                           y = y + 1),
                  col = "#7570b3", size = .5) +
        geom_line(data = density_obs[["LUI"]] %>% 
                    mutate(y = y / (.75 * max_y),
                           y = y + 1),
                  col = "white", size = 1.25) +
        geom_line(data = density_obs[["LUI"]] %>% 
                    mutate(y = y / (.75 * max_y),
                           y = y + 1),
                  col = "#d95f02", size = .5) +
        geom_segment(data = hdis_outer, aes(x = CI_low, xend = CI_high, 
                                            y = as.numeric(data)/3 + .1, yend = as.numeric(data)/3 + .1,
                                            col = data), size = .3) +
        geom_segment(data = hdis_core, aes(x = CI_low, xend = CI_high, 
                                           y = as.numeric(data)/3 + .1, yend = as.numeric(data)/3 + .1,
                                           col = data), size = 1.25) +
        geom_point(data = est, aes(x = MAP, y = as.numeric(data)/3 + .1), 
                   col = "white", size = 2) +
        geom_point(data = est, aes(x = MAP, y = as.numeric(data)/3 + .1, 
                                   col = data), size = 1.25) +
        scale_color_manual(values = c(nullmodel = "#7570b3", 
                                      observed = "#d95f02"),
                           guide = F) + 
        theme(axis.title = element_blank(),
              axis.line.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank()) +
        scale_y_continuous(limits = c(0.2 , 3))

      posterior_plots_swnet_def[[res]][["LUI"]] <- p
      
      range_x_swnet[[res]]$min <- min(density_null1$x, 
                                      density_obs[["LUI"]]$x)
      range_x_swnet[[res]]$max <- max(density_null1$x, 
                                      density_obs[["LUI"]]$x)
      
      
      comb_nullmodel_int <- c()
      
      for (i in 1:length(fit_swnet_null[[datatype.def]])){
        
        comb_nullmodel_int <- c(comb_nullmodel_int,
                                rowMeans(filter(fit_swnet_null[[datatype.def]][[i]],
                                                response == res)[c("lui_comb.a.1", 
                                                                   "lui_comb.a.2", 
                                                                   "lui_comb.a.3")]))
      }
      
      
      
      sd_target <- sd(nw.metrics.swnet[[datatype.def]][[res]])
      m_target <- mean(nw.metrics.swnet[[datatype.def]][[res]])
      
      hdis_core_int <- rbind(nullmodel = hdi(comb_nullmodel_int * sd_target + m_target, 
                                             ci = 1 - sqrt(0.05)), # argument: chance to be outside of both core areas (no overlap) is alpha * alpha
                             observed = hdi(rowMeans(stats_b_swnet[[datatype.def]][[res]][["lui_comb"]]$a * sd_target + m_target), ci = 1 - sqrt(0.05))) %>% 
        rownames_to_column("data") %>% 
        mutate(data = factor(data, levels = .$data))
      
      hdis_outer_int <- rbind(nullmodel = hdi(comb_nullmodel_int * sd_target + m_target, ci = .95), 
                              observed = hdi(rowMeans(stats_b_swnet[[datatype.def]][[res]][["lui_comb"]]$a * sd_target + m_target), ci = .95)) %>% 
        rownames_to_column("data") %>% 
        mutate(data = factor(data, levels = .$data))
      
      est_int <- rbind(nullmodel = map_estimate(comb_nullmodel_int * sd_target + m_target), 
                       observed = map_estimate(rowMeans(stats_b_swnet[[datatype.def]][[res]][["lui_comb"]]$a * sd_target + m_target))) %>% 
        rownames_to_column("data") %>% 
        mutate(data = factor(data, levels = .$data))
      
      
      posterior_plots_swnet_int[[res]] <-
        ggplot() +
        geom_segment(data = hdis_outer_int, aes(y = CI_low, yend = CI_high, 
                                                x = data, xend = data,
                                                col = data), size = .3) +
        geom_segment(data = hdis_core_int, aes(y = CI_low, yend = CI_high, 
                                               x = data, xend = data,
                                               col = data), size = 1.25) +
        geom_point(data = est_int, aes(y = MAP, x = data), 
                   col = "white", size = 2) +
        geom_point(data = est_int, aes(y = MAP, x = data, 
                                       col = data), size = 1.25) +
        scale_color_manual(values = c(nullmodel = "#7570b3", 
                                      observed = "#d95f02"),
                           guide = F) + 
        theme(axis.title = element_blank(),
              axis.line.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank()) 
      
      
      # LUI split ----------------------------------------------------------------.
      
      for (mfg in 1:3){
        
        p <- ggplot(data = density_obs[[mfg_names[mfg]]], 
                    aes(x, y)) +
          xlab(res)
        
        
        
        density_null1 <- data.frame(density(comb_nullmodel[[mfg_names[mfg]]])) %>% 
          filter(y > y_min_thresh)
        
        
        for (i in 1:n.nullshades){
          
          p <-
            p +
            geom_ribbon(data = data.frame(density(filter(fit_swnet_null[[datatype.def]][[i]],
                                                         response == res)[ , paste0("lui_split.b.", mfg)])) %>% 
                          filter(y > y_min_thresh) %>% 
                          mutate(y = y / (.75 * max_y),
                                 y = y + 1), 
                        aes(ymin = 1, ymax = y),
                        fill = "#7570b3", alpha = .1) 
          
        }
        
        
        hdis_core <- rbind(nullmodel = hdi(comb_nullmodel[[mfg_names[mfg]]], ci = 1 - sqrt(0.05)), 
                           observed = hdi(stats_b_swnet[[datatype.def]][[res]][["lui_split"]]$b[, mfg], ci = 1 - sqrt(0.05))) %>% 
          rownames_to_column("data") %>% 
          mutate(data = factor(data, levels = .$data))
        
        hdis_outer <- rbind(nullmodel = hdi(comb_nullmodel[[mfg_names[mfg]]], ci = .95), 
                            observed = hdi(stats_b_swnet[[datatype.def]][[res]][["lui_split"]]$b[, mfg], ci = .95)) %>% 
          rownames_to_column("data") %>% 
          mutate(data = factor(data, levels = .$data))
        
        est <- rbind(nullmodel = map_estimate(comb_nullmodel[[mfg_names[mfg]]]), 
                     observed = map_estimate(stats_b_swnet[[datatype.def]][[res]][["lui_split"]]$b[, mfg])) %>% 
          rownames_to_column("data") %>% 
          mutate(data = factor(data, levels = .$data))
        
        p <-
          p +
          geom_vline(xintercept = 0, lty = "dashed") +
          geom_line(data = density_null1 %>% 
                      mutate(y = y / (.75 * max_y),
                             y = y + 1),
                    col = "white", size = 1.25) +
          geom_line(data = density_null1 %>% 
                      mutate(y = y / (.75 * max_y),
                             y = y + 1),
                    col = "#7570b3", size = .5) +
          geom_line(data = density_obs[[mfg_names[mfg]]] %>% 
                      mutate(y = y / (.75 * max_y),
                             y = y + 1),
                    col = "white", size = 1.25) +
          geom_line(data = density_obs[[mfg_names[mfg]]] %>% 
                      mutate(y = y / (.75 * max_y),
                             y = y + 1),
                    col = "#d95f02", size = .5) +
          geom_segment(data = hdis_outer, aes(x = CI_low, xend = CI_high, 
                                              y = as.numeric(data)/3 + .1, yend = as.numeric(data)/3 + .1,
                                              col = data), size = .3) +
          geom_segment(data = hdis_core, aes(x = CI_low, xend = CI_high, 
                                             y = as.numeric(data)/3 + .1, yend = as.numeric(data)/3 + .1,
                                             col = data), size = 1.25) +
          geom_point(data = est, aes(x = MAP, y = as.numeric(data)/3 + .1), 
                     col = "white", size = 2) +
          geom_point(data = est, aes(x = MAP, y = as.numeric(data)/3 + .1, 
                                     col = data), size = 1.25) +
          scale_color_manual(values = c(nullmodel = "#7570b3", 
                                        observed = "#d95f02"),
                             guide = F) + 
          theme(axis.title = element_blank(),
                axis.line.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank()) +
          scale_y_continuous(limits = c(0.2 , 3))

        posterior_plots_swnet_def[[res]][[mfg_names[mfg]]] <- p
        
        range_x_swnet[[res]]$min <- c(range_x_swnet[[res]]$min, min(density_null1$x, 
                                                                    density_obs[[mfg_names[mfg]]]$x))
        range_x_swnet[[res]]$max <- c(range_x_swnet[[res]]$max, max(density_null1$x, 
                                                                    density_obs[[mfg_names[mfg]]]$x))
        
      }
      
      
      # y axes plots -------------------------------------------------------------.
      
      y_axes_swnet[[res]] <-
        ggplot() + 
        scale_y_continuous(limits = c(0.2 , 3),
                           breaks = c(0, floor(max_y)) / (0.75 * max_y) + 1,
                           labels = c(0, floor(max_y)),
                           position = "right") +
        scale_x_continuous(limits = c(0, 0)) +
        annotate(geom = "segment", x = Inf, xend = Inf, y = 1, yend = 1 / (0.75) + 1) +
        theme(axis.line.y = element_blank(),
              axis.line.x = element_blank(),
              axis.text.x = element_blank(),
              axis.title = element_blank(),
              axis.ticks.x = element_blank(),
              plot.margin = margin(r = 0))
    }
  }
  
  axis.text.size <- 8
  axis.title.size <- 10
  label.size <- 12
  
  posterior_plots_swnet_int <- lapply(posterior_plots_swnet_int, function(x) x + 
                                        theme(axis.text = element_text(size = axis.text.size))+
                                        scale_y_continuous(position = "right")) 
  
  y_axes_swnet <- lapply(y_axes_swnet, function(x) x + 
                           theme(axis.text = element_text(size = axis.text.size)))
  
  
  legend <-
    ggplot(data.frame(var = factor(c("observed", "null model"), 
                                   levels = c("observed", "null model")),
                      x = NA, y = NA)) +
    geom_area(aes(fill = var, x = x, y = y)) +
    scale_fill_manual(values = c(observed = "#d95f02",
                                 "null model" = "#7570b3"))  +
    theme_nothing() +
    theme(legend.position = c(.5, .5),
          legend.justification = c(.5, .5),
          legend.title = element_blank(),
          legend.background = element_rect(fill = col.strip,
                                           colour = NA),
          legend.text = element_text(size = axis.title.size),
          legend.key.height = unit(2, "mm"))
  
  empty <- ggplot() + geom_blank()
  
  
  range_x <- c(min(c(range_x_swnet$nested$min, range_x_swnet$qmod$min, range_x_swnet$rob.A$min)),
               max(c(range_x_swnet$nested$max, range_x_swnet$qmod$max, range_x_swnet$rob.A$max)))
  
  
  rel_heights <- c(.22, .22, .9, .9, .9, .32, .18)
  rel_heights <- rel_heights / sum(rel_heights)
  
  rel_widths <- c(.26, 1, 1, 1, 1) # metric label, 4 figures (land-use components)
  rel_widths <- rel_widths / sum(rel_widths)
  rel_widths_sub <- c(sum(rel_widths[2:5]), rel_widths[2] / 3, rel_widths[2] / 4.1) # core figure, probability density axis, probability density title
  
  rel_widths_overall <- c(4.7, .2, 1)
  
  
  col.strip <- "grey80"
  
  
  p <- plot_grid(
    # Slope strip:
    plot_grid(plot_grid(ggplot_labeller(expression(bold("a")), 0, label.size, NA),
                        ggplot_labeller("Slope", 0, label.size, col.strip),
                        nrow = 1, rel_widths = c(rel_widths[1], sum(rel_widths_sub))),
              
              plot_grid(
                # core figure (slope):
                plot_grid(plot_grid(empty,
                                    ggplot_labeller("LUI", 0, label.size, col.strip) +
                                      theme(plot.background = element_rect(fill = "transparent", colour = NA, size = NA)),
                                    ggplot_labeller("Mowing", 0, label.size, col.strip) +
                                      theme(plot.background = element_rect(fill = "transparent", colour = NA, size = NA)),
                                    ggplot_labeller("Grazing", 0, label.size, col.strip) +
                                      theme(plot.background = element_rect(fill = "transparent", colour = NA, size = NA)),
                                    ggplot_labeller("Fertilisation", 0, label.size, col.strip) +
                                      theme(plot.background = element_rect(fill = "transparent", colour = NA, size = NA)),
                                    nrow = 1, rel_widths = rel_widths),
                          plot_grid(plotlist = c(list(ggplot_labeller("Modularity", 90, label.size, col.strip) +
                                                        theme(plot.background = element_rect(fill = "transparent", colour = NA, size = NA))),
                                                 lapply(posterior_plots_swnet_def$qmod,
                                                        function(x) x +
                                                          scale_x_continuous(limits = range_x * 1.05) +
                                                          theme_nothing())[c(1,2,4,3)]),
                                    nrow = 1, align = "h", rel_widths = rel_widths),
                          plot_grid(plotlist = c(list(ggplot_labeller("Nestedness", 90, label.size, col.strip) +
                                                        theme(plot.background = element_rect(fill = "transparent", colour = NA, size = NA))),
                                                 lapply(posterior_plots_swnet_def$nested,
                                                        function(x) x +
                                                          scale_x_continuous(limits = range_x * 1.05) +
                                                          theme_nothing())[c(1,2,4,3)]),
                                    nrow = 1, align = "h", rel_widths = rel_widths),
                          plot_grid(plotlist = c(list(ggplot_labeller("Robustness", 90, label.size, col.strip) +
                                                        theme(plot.background = element_rect(fill = "transparent", colour = NA, size = NA))),
                                                 lapply(posterior_plots_swnet_def$rob.A,
                                                        function(x) x +
                                                          scale_x_continuous(limits = range_x * 1.05) +
                                                          theme_nothing())[c(1,2,4,3)]),
                                    nrow = 1, align = "h", rel_widths = rel_widths),
                          plot_grid(plotlist =  c(list(empty),
                                                  rep(list(posterior_plots_swnet_def$rob.A$LUI +
                                                             scale_x_continuous(limits = range_x * 1.05) +
                                                             scale_y_continuous(limits = c(0, 0)) +
                                                             theme(axis.text.x = element_text(size = axis.text.size,
                                                                                              angle = 90, hjust = 1, vjust = .5),
                                                                   plot.background = element_rect(fill = "transparent", colour = NA, size = NA),
                                                                   panel.background = element_rect(fill = "transparent", colour = NA, size = NA))), 4)),
                                    nrow = 1, rel_widths = rel_widths),
                          plot_grid(plotlist =  list(empty,
                                                     ggplot_labeller("Standardised model coefficient", 0, axis.title.size, NA)),
                                    nrow = 1, rel_widths = c(rel_widths[1], sum(rel_widths[2:5]))),
                          ncol = 1, rel_heights = rel_heights[-1]),
                
                # y-axis density:
                plot_grid(empty, 
                          y_axes_swnet$qmod,
                          y_axes_swnet$nested,
                          y_axes_swnet$rob.A,
                          empty, empty,
                          ncol = 1, rel_heights = rel_heights[-1], align = "v"),
                
                # y-axis density label:
                plot_grid(empty,
                          ggplot_labeller("Probability density", -90, axis.title.size, NA) +
                            theme(plot.margin = margin(r = 4)),
                          empty, empty,
                          ncol = 1,
                          rel_heights = c(rel_heights[2], sum(rel_heights[3:5]), 
                                          rel_heights[6:7])),
                
                nrow = 1, rel_widths = rel_widths_sub),
              ncol = 1, rel_heights = c(rel_heights[1], sum(rel_heights[-1]))
    ),
    
    #-------------------------------------------------------------------------------.
    ggplot() + 
      annotate(geom = "segment", x = 1, xend = 1, y = -.8, yend = Inf, lty = 3) +
      scale_y_continuous(limits = c(-1, 1)) +
      theme_nothing(),
    #-------------------------------------------------------------------------------.
    
    plot_grid(
      plot_grid(ggplot_labeller(expression(bold("b")), 0, label.size, NA),
                ggplot_labeller("Intercept", 0, label.size, col.strip),
                nrow = 1, rel_widths = c(1, 5 * rel_widths_overall[3])),
      plot_grid(plot_grid(plot_grid(empty,
                                    posterior_plots_swnet_int$qmod,
                                    posterior_plots_swnet_int$nested,
                                    posterior_plots_swnet_int$rob.A,
                                    ncol = 1, align = "v", rel_heights = rel_heights[2:5]),
                          plot_grid(empty,
                                    ggplot_labeller("Rescaled model coefficient", -90, axis.title.size, NA),
                                    rel_heights = c(rel_heights[2], sum(rel_heights[3:5])),
                                    ncol = 1),
                          nrow = 1, rel_widths = c(5, 1)),
                legend, 
                ncol = 1,
                rel_heights = c(sum(rel_heights[3:5]), sum(rel_heights[6:7]))),
      ncol = 1, rel_heights = c(rel_heights[1], sum(rel_heights[-1]))
    ),
    rel_widths = rel_widths_overall, nrow = 1)
  
}

# plotting (different corrections) #############################################
################################################################################.

y_min_thresh <- 0.01

posterior_plots_swnet_corr <- list()
range_x_swnet_corr <- list()

corrections <- c("fw", "pw", "mr", "bm",
                 "fw.pw", "fw.mr", "fw.bm", "pw.mr", "pw.bm", "mr.bm",
                 "fw.pw.mr", "fw.pw.bm", "fw.mr.bm", "pw.mr.bm",
                 "fw.pw.mr.bm", "raw")

for (datatype_i in corrections[corrections %in% names(stats_b_swnet)]){
  for (res in responses){
    if (res %in% levels(fit_swnet_null[[datatype_i]]$data_00001$response)){
      range_x_swnet_corr[[datatype_i]][[res]] <- c(0, 0)
      
      # LUI comb -----------------------------------------------------------------.
      density_obs <- data.frame(density(stats_b_swnet[[datatype_i]][[res]][["lui_comb"]]$b)) %>% 
        filter(y > y_min_thresh)
      p <- ggplot(data = density_obs, 
                  aes(x, y)) +
        xlab(res)
      
      comb_nullmodel <- c()

      for (i in 1:length(fit_swnet_null[[datatype_i]])){
        
        comb_nullmodel <- c(comb_nullmodel,
                            filter(fit_swnet_null[[datatype_i]][[i]],
                                   response == res)$lui_comb.b)
      }
      
      
      
      max_y <- max(density(comb_nullmodel)$y, 
                   p$data$y)
      
      density_null <- data.frame(density(comb_nullmodel)) %>% 
        filter(y > y_min_thresh) 
      
      
      
      hdis_core <- rbind(nullmodel = hdi(comb_nullmodel, ci = 1 - sqrt(0.05)), # argument: chance to be outside of both core areas (no overlap) is alpha * alpha
                         observed = hdi(stats_b_swnet[[datatype_i]][[res]][["lui_comb"]]$b, ci = 1 - sqrt(0.05))) %>% 
        rownames_to_column("data") %>% 
        mutate(data = factor(data, levels = .$data))
      
      hdis_outer <- rbind(nullmodel = hdi(comb_nullmodel, ci = .95), 
                          observed = hdi(stats_b_swnet[[datatype_i]][[res]][["lui_comb"]]$b, ci = .95)) %>% 
        rownames_to_column("data") %>% 
        mutate(data = factor(data, levels = .$data))
      
      est <- rbind(nullmodel = map_estimate(comb_nullmodel), 
                   observed = map_estimate(stats_b_swnet[[datatype_i]][[res]][["lui_comb"]]$b)) %>% 
        rownames_to_column("data") %>% 
        mutate(data = factor(data, levels = .$data))
      
      p <-
        p +
        geom_vline(xintercept = 0, lty = "dashed") +
        geom_line(data = density_null %>% 
                    mutate(y = y / (.75 * max_y),
                           y = y + 1),
                  col = "#7570b3", size = .75) +
        geom_line(data = density_obs %>% 
                    mutate(y = y / (.75 * max_y),
                           y = y + 1),
                  col = "#d95f02", size = .75) +
        geom_segment(data = hdis_outer, aes(x = CI_low, xend = CI_high, 
                                            y = as.numeric(data)/3 + .1, yend = as.numeric(data)/3 + .1,
                                            col = data), size = .3) +
        geom_segment(data = hdis_core, aes(x = CI_low, xend = CI_high, 
                                           y = as.numeric(data)/3 + .1, yend = as.numeric(data)/3 + .1,
                                           col = data), size = 1.25) +
        geom_point(data = est, aes(x = MAP, y = as.numeric(data)/3 + .1, 
                                   col = data), size = 1.75) +
        scale_color_manual(values = c(nullmodel = "#7570b3", 
                                      observed = "#d95f02"),
                           guide = F) + 
        theme(axis.title = element_blank(),
              axis.line.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank()) +
        scale_y_continuous(limits = c(0.2 , 3))
      # scale_x_continuous(limits = c(min_x + xrange_off, 
      #                               max_x - xrange_off))
      
      posterior_plots_swnet_corr[[res]][["LUI"]][[datatype_i]] <- p
      
      range_x_swnet_corr[[datatype_i]][[res]][1] <- min(density_null$x,
                                                        density_obs$x, 
                                                        range_x_swnet_corr[[datatype_i]][[res]][1])
      range_x_swnet_corr[[datatype_i]][[res]][2] <- max(density_null$x,
                                                        density_obs$x,
                                                        range_x_swnet_corr[[datatype_i]][[res]][2])
      
      
      # LUI split ----------------------------------------------------------------.
      mfg_names <- c("Mowing", "Fertilisation", "Grazing")
      for (mfg in 1:3){
        
        density_obs <- data.frame(density(stats_b_swnet[[datatype_i]][[res]][["lui_split"]]$b[, mfg])) %>% 
          filter(y > y_min_thresh)
        p <- ggplot(data = density_obs, 
                    aes(x, y)) +
          xlab(res)
        
        
        comb_nullmodel <- c()
        
        
        for (i in 1:length(fit_swnet_null[[datatype_i]])){
          
          comb_nullmodel <- c(comb_nullmodel,
                              filter(fit_swnet_null[[datatype_i]][[i]],
                                     response == res)[ , paste0("lui_split.b.", mfg)])
          
        }
        
        max_y <- max(density(comb_nullmodel)$y, 
                     p$data$y)
        
        density_null <- data.frame(density(comb_nullmodel)) %>% 
          filter(y > y_min_thresh)
        
        
        
        hdis_core <- rbind(nullmodel = hdi(comb_nullmodel, ci = 1 - sqrt(0.05)), 
                           observed = hdi(stats_b_swnet[[datatype_i]][[res]][["lui_split"]]$b[, mfg], ci = 1 - sqrt(0.05))) %>% 
          rownames_to_column("data") %>% 
          mutate(data = factor(data, levels = .$data))
        
        hdis_outer <- rbind(nullmodel = hdi(comb_nullmodel, ci = .95), 
                            observed = hdi(stats_b_swnet[[datatype_i]][[res]][["lui_split"]]$b[, mfg], ci = .95)) %>% 
          rownames_to_column("data") %>% 
          mutate(data = factor(data, levels = .$data))
        
        est <- rbind(nullmodel = map_estimate(comb_nullmodel), 
                     observed = map_estimate(stats_b_swnet[[datatype_i]][[res]][["lui_split"]]$b[, mfg])) %>% 
          rownames_to_column("data") %>% 
          mutate(data = factor(data, levels = .$data))
        
        p <-
          p +
          geom_vline(xintercept = 0, lty = "dashed") +
          geom_line(data = density_null %>% 
                      mutate(y = y / (.75 * max_y),
                             y = y + 1),
                    col = "#7570b3", size = .75) +
          geom_line(data = density_obs %>% 
                      mutate(y = y / (.75 * max_y),
                             y = y + 1),
                    col = "#d95f02", size = .75) +
          geom_segment(data = hdis_outer, aes(x = CI_low, xend = CI_high, 
                                              y = as.numeric(data)/3 + .1, yend = as.numeric(data)/3 + .1,
                                              col = data), size = .3) +
          geom_segment(data = hdis_core, aes(x = CI_low, xend = CI_high, 
                                             y = as.numeric(data)/3 + .1, yend = as.numeric(data)/3 + .1,
                                             col = data), size = 1.25) +
          geom_point(data = est, aes(x = MAP, y = as.numeric(data)/3 + .1, 
                                     col = data), size = 1.75) +
          scale_color_manual(values = c(nullmodel = "#7570b3", 
                                        observed = "#d95f02"),
                             guide = F) + 
          theme(axis.title = element_blank(),
                axis.line.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank()) +
          scale_y_continuous(limits = c(0.2 , 3))

        posterior_plots_swnet_corr[[res]][[mfg_names[mfg]]][[datatype_i]] <- p
        
        range_x_swnet_corr[[datatype_i]][[res]][1] <- min(density_null$x,
                                                          density_obs$x, 
                                                          range_x_swnet_corr[[datatype_i]][[res]][1])
        range_x_swnet_corr[[datatype_i]][[res]][2] <- max(density_null$x,
                                                          density_obs$x,
                                                          range_x_swnet_corr[[datatype_i]][[res]][2])
        
      }
      
    }
  }
}


res <- "qmod"
var <- "LUI"


range_x_swnet_corr_sel <- list()
for (res in responses){
  range_x_swnet_corr_sel[[res]] <-
    c(min(sapply(range_x_swnet_corr[!names(range_x_swnet_corr) %in% 
                                      c("all.plants", "no.poly", "pres.abs")],
                 function(x) x[[res]][1])),
      max(sapply(range_x_swnet_corr[!names(range_x_swnet_corr) %in% 
                                      c("all.plants", "no.poly", "pres.abs")],
                 function(x) x[[res]][2])))
}

f.plotadjust <- function(x){
  x + theme(plot.background = element_rect(fill = "transparent", colour = NA),
            panel.background = element_rect(fill = "transparent", colour = NA),
            axis.text.x = element_blank(),
            plot.margin = margin(0, 0, 0, 0)) +
    scale_x_continuous(limits = range_x_swnet_corr_sel[[res]])
}

swnet_corrplots <- list()
for (res in responses){
  for (var in c("LUI", "Mowing", "Grazing", "Fertilisation")){
    label <- var
    
    row1 <- plot_grid(plotlist = lapply(list(ggplot() + ggtitle(label) + 
                                               theme(plot.title = element_text(size = 25),
                                                     axis.line = element_blank(),
                                                     axis.ticks = element_blank(),
                                                     axis.text = element_blank()), 
                                             NULL, 
                                             posterior_plots_swnet_corr[[res]][[var]]$fw, 
                                             NULL, 
                                             NULL),
                                        f.plotadjust),
                      nrow = 1)
    row2 <- plot_grid(plotlist = lapply(list(NULL, 
                                             posterior_plots_swnet_corr[[res]][[var]]$fw.bm, 
                                             posterior_plots_swnet_corr[[res]][[var]]$fw.pw.bm,  
                                             posterior_plots_swnet_corr[[res]][[var]]$fw.pw,  
                                             posterior_plots_swnet_corr[[res]][[var]]$fw.mr),
                                        f.plotadjust),
                      nrow = 1)
    row3 <- plot_grid(plotlist = lapply(list(posterior_plots_swnet_corr[[res]][[var]]$bm,
                                             posterior_plots_swnet_corr[[res]][[var]]$fw.mr.bm,
                                             posterior_plots_swnet_corr[[res]][[var]]$fw.pw.mr.bm, 
                                             posterior_plots_swnet_corr[[res]][[var]]$fw.pw.mr,
                                             posterior_plots_swnet_corr[[res]][[var]]$pw),
                                        f.plotadjust),
                      nrow = 1)
    row4 <- plot_grid(plotlist = lapply(list(NULL,
                                             posterior_plots_swnet_corr[[res]][[var]]$mr.bm,
                                             posterior_plots_swnet_corr[[res]][[var]]$pw.mr.bm,
                                             posterior_plots_swnet_corr[[res]][[var]]$pw.mr,
                                             NULL),
                                        f.plotadjust),
                      nrow = 1)
    row5 <- plot_grid(plotlist = lapply(list(posterior_plots_swnet_corr[[res]][[var]]$raw,
                                             NULL,
                                             posterior_plots_swnet_corr[[res]][[var]]$mr,
                                             posterior_plots_swnet_corr[[res]][[var]]$pw.bm,
                                             NULL), 
                                        f.plotadjust),
                      nrow = 1)
    row6 <- plot_grid(NULL, NULL,
                      ggplot() + geom_blank() +
                        scale_x_continuous(limits = range_x_swnet_corr_sel[[res]]) +
                        theme(axis.line.x = element_blank(),
                              axis.ticks.x = element_blank(),
                              axis.text.x = element_text(angle = 90, size = 16,
                                                         hjust = 1, vjust = .5),
                              plot.background = element_rect(fill = "transparent", colour = NA),
                              panel.background = element_rect(fill = "transparent", colour = NA),
                              plot.margin = margin(0, 0, 0, 0)), 
                      NULL, NULL,
                      nrow = 1)
    
    swnet_corrplots[[res]][[var]] <- plot_grid(row1, row2, row3, row4, row5, row6, ncol = 1,
                      rel_heights = c(1, 1, 1, 1, 1, .5))
  }
}

separator <- ggplot() + 
  annotate(geom = "segment", x = 1, xend = 1, y = -Inf, yend = Inf, lty = 3) +
  scale_y_continuous(limits = c(-1, 1)) +
  theme_nothing()

for (res in responses){
  p <- plot_grid(
    plot_grid(swnet_corrplots[[res]]$LUI, separator, swnet_corrplots[[res]]$Mowing,
              nrow = 1, rel_widths = c(20, 1, 20)),
    
    ggplot() + 
      annotate(geom = "segment", x = Inf, xend = -Inf, y = 1, yend = 1, lty = 3) +
      scale_x_continuous(limits = c(-1, 1)) +
      theme_nothing(),
    
    plot_grid(swnet_corrplots[[res]]$Grazing, separator, swnet_corrplots[[res]]$Fertilisation,
              nrow = 1, rel_widths = c(20, 1, 20)),
    
    ncol = 1, rel_heights = c(20, 1, 20)
  )
}


# investigate the modularity pattern ###########################################
################################################################################.

extract.modules <- function(pl_i){
  
  int_matrix_plot <- int_global_swnet_raw %>% 
    filter(PlotID == pl_i) %>% 
    select(-PlotID) %>% 
    spread(Herb_Name_std, interaction, fill = 0) %>% 
    column_to_rownames("Host_Name_std") %>% 
    as.data.frame()
  
  
  modules.target <- computeModules(int_matrix_plot)@modules
  
  colnames(modules.target) <- c("", "", rownames(int_matrix_plot), colnames(int_matrix_plot))
  rownames(modules.target) <- c("", 1:(nrow(modules.target) - 1))
  modules.target <- modules.target[-1, -c(1:2)]
  for (i in 1:nrow(modules.target)){
    modules.target[i, ][modules.target[i, ] > 0] <- i
  }
  out <- data.frame(species = colnames(modules.target),
                    module = colSums(modules.target),
                    stringsAsFactors = F)
  rownames(out) <- NULL
  out <- out %>% 
    column_to_rownames("species") %>% 
    t
  
  
  out <- as.data.frame(out)
  out$PlotID <- pl_i
  out %>% 
    select(PlotID, everything())
  
}

d.modules <- mclapply(rownames(swnet_ct), extract.modules)
d.modules <- do.call(bind_rows, d.modules) %>%
  column_to_rownames("PlotID")


f.modularity.area <- function(pl_i) {
  
  d.modules.target <-
    d.modules[pl_i, ] %>% 
    gather(species, module) %>% 
    filter(!is.na(module))
  
  int_global_swnet_raw %>% 
    filter(PlotID == pl_i) %>% 
    select(-PlotID) %>% 
    left_join(d.modules.target %>%
                rename(module.host = "module"),
              by = c("Host_Name_std" = "species")) %>% 
    left_join(d.modules.target %>%
                rename(module.herb = "module"),
              by = c("Herb_Name_std" = "species")) %>% 
    mutate(in.out = ifelse(module.host == module.herb, paste0("in.", module.host), "out"),
           area.tot = length(unique(Host_Name_std)) * length(unique(Herb_Name_std))) %>% 
    left_join(plant_agg_g %>% select(-c(PlotID, Cover)) %>% distinct(),
              by = "Host_Name_std") %>% 
    group_by(in.out) %>% 
    mutate(max = max(table(Host_Order)/sum(table(Host_Order))),
           what = ifelse(max > .5,  # if > 50 % of interactions in a module are belonging to one Order
                         names(table(Host_Order)[max(table(Host_Order)) == table(Host_Order)]),
                         "diverse")) %>% 
    summarise(area = length(unique(Host_Name_std)) * length(unique(Herb_Name_std)),
              area.tot = unique(area.tot),
              what = unique(what)) %>% 
    ungroup() %>% 
    mutate(area = ifelse(in.out == "out", 0, area),
           area = ifelse(in.out == "out", area.tot - sum(area), area),
           what = ifelse(in.out == "out", "out", what)) %>% 
    group_by(what) %>% 
    summarise(area.prop = sum(area) / unique(area.tot)) %>% 
    ungroup() %>% 
    mutate(PlotID = pl_i)
}

mod.area <- lapply(rownames(swnet_ct), f.modularity.area)

mod.area <- do.call(rbind, mod.area)


mod.area %>% 
  mutate(what = ifelse(what %in% c("diverse", "Poales", "out"), what, "rest")) %>% 
  group_by(PlotID, what) %>% 
  summarise(area.prop = sum(area.prop)) %>% 
  ungroup() %>% 
  left_join(lui.glob, by = "PlotID") %>%
  select(M_std, F_std, G_std, what, area.prop) %>% 
  gather(var, int, -c(what, area.prop)) %>% 
  mutate(var = factor(var, levels = c("M_std", "F_std", "G_std"))) %>% 
  ggplot(aes(x = int, y = area.prop, col = what)) +
  geom_point() + 
  stat_smooth(method = "lm") +
  facet_grid(~ var, scales = "free_x") 

mod.area.in <-
  mod.area %>% 
  filter(what != "out") %>% 
  group_by(PlotID) %>% 
  mutate(area.prop.in = area.prop / sum(area.prop)) %>% 
  mutate(what = ifelse(what %in% c("diverse", "Poales"), what, "rest")) %>% 
  group_by(PlotID, what) %>% 
  summarise(area.prop.in = sum(area.prop.in)) %>% 
  ungroup() %>% 
  spread(what, area.prop.in, fill = 0) %>% 
  gather(what, area.prop.in, -PlotID)

# models:

lm.area <- mod.area.in %>% 
  left_join(lui.glob, by = "PlotID") %>%
  group_split(what) %>% 
  setNames(c("diverse", "Poales", "rest")) %>%
  map(~ lmer(asin(sqrt(area.prop.in)) ~ M_z + G_z + F_z  + (1 | Exploratory), data = .))


lapply(lm.area, modelplots, expl.plots = F) # OK for Poales (not for diverse, rest!)
lapply(lm.area, summary)

pval <- lapply(lm.area, function(x) summary(x)$coefficients[c("M_z", "G_z", "F_z"), 5]) %>% 
  do.call(rbind, .) %>% 
  as.data.frame() %>% 
  rownames_to_column("what") %>% 
  gather(var, pval, -what) %>% 
  mutate(var = gsub("M_z", "Mowing", var),
         var = gsub("G_z", "Grazing", var),
         var = gsub("F_z", "Fertilisation", var))

# plotting:

colors.plantgroups <- c(Poales = "#A3ED5A",
                        diverse = "#EB4B98",
                        rest = "#383736") 

labels.plantgroups <- c(Poales = "Grasses",
                        diverse = "Diverse",
                        rest = "Rest") 

var.labels <- c(M_std = "Mowing", G_std = "Grazing", F_std = "Fertilisation")

mod.area.in %>% 
  left_join(lui.glob, by = "PlotID") %>%
  select(M_std, F_std, G_std, what, area.prop.in) %>% 
  gather(var, int, -c(what, area.prop.in)) %>% 
  mutate(var = var.labels[var]) %>% 
  left_join(pval, by = c("what", "var")) %>% 
  mutate(what = factor(what, levels = names(labels.plantgroups)),
         var = factor(var, levels = var.labels),
         sign = ifelse(pval > .05, "ns", "sign")) %>% 
  ggplot(aes(x = int, y = area.prop.in, col = what)) +
  geom_point(size = .5) + 
  stat_smooth(aes(lty = sign), method = "lm", se = F) +
  facet_grid(~ var, scales = "free_x") +
  scale_color_manual(values = colors.plantgroups, labels = labels.plantgroups,
                     name = "Group") +
  scale_linetype_manual(values = c(ns = 2, sign = 1), guide = F) +
  xlab("Intensity") +
  ylab("Module area occupied") +
  scale_y_continuous(breaks = c(0, .25, .5, .75, 1),
                     labels = paste0(c(0, 25, 50, 75, 100), "%"))


# compare interactions in modules to interactions outside modules --------------.

f.modularity.density.in.out <- function(pl_i) {
  
  d.modules.target <- d.modules[pl_i, ] %>%
    gather(species, module) %>%
    filter(!is.na(module))
  
  int_global_swnet_raw %>% 
    filter(PlotID == pl_i) %>% 
    select(-PlotID) %>% 
    left_join(d.modules.target %>%
                rename(module.host = "module"),
              by = c("Host_Name_std" = "species")) %>% 
    left_join(d.modules.target %>%
                rename(module.herb = "module"),
              by = c("Herb_Name_std" = "species")) %>% 
    mutate(in.out = ifelse(module.host == module.herb, "in", "out"),
           module = module.herb,
           area.tot = length(unique(Host_Name_std)) * length(unique(Herb_Name_std)),
           host.length.tot = length(unique(Host_Name_std)),
           herb.length.tot = length(unique(Herb_Name_std))) %>% 
    left_join(plant_agg_g %>% select(-c(PlotID, Cover)) %>% distinct(),
              by = "Host_Name_std") %>% 
    group_by(in.out, module) %>% 
    mutate(max = max(table(Host_Order)/sum(table(Host_Order))),
           what = ifelse(max > .5, 
                         names(table(Host_Order)[max(table(Host_Order)) == table(Host_Order)]),
                         "diverse")) %>% 
    summarise(area = length(unique(Host_Name_std)) * length(unique(Herb_Name_std)),
              area.int = n(),
              area.tot = unique(area.tot),
              what = unique(what),
              area.out = length(unique(Herb_Name_std)) * unique(host.length.tot) - area) %>% 
    group_by(module) %>% 
    summarise(area.in = area[in.out == "in"],
              area.in.int = area.int[in.out == "in"],
              area.out = area.out[in.out == "in"],
              area.out.int = ifelse(is.null(area.int[in.out == "out"]), 0, 
                                    area.int[in.out == "out"]),
              what = what[in.out == "in"]) %>% 
    ungroup() %>% 
    mutate(area.out.int = ifelse(is.na(area.out.int), 0, area.out.int)) %>% 
    group_by(what) %>% 
    summarise(dens.in = sum(area.in.int) / sum(area.in),
              dens.out = sum(area.out.int) / sum(area.out),
              size = sum(area.in + area.out)) %>% 
    ungroup() %>% 
    mutate(PlotID = pl_i)
  
  
}

d.density <- mclapply(rownames(swnet_ct), f.modularity.density.in.out)
d.density <- do.call(rbind, d.density)

d.density %>% 
  mutate(diff.in.out = dens.in - dens.out,
         what = ifelse(what %in% c("diverse", "Poales"), what, "rest"),
         what = factor(what, levels = names(labels.plantgroups))) %>%
  ggplot(aes(x = what, y = diff.in.out, col = what)) +
  geom_point(aes(size = size), alpha = .5, position = position_jitter(.3)) +
  geom_boxplot(fill = "grey80", alpha = .5) +
  scale_color_manual(values = colors.plantgroups, guide = F) +
  theme(axis.title.x = element_blank()) +
  scale_size_continuous(range = c(.25, 3), name = "Module size") +
  scale_x_discrete(labels = labels.plantgroups) +
  ylab("Connectance difference")


# robustness & specialization along grazing gradient ###########################
################################################################################.

p1 <- nw.metrics.swnet$raw %>% 
  left_join(lui.glob) %>% 
  ggplot(aes(x = G_std, y = rob.A)) +
  geom_point(size = .5) + 
  stat_smooth(method = "loess", col = "olivedrab",
              fill = "olivedrab") +
  ylab("Robustness") +
  xlab("Grazing intensity")



# sharp increase in robustness with little grazing, then again a bit of a drop


# specialisation along gradient

spec_cat_g <- int_global_swnet_raw %>% 
  # filter(Herb_Name_std == "Acalypta parvula") %>% View
  left_join(plant_agg_g %>% select(-c(PlotID, Cover)) %>% distinct(),
            by = "Host_Name_std") %>% 
  group_by(Herb_Name_std) %>% 
  summarise(n_spec = length(unique(Host_Species)),
            n_gen = length(unique(Host_Genus)),
            n_fam = length(unique(Host_Family))) %>% 
  ungroup() %>% 
  mutate(specialisation = ifelse(n_spec == 1 & n_gen == 1, "mono_spec",
                                 ifelse(n_gen == 1, "mono_gen",
                                        ifelse(n_fam == 1, "mono_fam",
                                               "poly"))))

spec_spec <- swnet_agg %>% 
  inner_join(spec_cat_g, by = "Herb_Name_std") %>% 
  group_by(PlotID, specialisation) %>% 
  summarise(NumberAdults = sum(NumberAdults),
            NumberSpecies = length(unique(Herb_Name_std))) %>% 
  ungroup() %>% 
  group_by(PlotID) %>% 
  mutate(spec_perc_abund = NumberAdults / sum(NumberAdults) * 100,
         spec_perc_rich = NumberSpecies / sum(NumberSpecies) * 100) %>% 
  ungroup() %>% 
  filter(specialisation == "mono_spec")

p2 <- spec_spec %>% 
  left_join(lui.glob, by = "PlotID") %>% 
  ggplot(aes(x = G_std, y = spec_perc_rich)) +
  geom_point(size = .5) + 
  stat_smooth(method = "loess", col = "olivedrab",
              fill = "olivedrab") +
  ylab("Percentage monophagous species") +
  xlab("Grazing intensity") +
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20, 25),
                     labels = paste0(c(0, 5, 10, 15, 20, 25), "%"),
                     limits = c(0, NA))

plot_grid(plot_grid(ggplot_labeller(expression(bold("a")), 0, label.size, NA), empty, 
                    ggplot_labeller(expression(bold("b")), 0, label.size, NA), empty, nrow = 1,
                    rel_widths = c(1, 5, 1, 5)),
          plot_grid(p1, p2, nrow = 1), nrow = 2, rel_heights = c(1, 10))

#  FORESTS / WINDOW TRAPS -------------------------------- #####################
################################################################################.

#  prepare data and interaction matrices #######################################
################################################################################.

# herbivores -------------------------------------------------------------------.

herb_list_wintr <- wintr_agg %>% 
  select(Herb_Name_std) %>% 
  arrange(Herb_Name_std) %>% 
  unique()

wintr_ct <- wintr_agg %>% 
  select(PlotID, Herb_Name_std, NumberAdults) %>% 
  spread(Herb_Name_std, NumberAdults, fill = 0) %>% 
  as.data.frame() %>% 
  column_to_rownames("PlotID")

# weight omnivore and mostly carnivore/fungivore species less:
feeding_weights <- herb_list_wintr %>% 
  left_join(herb_splist, by = "Herb_Name_std") %>% 
  mutate(feeding_weight = ifelse(Feeding_guild_short == "h", 1,
                                 ifelse(Feeding_guild_short == "o",
                                        0.5, 0.25))) %>% 
  select(Herb_Name_std, feeding_weight)
wintr_ct_fw <- as.data.frame(t(feeding_weights$feeding_weight * t(wintr_ct)))

# weight with metabolic rate:
lw_wintr <- herb_list_wintr$Herb_Name_std %>% 
  future_map(f_lw, regs_mass_length) %>% 
  do.call(rbind, .)

lw_wintr <- herb_list_wintr %>% 
  left_join(lw_wintr, by = "Herb_Name_std") %>% 
  left_join(select(herb_splist, Herb_Name_std, Body_Size), by = "Herb_Name_std") 

lw_wintr <- lw_wintr %>% 
  mutate(weight = 10 ^ a * (Body_Size ^ b),
         metabolic_rel = weight ^ 0.75)

wintr_ct_mr <- 
  wintr_agg %>% 
  filter(Herb_Name_std %in% herb_list_wintr$Herb_Name_std) %>%
  left_join(lw_wintr, by = "Herb_Name_std") %>%
  mutate(metabolic_entities = NumberAdults * metabolic_rel) %>% 
  select(PlotID, Herb_Name_std, metabolic_entities) %>% 
  spread(Herb_Name_std, metabolic_entities, fill = 0) %>% 
  as.data.frame() %>% 
  column_to_rownames("PlotID")

# combine feeding weighting and metabolic rate:
wintr_ct_fw_mr <- as.data.frame(t(feeding_weights$feeding_weight * t(wintr_ct_mr)))


# plants -----------------------------------------------------------------------.

plant_ct_f <- plant_agg_f %>% 
  filter(PlotID %in% rownames(wintr_ct)) %>% 
  select(PlotID, Host_Name_std, Cover) %>% 
  spread(Host_Name_std, Cover, fill = 0) %>% 
  as.data.frame() %>% 
  column_to_rownames("PlotID")

# biomass instead of cover:
plant_ct_f_bm <- plant_agg_f %>% 
  filter(PlotID %in% rownames(wintr_ct)) %>% 
  select(PlotID, Host_Name_std, biomass) %>% 
  spread(Host_Name_std, biomass, fill = 0) %>% 
  as.data.frame() %>% 
  column_to_rownames("PlotID")


# all plants present in very low numbers
plant_ct_f_nozero <- plant_ct_f
plant_ct_f_nozero[plant_ct_f_nozero == 0] <- 10^-25

# global interaction matrix ----------------------------------------------------.

# weighting all taxonomic level information the same
tax_weights_nw <- data.frame(group = ordered(host.rank.levels, 
                                             levels = host.rank.levels))
tax_weights_nw$weight <- 1
poly_weight_nw <- 1


int_mat_wintr <- mclapply(herb_list_wintr$Herb_Name_std, 
                          construct_int_matrix, plant_agg_f, 
                          tax_weights_nw, poly_weight_nw)
int_mat_wintr <- do.call(bind_rows, int_mat_wintr)
int_mat_wintr[is.na(int_mat_wintr)] <- 0
int_mat_wintr <- int_mat_wintr[ , sort(names(int_mat_wintr))]
int_mat_wintr <- int_mat_wintr %>% 
  mutate(ID = paste(Herb_Name_std, stage, sep = "_")) %>% 
  column_to_rownames("ID") %>% 
  select(-c(Herb_Name_std, stage))


# weighting  different taxonomic level information differently
tax_weights_pw <- data.frame(group = ordered(host.rank.levels, 
                                             levels = host.rank.levels))
tax_weights_pw <- tax_weights_pw %>% 
  mutate(weight = ifelse(group < "Family", 0.25,
                         ifelse(group < "Genus", 0.5, 1)))
poly_weight_pw <- 0.25

int_mat_wintr_pw <- mclapply(herb_list_wintr$Herb_Name_std, 
                             construct_int_matrix, plant_agg_f, 
                             tax_weights_pw, poly_weight_pw)
int_mat_wintr_pw <- do.call(bind_rows, int_mat_wintr_pw)
int_mat_wintr_pw[is.na(int_mat_wintr_pw)] <- 0
int_mat_wintr_pw <- int_mat_wintr_pw[ , sort(names(int_mat_wintr_pw))]
int_mat_wintr_pw <- int_mat_wintr_pw %>% 
  mutate(ID = paste(Herb_Name_std, stage, sep = "_")) %>% 
  column_to_rownames("ID") %>% 
  select(-c(Herb_Name_std, stage))

# no polyphagous interactions
tax_weights_no_poly <- data.frame(group = ordered(host.rank.levels, 
                                                  levels = host.rank.levels))
tax_weights_no_poly <- tax_weights_no_poly %>% 
  mutate(weight = ifelse(group < "Family", 0 ,1))
int_mat_wintr_no_poly <- mclapply(herb_list_wintr$Herb_Name_std, 
                                  construct_int_matrix, plant_agg_f, 
                                  tax_weights_no_poly, poly_weight = 0,
                                  poly_all = F, poly_sub = F)
int_mat_wintr_no_poly <- do.call(bind_rows, int_mat_wintr_no_poly)
# add null columns:
nullcols  <-  unique(plant_agg_f$Host_Name_std[!plant_agg_f$Host_Name_std %in% names(int_mat_wintr_no_poly)])
int_mat_wintr_no_poly[, nullcols] <- NA; rm(nullcols)
int_mat_wintr_no_poly[is.na(int_mat_wintr_no_poly)] <- 0
int_mat_wintr_no_poly <- int_mat_wintr_no_poly[ , sort(names(int_mat_wintr_no_poly))]
int_mat_wintr_no_poly <- int_mat_wintr_no_poly %>% 
  mutate(ID = paste(Herb_Name_std, stage, sep = "_")) %>% 
  column_to_rownames("ID") %>% 
  select(-c(Herb_Name_std, stage))

# calculate global interaction matrices ########################################
################################################################################.

int_global_wintr_raw <- rownames(wintr_ct) %>% 
  map(f.int.global, int_mat_wintr, wintr_ct, plant_ct_f, stage_weights) %>% 
  do.call(rbind, .)

int_global_wintr_fw <- rownames(wintr_ct_fw) %>% 
  map(f.int.global, int_mat_wintr, wintr_ct_fw, plant_ct_f, stage_weights) %>% 
  do.call(rbind, .)

int_global_wintr_pw <- rownames(wintr_ct) %>% 
  map(f.int.global, int_mat_wintr_pw, wintr_ct, plant_ct_f, stage_weights) %>% 
  do.call(rbind, .)

int_global_wintr_mr <- rownames(wintr_ct_mr) %>% 
  map(f.int.global, int_mat_wintr, wintr_ct_mr, plant_ct_f, stage_weights) %>% 
  do.call(rbind, .)

int_global_wintr_bm <- rownames(wintr_ct) %>% 
  map(f.int.global, int_mat_wintr, wintr_ct, plant_ct_f_bm, stage_weights) %>% 
  do.call(rbind, .)

int_global_wintr_fw_pw <- rownames(wintr_ct_fw) %>% 
  map(f.int.global, int_mat_wintr_pw, wintr_ct_fw, plant_ct_f, stage_weights) %>% 
  do.call(rbind, .)

int_global_wintr_fw_mr <- rownames(wintr_ct_fw_mr) %>% 
  map(f.int.global, int_mat_wintr, wintr_ct_fw_mr, plant_ct_f, stage_weights) %>% 
  do.call(rbind, .)

int_global_wintr_fw_bm <- rownames(wintr_ct_fw) %>% 
  map(f.int.global, int_mat_wintr, wintr_ct_fw, plant_ct_f_bm, stage_weights) %>% 
  do.call(rbind, .)

int_global_wintr_pw_mr <- rownames(wintr_ct_mr) %>% 
  map(f.int.global, int_mat_wintr_pw, wintr_ct_mr, plant_ct_f, stage_weights) %>% 
  do.call(rbind, .)

int_global_wintr_pw_bm <- rownames(wintr_ct) %>% 
  map(f.int.global, int_mat_wintr_pw, wintr_ct, plant_ct_f_bm, stage_weights) %>% 
  do.call(rbind, .)

int_global_wintr_mr_bm <- rownames(wintr_ct_mr) %>% 
  map(f.int.global, int_mat_wintr, wintr_ct_mr, plant_ct_f_bm, stage_weights) %>% 
  do.call(rbind, .)

int_global_wintr_fw_pw_mr <- rownames(wintr_ct_fw_mr) %>% 
  map(f.int.global, int_mat_wintr_pw, wintr_ct_fw_mr, plant_ct_f, stage_weights) %>% 
  do.call(rbind, .)

int_global_wintr_fw_pw_bm <- rownames(wintr_ct_fw) %>% 
  map(f.int.global, int_mat_wintr_pw, wintr_ct_fw, plant_ct_f_bm, stage_weights) %>% 
  do.call(rbind, .)

int_global_wintr_fw_mr_bm <- rownames(wintr_ct_fw_mr) %>% 
  map(f.int.global, int_mat_wintr, wintr_ct_fw_mr, plant_ct_f_bm, stage_weights) %>% 
  do.call(rbind, .)

int_global_wintr_pw_mr_bm <- rownames(wintr_ct_mr) %>% 
  map(f.int.global, int_mat_wintr_pw, wintr_ct_mr, plant_ct_f_bm, stage_weights) %>% 
  do.call(rbind, .)

int_global_wintr_fw_pw_mr_bm <- rownames(wintr_ct_fw_mr) %>% 
  map(f.int.global, int_mat_wintr_pw, wintr_ct_fw_mr, plant_ct_f_bm, stage_weights) %>% 
  do.call(rbind, .)

int_global_wintr_all_plants <- rownames(wintr_ct) %>% 
  map(f.int.global, int_mat_wintr, wintr_ct, plant_ct_f_nozero, stage_weights) %>% 
  do.call(rbind, .)

int_global_wintr_no_poly <- rownames(wintr_ct) %>% 
  map(f.int.global, int_mat_wintr_no_poly, wintr_ct, plant_ct_f, stage_weights) %>% 
  do.call(rbind, .)

int_global_wintr_pres_abs <- int_global_wintr_raw %>% 
  mutate(interaction = 1)


# calculate  metrics ###########################################################
################################################################################.

nw.metrics.wintr <- list()
mod <- T

# "raw" metrics ----------------------------------------------------------------.
tmp <- mclapply(rownames(wintr_ct), f.nw.metrics,
                int_global_wintr_raw,
                mod = mod)
nw.metrics.wintr[["raw"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# correction: feeding type -----------------------------------------------------.
tmp <- mclapply(rownames(wintr_ct_fw), f.nw.metrics,
                int_global_wintr_fw,
                mod = mod)
nw.metrics.wintr[["fw"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# correction: taxonomic information --------------------------------------------.
tmp<- mclapply(rownames(wintr_ct), f.nw.metrics,
               int_global_wintr_pw,
               mod = mod)
nw.metrics.wintr[["pw"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# correction: metabolic rate ---------------------------------------------------.
tmp<- mclapply(rownames(wintr_ct_mr), f.nw.metrics,
               int_global_wintr_mr,
               mod = mod)
nw.metrics.wintr[["mr"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# correction: plant biomass ----------------------------------------------------.
tmp<- mclapply(rownames(wintr_ct), f.nw.metrics,
               int_global_wintr_bm,
               mod = mod)
nw.metrics.wintr[["bm"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# correction: feeding type, taxonomic information ------------------------------.
tmp <- mclapply(rownames(wintr_ct_fw), f.nw.metrics,
                int_global_wintr_fw_pw,
                mod = mod)
nw.metrics.wintr[["fw.pw"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# correction: feeding type, metabolic rate -------------------------------------.
tmp <- mclapply(rownames(wintr_ct_fw_mr), f.nw.metrics,
                int_global_wintr_fw_mr,
                mod = mod)
nw.metrics.wintr[["fw.mr"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# correction: feeding type, plant biomass --------------------------------------.
tmp <- mclapply(rownames(wintr_ct_fw), f.nw.metrics,
                int_global_wintr_fw_bm,
                mod = mod)
nw.metrics.wintr[["fw.bm"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# correction: taxonomic information, metabolic rate ----------------------------.
tmp <- mclapply(rownames(wintr_ct_mr), f.nw.metrics,
                int_global_wintr_pw_mr,
                mod = mod)
nw.metrics.wintr[["pw.mr"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# correction: taxonomic information, plant biomass -----------------------------.
tmp <- mclapply(rownames(wintr_ct), f.nw.metrics,
                int_global_wintr_pw_bm,
                mod = mod)
nw.metrics.wintr[["pw.bm"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# correction: metabolic rate, plant biomass ------------------------------------.
tmp <- mclapply(rownames(wintr_ct_mr), f.nw.metrics,
                int_global_wintr_mr_bm,
                mod = mod)
nw.metrics.wintr[["mr.bm"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# correction: feeding type, taxonomic information, metabolic rate --------------.
tmp <- mclapply(rownames(wintr_ct_fw_mr), f.nw.metrics,
                int_global_wintr_fw_pw_mr,
                mod = mod)
nw.metrics.wintr[["fw.pw.mr"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# correction: feeding type, taxonomic information, plant biomass ---------------.
tmp <- mclapply(rownames(wintr_ct_fw), f.nw.metrics,
                int_global_wintr_fw_pw_bm,
                mod = mod)
nw.metrics.wintr[["fw.pw.bm"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# correction: feeding type, metabolic rate, plant biomass ----------------------.
tmp <- mclapply(rownames(wintr_ct_fw_mr), f.nw.metrics,
                int_global_wintr_fw_mr_bm,
                mod = mod)
nw.metrics.wintr[["fw.mr.bm"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# correction: taxonomic information, metabolic rate, plant biomass -------------.
tmp <- mclapply(rownames(wintr_ct_mr), f.nw.metrics,
                int_global_wintr_pw_mr_bm,
                mod = mod)
nw.metrics.wintr[["pw.mr.bm"]] <- do.call(rbind.data.frame, tmp); rm(tmp)

# correction: feeding type, taxonomic information, metabolic rate, plant biomass
tmp <- mclapply(rownames(wintr_ct_fw_mr), f.nw.metrics,
                int_global_wintr_fw_pw_mr_bm,
                mod = mod)
nw.metrics.wintr[["fw.pw.mr.bm"]] <- do.call(rbind.data.frame, tmp); rm(tmp)


# no polyphagous species -------------------------------------------------------.
tmp <- mclapply(rownames(wintr_ct), f.nw.metrics,
                int_global_wintr_no_poly,
                mod = mod)
nw.metrics.wintr[["no.poly"]] <- do.call(rbind.data.frame, tmp); rm(tmp)


# all plants present -----------------------------------------------------------.
tmp <- mclapply(rownames(wintr_ct), f.nw.metrics,
                int_global_wintr_all_plants,
                mod = mod)
nw.metrics.wintr[["all.plants"]] <- do.call(rbind.data.frame, tmp); rm(tmp)


# PRESENCE - ABSENCE -----------------------------------------------------------.
tmp <- mclapply(rownames(wintr_ct), f.nw.metrics,
                int_global_wintr_pres_abs,
                mod = mod)
nw.metrics.wintr[["pres.abs"]] <- do.call(rbind.data.frame, tmp); rm(tmp)


# hierarchical models ##########################################################
################################################################################.

# define function
hier_mod_wintr <- function(data){

  data <- data %>%
    left_join(formi, by = "PlotID")
  
  out <- list()
  for (res in responses[responses %in% names(data)]){
    
    data$response <- (data[, res] - mean(data[, res]))/sd(data[, res])
    
    data.comb <- list(J = length(unique(data$Exploratory)),
                      N = nrow(data),
                      Exploratory = as.integer(as.factor(data$Exploratory)),
                      ForMI_z = as.vector(data$ForMI_z),
                      y = as.vector(data$response))
    
    
    fit <- sampling(object = stan_formi_mod, data = data.comb,
                    chains = 4, warmup = 500, iter = 5500, control = list(adapt_delta = 0.99))
    out[[res]] <- list(formi_comb = extract(fit)[c("a", "b", "mu_a", "sigma_a", "sigma_b", "sigma_y")])
    
    
    
    data.split <- list(J = length(unique(data$Exploratory)),
                       N = nrow(data),
                       Exploratory = as.integer(as.factor(data$Exploratory)),
                       Inonat_z = as.vector(data$Inonat_z),
                       Idwcut_z = as.vector(data$Idwcut_z),
                       Iharv_z = as.vector(data$Iharv_z),
                       y = as.vector(data$response))
    
    
    fit <- sampling(object = stan_formisplit_mod, data = data.split,
                    chains = 4, warmup = 500, iter = 5500, control = list(adapt_delta = 0.99))
    
    out[[res]] <- c(out[[res]], list(formi_split = extract(fit)[c("a", "b", "mu_a", "sigma_a", "sigma_b_non", "sigma_b_cut", "sigma_b_harv", "sigma_y")]))
  }
  
  out
}

# apply function
stats_b_wintr <- lapply(nw.metrics.wintr, hier_mod_wintr)
names(stats_b_wintr) <- names(nw.metrics.wintr)

# null models: on external cluster  ############################################
################################################################################.

# data export ------------------------------------------------------------------.

plots <- rownames(wintr_ct)

save(
     int_global_wintr_raw, 
     int_global_wintr_fw, int_global_wintr_pw,
     int_global_wintr_mr, int_global_wintr_bm, int_global_wintr_fw_pw,
     int_global_wintr_fw_mr, int_global_wintr_fw_bm, int_global_wintr_pw_mr,
     int_global_wintr_pw_bm, int_global_wintr_mr_bm, int_global_wintr_fw_pw_mr,
     int_global_wintr_fw_pw_bm, int_global_wintr_fw_mr_bm,
     int_global_wintr_pw_mr_bm, int_global_wintr_fw_pw_mr_bm,
     int_global_wintr_no_poly, int_global_wintr_all_plants,
     int_global_wintr_pres_abs,
     plots, 
     formi, nw.metrics.wintr,
     file = "wintr_data.Rdata")

# data import ------------------------------------------------------------------.
# after analysis

fit_wintr_null <- list()

loop.files <- list.files("dir_nullmodels",
                         full.names = T)
loop.files <- loop.files[grepl("wintr", loop.files)]
loop.names <- gsub("wintr_", "", list.files("dir_nullmodels",
                                            full.names = F))
loop.names <- loop.names[!grepl("swnet", loop.names) & !grepl("system", loop.names) &
                           !grepl("150", loop.names)]


for (j in 1:length(loop.files)){
  tmp.files <- list.files(loop.files[j], full.names = T)
  if(length(tmp.files) > 0){
    tmp.names <- gsub(".txt", "", list.files(loop.files[j], full.names = F))
    
    # set values for the colClasses argument in read.csv. Improves speed (maybe? not too much)
    col.test <- read.table(tmp.files[1])
    cols <- ifelse(names(col.test) %in% c("response", "null.model", 
                                          "formi_comb.a.1", "formi_comb.a.2", "formi_comb.a.3",
                                          "formi_comb.b", 
                                          "formi_split.a.1", "formi_split.a.2", "formi_split.a.3",
                                          "formi_split.b.1","formi_split.b.2", 
                                          "formi_split.b.3"),
                   ifelse(names(col.test) %in% c("response", "null.model"), "factor", "numeric"), 
                   "NULL")
    names(cols) <- names(col.test)
    
    for (i in 1:length(tmp.files)){
      fit_wintr_null[[gsub("_", ".", loop.names[j])]][[tmp.names[i]]] <- read.table(tmp.files[i], colClasses = cols)
      
    }
  }

}

# plotting main figure #########################################################
################################################################################.

y_min_thresh <- 0.01

n.nullshades <- 100

for (datatype.def in c("raw", "no.poly", "pres.abs", "all.plants")){
  
  posterior_plots_wintr_def <- list()
  posterior_plots_wintr_int <- list()
  y_axes_wintr <- list()
  range_x_wintr <- list()
  
  for (res in responses){
    if (res %in% levels(fit_wintr_null[[datatype.def]]$data_00001$response)){
      
      # gather all data ----------------------------------------------------------.
      
      comb_nullmodel <- list()
      density_obs <- list()
      
      
      density_obs[["ForMI"]] <- data.frame(density(stats_b_wintr[[datatype.def]][[res]][["formi_comb"]]$b)) %>% 
        filter(y > y_min_thresh)
      
      for (i in 1:length(fit_wintr_null[[datatype.def]])){
        
        comb_nullmodel[["ForMI"]] <- c(comb_nullmodel[["ForMI"]],
                                       filter(fit_wintr_null[[datatype.def]][[i]],
                                              response == res)$formi_comb.b)
      }
      
      
      InIdIh_names <- c("Inonat", "Idwcut", "Iharv")
      for (InIdIh in 1:3){
        
        density_obs[[InIdIh_names[InIdIh]]] <- data.frame(density(stats_b_wintr[[datatype.def]][[res]][["formi_split"]]$b[, InIdIh])) %>% 
          filter(y > y_min_thresh)
        
        
        
        for (i in 1:length(fit_wintr_null[[datatype.def]])){
          
          comb_nullmodel[[InIdIh_names[InIdIh]]] <- c(comb_nullmodel[[InIdIh_names[InIdIh]]],
                                                      filter(fit_wintr_null[[datatype.def]][[i]],
                                                             response == res)[ , paste0("formi_split.b.", InIdIh)])
          
        }
      }
      
      max_y <- max(c(sapply(comb_nullmodel, function(x) max(density(x)$y)),
                     sapply(density_obs, function(x) max(x$y))))
      
      
      # ForMI comb -----------------------------------------------------------------.
      p <- ggplot(data = density_obs[["ForMI"]], 
                  aes(x, y)) +
        xlab(res)
      
      
      density_null1 <- data.frame(density(comb_nullmodel[["ForMI"]])) %>% 
        filter(y > y_min_thresh) 
      
      
      for (i in 1:n.nullshades){
        d_target <- data.frame(density(filter(fit_wintr_null[[datatype.def]][[i]],
                                              response == res)$formi_comb.b)) %>% 
          filter(y > y_min_thresh) %>%
          mutate(y = y / (.75 * max_y),
                 y = y + 1)
        
        p <-
          p +
          geom_ribbon(data = d_target, 
                      aes(ymin = 1, ymax = y),
                      fill = "#7570b3", alpha = .1) 
        
      }
      
      hdis_core <- rbind(nullmodel = hdi(comb_nullmodel[["ForMI"]], ci = 1 - sqrt(0.05)), # argument: chance to be outside of both core areas (no overlap) is alpha * alpha
                         observed = hdi(stats_b_wintr[[datatype.def]][[res]][["formi_comb"]]$b, ci = 1 - sqrt(0.05))) %>% 
        rownames_to_column("data") %>% 
        mutate(data = factor(data, levels = .$data))
      
      hdis_outer <- rbind(nullmodel = hdi(comb_nullmodel[["ForMI"]], ci = .95), 
                          observed = hdi(stats_b_wintr[[datatype.def]][[res]][["formi_comb"]]$b, ci = .95)) %>% 
        rownames_to_column("data") %>% 
        mutate(data = factor(data, levels = .$data))
      
      est <- rbind(nullmodel = map_estimate(comb_nullmodel[["ForMI"]]), 
                   observed = map_estimate(stats_b_wintr[[datatype.def]][[res]][["formi_comb"]]$b)) %>% 
        rownames_to_column("data") %>% 
        mutate(data = factor(data, levels = .$data))
      
      p <-
        p +
        geom_vline(xintercept = 0, lty = "dashed") +
        geom_line(data = density_null1 %>% 
                    mutate(y = y / (.75 * max_y),
                           y = y + 1),
                  col = "white", size = 1.25) +
        geom_line(data = density_null1 %>% 
                    mutate(y = y / (.75 * max_y),
                           y = y + 1),
                  col = "#7570b3", size = .5) +
        geom_line(data = density_obs[["ForMI"]] %>% 
                    mutate(y = y / (.75 * max_y),
                           y = y + 1),
                  col = "white", size = 1.25) +
        geom_line(data = density_obs[["ForMI"]] %>% 
                    mutate(y = y / (.75 * max_y),
                           y = y + 1),
                  col = "#d95f02", size = .5) +
        geom_segment(data = hdis_outer, aes(x = CI_low, xend = CI_high, 
                                            y = as.numeric(data)/3 + .1, yend = as.numeric(data)/3 + .1,
                                            col = data), size = .3) +
        geom_segment(data = hdis_core, aes(x = CI_low, xend = CI_high, 
                                           y = as.numeric(data)/3 + .1, yend = as.numeric(data)/3 + .1,
                                           col = data), size = 1.25) +
        geom_point(data = est, aes(x = MAP, y = as.numeric(data)/3 + .1), 
                   col = "white", size = 2) +
        geom_point(data = est, aes(x = MAP, y = as.numeric(data)/3 + .1, 
                                   col = data), size = 1.25) +
        scale_color_manual(values = c(nullmodel = "#7570b3", 
                                      observed = "#d95f02"),
                           guide = F) + 
        theme(axis.title = element_blank(),
              axis.line.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank()) +
        scale_y_continuous(limits = c(0.2 , 3))

      posterior_plots_wintr_def[[res]][["ForMI"]] <- p
      
      range_x_wintr[[res]]$min <- min(density_null1$x, 
                                      density_obs[["ForMI"]]$x)
      range_x_wintr[[res]]$max <- max(density_null1$x, 
                                      density_obs[["ForMI"]]$x)
      
      
      # intercept: 
      
      comb_nullmodel_int <- c()
      
      for (i in 1:length(fit_wintr_null[[datatype.def]])){
        
        comb_nullmodel_int <- c(comb_nullmodel_int,
                                rowMeans(filter(fit_wintr_null[[datatype.def]][[i]],
                                                response == res)[c("formi_comb.a.1", 
                                                                   "formi_comb.a.2", 
                                                                   "formi_comb.a.3")]))
      }
      
      
      
      sd_target <- sd(nw.metrics.wintr[[datatype.def]][[res]])
      m_target <- mean(nw.metrics.wintr[[datatype.def]][[res]])
      
      hdis_core_int <- rbind(nullmodel = hdi(comb_nullmodel_int * sd_target + m_target, 
                                             ci = 1 - sqrt(0.05)), # argument: chance to be outside of both core areas (no overlap) is alpha * alpha
                             observed = hdi(rowMeans(stats_b_wintr[[datatype.def]][[res]][["formi_comb"]]$a * sd_target + m_target), ci = 1 - sqrt(0.05))) %>% 
        rownames_to_column("data") %>% 
        mutate(data = factor(data, levels = .$data))
      
      hdis_outer_int <- rbind(nullmodel = hdi(comb_nullmodel_int * sd_target + m_target, ci = .95), 
                              observed = hdi(rowMeans(stats_b_wintr[[datatype.def]][[res]][["formi_comb"]]$a * sd_target + m_target), ci = .95)) %>% 
        rownames_to_column("data") %>% 
        mutate(data = factor(data, levels = .$data))
      
      est_int <- rbind(nullmodel = map_estimate(comb_nullmodel_int * sd_target + m_target), 
                       observed = map_estimate(rowMeans(stats_b_wintr[[datatype.def]][[res]][["formi_comb"]]$a * sd_target + m_target))) %>% 
        rownames_to_column("data") %>% 
        mutate(data = factor(data, levels = .$data))
      
      
      posterior_plots_wintr_int[[res]] <-
        ggplot() +
        geom_segment(data = hdis_outer_int, aes(y = CI_low, yend = CI_high, 
                                                x = data, xend = data,
                                                col = data), size = .3) +
        geom_segment(data = hdis_core_int, aes(y = CI_low, yend = CI_high, 
                                               x = data, xend = data,
                                               col = data), size = 1.25) +
        geom_point(data = est_int, aes(y = MAP, x = data), 
                   col = "white", size = 2) +
        geom_point(data = est_int, aes(y = MAP, x = data, 
                                       col = data), size = 1.25) +
        scale_color_manual(values = c(nullmodel = "#7570b3", 
                                      observed = "#d95f02"),
                           guide = F) + 
        theme(axis.title = element_blank(),
              axis.line.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank()) 
      
      
      # ForMI split ----------------------------------------------------------------.
      
      for (InIdIh in 1:3){
        
        p <- ggplot(data = density_obs[[InIdIh_names[InIdIh]]], 
                    aes(x, y)) +
          xlab(res)
        
        
        
        density_null1 <- data.frame(density(comb_nullmodel[[InIdIh_names[InIdIh]]])) %>% 
          filter(y > y_min_thresh)
        
        
        for (i in 1:n.nullshades){
          
          p <-
            p +
            geom_ribbon(data = data.frame(density(filter(fit_wintr_null[[datatype.def]][[i]],
                                                         response == res)[ , paste0("formi_split.b.", InIdIh)])) %>% 
                          filter(y > y_min_thresh) %>% 
                          mutate(y = y / (.75 * max_y),
                                 y = y + 1), 
                        aes(ymin = 1, ymax = y),
                        fill = "#7570b3", alpha = .1) 
          
        }
        
        
        hdis_core <- rbind(nullmodel = hdi(comb_nullmodel[[InIdIh_names[InIdIh]]], ci = 1 - sqrt(0.05)), 
                           observed = hdi(stats_b_wintr[[datatype.def]][[res]][["formi_split"]]$b[, InIdIh], ci = 1 - sqrt(0.05))) %>% 
          rownames_to_column("data") %>% 
          mutate(data = factor(data, levels = .$data))
        
        hdis_outer <- rbind(nullmodel = hdi(comb_nullmodel[[InIdIh_names[InIdIh]]], ci = .95), 
                            observed = hdi(stats_b_wintr[[datatype.def]][[res]][["formi_split"]]$b[, InIdIh], ci = .95)) %>% 
          rownames_to_column("data") %>% 
          mutate(data = factor(data, levels = .$data))
        
        est <- rbind(nullmodel = map_estimate(comb_nullmodel[[InIdIh_names[InIdIh]]]), 
                     observed = map_estimate(stats_b_wintr[[datatype.def]][[res]][["formi_split"]]$b[, InIdIh])) %>% 
          rownames_to_column("data") %>% 
          mutate(data = factor(data, levels = .$data))
        
        p <-
          p +
          geom_vline(xintercept = 0, lty = "dashed") +
          geom_line(data = density_null1 %>% 
                      mutate(y = y / (.75 * max_y),
                             y = y + 1),
                    col = "white", size = 1.25) +
          geom_line(data = density_null1 %>% 
                      mutate(y = y / (.75 * max_y),
                             y = y + 1),
                    col = "#7570b3", size = .5) +
          geom_line(data = density_obs[[InIdIh_names[InIdIh]]] %>% 
                      mutate(y = y / (.75 * max_y),
                             y = y + 1),
                    col = "white", size = 1.25) +
          geom_line(data = density_obs[[InIdIh_names[InIdIh]]] %>% 
                      mutate(y = y / (.75 * max_y),
                             y = y + 1),
                    col = "#d95f02", size = .5) +
          geom_segment(data = hdis_outer, aes(x = CI_low, xend = CI_high, 
                                              y = as.numeric(data)/3 + .1, yend = as.numeric(data)/3 + .1,
                                              col = data), size = .3) +
          geom_segment(data = hdis_core, aes(x = CI_low, xend = CI_high, 
                                             y = as.numeric(data)/3 + .1, yend = as.numeric(data)/3 + .1,
                                             col = data), size = 1.25) +
          geom_point(data = est, aes(x = MAP, y = as.numeric(data)/3 + .1), 
                     col = "white", size = 2) +
          geom_point(data = est, aes(x = MAP, y = as.numeric(data)/3 + .1, 
                                     col = data), size = 1.25) +
          scale_color_manual(values = c(nullmodel = "#7570b3", 
                                        observed = "#d95f02"),
                             guide = F) + 
          theme(axis.title = element_blank(),
                axis.line.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank()) +
          scale_y_continuous(limits = c(0.2 , 3))

        posterior_plots_wintr_def[[res]][[InIdIh_names[InIdIh]]] <- p
        
        range_x_wintr[[res]]$min <- c(range_x_wintr[[res]]$min, min(density_null1$x, 
                                                                    density_obs[[InIdIh_names[InIdIh]]]$x))
        range_x_wintr[[res]]$max <- c(range_x_wintr[[res]]$max, max(density_null1$x, 
                                                                    density_obs[[InIdIh_names[InIdIh]]]$x))
        
      }
      
      
      # y axes plots -------------------------------------------------------------.
      
      y_axes_wintr[[res]] <-
        ggplot() + 
        scale_y_continuous(limits = c(0.2 , 3),
                           breaks = c(0, floor(max_y)) / (0.75 * max_y) + 1,
                           labels = c(0, floor(max_y)),
                           position = "right") +
        scale_x_continuous(limits = c(0, 0)) +
        annotate(geom = "segment", x = Inf, xend = Inf, y = 1, yend = 1 / (0.75) + 1) +
        theme(axis.line.y = element_blank(),
              axis.line.x = element_blank(),
              axis.text.x = element_blank(),
              axis.title = element_blank(),
              axis.ticks.x = element_blank(),
              plot.margin = margin(r = 0))
    }
  }
  
  axis.text.size <- 8
  axis.title.size <- 10
  label.size <- 12
  
  posterior_plots_wintr_int <- lapply(posterior_plots_wintr_int, function(x) x + 
                                        theme(axis.text = element_text(size = axis.text.size))+
                                        scale_y_continuous(position = "right")) 
  
  y_axes_wintr <- lapply(y_axes_wintr, function(x) x + 
                           theme(axis.text = element_text(size = axis.text.size)))
  
  
  legend <-
    ggplot(data.frame(var = factor(c("observed", "null model"), 
                                   levels = c("observed", "null model")),
                      x = NA, y = NA)) +
    geom_area(aes(fill = var, x = x, y = y)) +
    scale_fill_manual(values = c(observed = "#d95f02",
                                 "null model" = "#7570b3"))  +
    theme_nothing() +
    theme(legend.position = c(.5, .5),
          legend.justification = c(.5, .5),
          legend.title = element_blank(),
          legend.background = element_rect(fill = col.strip,
                                           colour = NA),
          legend.text = element_text(size = axis.title.size),
          legend.key.height = unit(2, "mm"))
  
  empty <- ggplot() + geom_blank()
  
  
  range_x <- c(min(c(range_x_wintr$nested$min, range_x_wintr$qmod$min, range_x_wintr$rob.A$min)),
               max(c(range_x_wintr$nested$max, range_x_wintr$qmod$max, range_x_wintr$rob.A$max)))
  
  
  rel_heights <- c(.22, .22, .9, .9, .9, .32, .18)
  rel_heights <- rel_heights / sum(rel_heights)
  
  rel_widths <- c(.26, 1, 1, 1, 1) # metric label, 4 figures (land-use components)
  rel_widths <- rel_widths / sum(rel_widths)
  rel_widths_sub <- c(sum(rel_widths[2:5]), rel_widths[2] / 3, rel_widths[2] / 4.1) # core figure, probability density axis, probability density title
  
  rel_widths_overall <- c(4.7, .2, 1)
  
  
  col.strip <- "grey80"
  
  
  p <- plot_grid(
    # Slope strip:
    plot_grid(plot_grid(ggplot_labeller(expression(bold("a")), 0, label.size, NA),
                        ggplot_labeller("Slope", 0, label.size, col.strip),
                        nrow = 1, rel_widths = c(rel_widths[1], sum(rel_widths_sub))),
              
              plot_grid(
                # core figure (slope):
                plot_grid(plot_grid(empty,
                                    ggplot_labeller("ForMI", 0, label.size, col.strip) +
                                      theme(plot.background = element_rect(fill = "transparent", colour = NA, size = NA)),
                                    ggplot_labeller("Inonat", 0, label.size, col.strip) +
                                      theme(plot.background = element_rect(fill = "transparent", colour = NA, size = NA)),
                                    ggplot_labeller("Iharv", 0, label.size, col.strip) +
                                      theme(plot.background = element_rect(fill = "transparent", colour = NA, size = NA)),
                                    ggplot_labeller("Idwcut", 0, label.size, col.strip) +
                                      theme(plot.background = element_rect(fill = "transparent", colour = NA, size = NA)),
                                    nrow = 1, rel_widths = rel_widths),
                          plot_grid(plotlist = c(list(ggplot_labeller("Modularity", 90, label.size, col.strip) +
                                                        theme(plot.background = element_rect(fill = "transparent", colour = NA, size = NA))),
                                                 lapply(posterior_plots_wintr_def$qmod,
                                                        function(x) x +
                                                          scale_x_continuous(limits = range_x * 1.05) +
                                                          theme_nothing())[c(1,2,4,3)]),
                                    nrow = 1, align = "h", rel_widths = rel_widths),
                          plot_grid(plotlist = c(list(ggplot_labeller("Nestedness", 90, label.size, col.strip) +
                                                        theme(plot.background = element_rect(fill = "transparent", colour = NA, size = NA))),
                                                 lapply(posterior_plots_wintr_def$nested,
                                                        function(x) x +
                                                          scale_x_continuous(limits = range_x * 1.05) +
                                                          theme_nothing())[c(1,2,4,3)]),
                                    nrow = 1, align = "h", rel_widths = rel_widths),
                          plot_grid(plotlist = c(list(ggplot_labeller("Robustness", 90, label.size, col.strip) +
                                                        theme(plot.background = element_rect(fill = "transparent", colour = NA, size = NA))),
                                                 lapply(posterior_plots_wintr_def$rob.A,
                                                        function(x) x +
                                                          scale_x_continuous(limits = range_x * 1.05) +
                                                          theme_nothing())[c(1,2,4,3)]),
                                    nrow = 1, align = "h", rel_widths = rel_widths),
                          plot_grid(plotlist =  c(list(empty),
                                                  rep(list(posterior_plots_wintr_def$rob.A$ForMI +
                                                             scale_x_continuous(limits = range_x * 1.05) +
                                                             scale_y_continuous(limits = c(0, 0)) +
                                                             theme(axis.text.x = element_text(size = axis.text.size,
                                                                                              angle = 90, hjust = 1, vjust = .5),
                                                                   plot.background = element_rect(fill = "transparent", colour = NA, size = NA),
                                                                   panel.background = element_rect(fill = "transparent", colour = NA, size = NA))), 4)),
                                    nrow = 1, rel_widths = rel_widths),
                          plot_grid(plotlist =  list(empty,
                                                     ggplot_labeller("Standardised model coefficient", 0, axis.title.size, NA)),
                                    nrow = 1, rel_widths = c(rel_widths[1], sum(rel_widths[2:5]))),
                          ncol = 1, rel_heights = rel_heights[-1]),
                
                # y-axis density:
                plot_grid(empty, 
                          y_axes_wintr$qmod,
                          y_axes_wintr$nested,
                          y_axes_wintr$rob.A,
                          empty, empty,
                          ncol = 1, rel_heights = rel_heights[-1], align = "v"),
                
                # y-axis density label:
                plot_grid(empty,
                          ggplot_labeller("Probability density", -90, axis.title.size, NA) +
                            theme(plot.margin = margin(r = 4)),
                          empty, empty,
                          ncol = 1,
                          rel_heights = c(rel_heights[2], sum(rel_heights[3:5]), 
                                          rel_heights[6:7])),
                
                nrow = 1, rel_widths = rel_widths_sub),
              ncol = 1, rel_heights = c(rel_heights[1], sum(rel_heights[-1]))
    ),
    
    #-------------------------------------------------------------------------------.
    ggplot() + 
      annotate(geom = "segment", x = 1, xend = 1, y = -.8, yend = Inf, lty = 3) +
      scale_y_continuous(limits = c(-1, 1)) +
      theme_nothing(),
    #-------------------------------------------------------------------------------.
    
    plot_grid(
      plot_grid(ggplot_labeller(expression(bold("b")), 0, label.size, NA),
                ggplot_labeller("Intercept", 0, label.size, col.strip),
                nrow = 1, rel_widths = c(1, 5 * rel_widths_overall[3])),
      plot_grid(plot_grid(plot_grid(empty,
                                    posterior_plots_wintr_int$qmod,
                                    posterior_plots_wintr_int$nested,
                                    posterior_plots_wintr_int$rob.A,
                                    ncol = 1, align = "v", rel_heights = rel_heights[2:5]),
                          plot_grid(empty,
                                    ggplot_labeller("Rescaled model coefficient", -90, axis.title.size, NA),
                                    rel_heights = c(rel_heights[2], sum(rel_heights[3:5])),
                                    ncol = 1),
                          nrow = 1, rel_widths = c(5, 1)),
                legend, 
                ncol = 1,
                rel_heights = c(sum(rel_heights[3:5]), sum(rel_heights[6:7]))),
      ncol = 1, rel_heights = c(rel_heights[1], sum(rel_heights[-1]))
    ),
    rel_widths = rel_widths_overall, nrow = 1)
  
}


# plotting (different corrections) #############################################
################################################################################.

y_min_thresh <- 0.01

posterior_plots_wintr_corr <- list()
range_x_wintr_corr <- list()

corrections <- c("fw", "pw", "mr", "bm",
                 "fw.pw", "fw.mr", "fw.bm", "pw.mr", "pw.bm", "mr.bm",
                 "fw.pw.mr", "fw.pw.bm", "fw.mr.bm", "pw.mr.bm",
                 "fw.pw.mr.bm", "raw")

for (datatype_i in corrections[corrections %in% names(stats_b_wintr)]){
  for (res in responses){
    if (res %in% levels(fit_wintr_null[[datatype_i]]$data_00001$response)){
      range_x_wintr_corr[[datatype_i]][[res]] <- c(0, 0)
      
      # LUI comb -----------------------------------------------------------------.
      density_obs <- data.frame(density(stats_b_wintr[[datatype_i]][[res]][["formi_comb"]]$b)) %>% 
        filter(y > y_min_thresh)
      p <- ggplot(data = density_obs, 
                  aes(x, y)) +
        xlab(res)
      
      comb_nullmodel <- c()
      
      for (i in 1:length(fit_wintr_null[[datatype_i]])){
        
        comb_nullmodel <- c(comb_nullmodel,
                            filter(fit_wintr_null[[datatype_i]][[i]],
                                   response == res)$formi_comb.b)
      }
      
      
      
      max_y <- max(density(comb_nullmodel)$y, 
                   p$data$y)
      
      density_null <- data.frame(density(comb_nullmodel)) %>% 
        filter(y > y_min_thresh) 
      
      
      
      hdis_core <- rbind(nullmodel = hdi(comb_nullmodel, ci = 1 - sqrt(0.05)), # argument: chance to be outside of both core areas (no overlap) is alpha * alpha
                         observed = hdi(stats_b_wintr[[datatype_i]][[res]][["formi_comb"]]$b, ci = 1 - sqrt(0.05))) %>% 
        rownames_to_column("data") %>% 
        mutate(data = factor(data, levels = .$data))
      
      hdis_outer <- rbind(nullmodel = hdi(comb_nullmodel, ci = .95), 
                          observed = hdi(stats_b_wintr[[datatype_i]][[res]][["formi_comb"]]$b, ci = .95)) %>% 
        rownames_to_column("data") %>% 
        mutate(data = factor(data, levels = .$data))
      
      est <- rbind(nullmodel = map_estimate(comb_nullmodel), 
                   observed = map_estimate(stats_b_wintr[[datatype_i]][[res]][["formi_comb"]]$b)) %>% 
        rownames_to_column("data") %>% 
        mutate(data = factor(data, levels = .$data))
      
      p <-
        p +
        geom_vline(xintercept = 0, lty = "dashed") +
        geom_line(data = density_null %>% 
                    mutate(y = y / (.75 * max_y),
                           y = y + 1),
                  col = "#7570b3", size = .75) +
        geom_line(data = density_obs %>% 
                    mutate(y = y / (.75 * max_y),
                           y = y + 1),
                  col = "#d95f02", size = .75) +
        geom_segment(data = hdis_outer, aes(x = CI_low, xend = CI_high, 
                                            y = as.numeric(data)/3 + .1, yend = as.numeric(data)/3 + .1,
                                            col = data), size = .3) +
        geom_segment(data = hdis_core, aes(x = CI_low, xend = CI_high, 
                                           y = as.numeric(data)/3 + .1, yend = as.numeric(data)/3 + .1,
                                           col = data), size = 1.25) +
        geom_point(data = est, aes(x = MAP, y = as.numeric(data)/3 + .1, 
                                   col = data), size = 1.75) +
        scale_color_manual(values = c(nullmodel = "#7570b3", 
                                      observed = "#d95f02"),
                           guide = F) + 
        theme(axis.title = element_blank(),
              axis.line.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank()) +
        scale_y_continuous(limits = c(0.2 , 3))

      posterior_plots_wintr_corr[[res]][["ForMI"]][[datatype_i]] <- p
      
      range_x_wintr_corr[[datatype_i]][[res]][1] <- min(density_null$x,
                                                        density_obs$x, 
                                                        range_x_wintr_corr[[datatype_i]][[res]][1])
      range_x_wintr_corr[[datatype_i]][[res]][2] <- max(density_null$x,
                                                        density_obs$x,
                                                        range_x_wintr_corr[[datatype_i]][[res]][2])
      
      
      # ForMI split ----------------------------------------------------------------.
      InIdIh_names <- c("Inonat", "Idwcut", "Iharv")
      for (InIdIh in 1:3){
        
        density_obs <- data.frame(density(stats_b_wintr[[datatype_i]][[res]][["formi_split"]]$b[, InIdIh])) %>% 
          filter(y > y_min_thresh)
        p <- ggplot(data = density_obs, 
                    aes(x, y)) +
          xlab(res)
        
        
        comb_nullmodel <- c()
        
        
        for (i in 1:length(fit_wintr_null[[datatype_i]])){
          
          comb_nullmodel <- c(comb_nullmodel,
                              filter(fit_wintr_null[[datatype_i]][[i]],
                                     response == res)[ , paste0("formi_split.b.", InIdIh)])
          
        }
        
        max_y <- max(density(comb_nullmodel)$y, 
                     p$data$y)
        
        density_null <- data.frame(density(comb_nullmodel)) %>% 
          filter(y > y_min_thresh)
        
        
        
        hdis_core <- rbind(nullmodel = hdi(comb_nullmodel, ci = 1 - sqrt(0.05)), 
                           observed = hdi(stats_b_wintr[[datatype_i]][[res]][["formi_split"]]$b[, InIdIh], ci = 1 - sqrt(0.05))) %>% 
          rownames_to_column("data") %>% 
          mutate(data = factor(data, levels = .$data))
        
        hdis_outer <- rbind(nullmodel = hdi(comb_nullmodel, ci = .95), 
                            observed = hdi(stats_b_wintr[[datatype_i]][[res]][["formi_split"]]$b[, InIdIh], ci = .95)) %>% 
          rownames_to_column("data") %>% 
          mutate(data = factor(data, levels = .$data))
        
        est <- rbind(nullmodel = map_estimate(comb_nullmodel), 
                     observed = map_estimate(stats_b_wintr[[datatype_i]][[res]][["formi_split"]]$b[, InIdIh])) %>% 
          rownames_to_column("data") %>% 
          mutate(data = factor(data, levels = .$data))
        
        p <-
          p +
          geom_vline(xintercept = 0, lty = "dashed") +
          geom_line(data = density_null %>% 
                      mutate(y = y / (.75 * max_y),
                             y = y + 1),
                    col = "#7570b3", size = .75) +
          geom_line(data = density_obs %>% 
                      mutate(y = y / (.75 * max_y),
                             y = y + 1),
                    col = "#d95f02", size = .75) +
          geom_segment(data = hdis_outer, aes(x = CI_low, xend = CI_high, 
                                              y = as.numeric(data)/3 + .1, yend = as.numeric(data)/3 + .1,
                                              col = data), size = .3) +
          geom_segment(data = hdis_core, aes(x = CI_low, xend = CI_high, 
                                             y = as.numeric(data)/3 + .1, yend = as.numeric(data)/3 + .1,
                                             col = data), size = 1.25) +
          geom_point(data = est, aes(x = MAP, y = as.numeric(data)/3 + .1, 
                                     col = data), size = 1.75) +
          scale_color_manual(values = c(nullmodel = "#7570b3", 
                                        observed = "#d95f02"),
                             guide = F) + 
          theme(axis.title = element_blank(),
                axis.line.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank()) +
          scale_y_continuous(limits = c(0.2 , 3))

        posterior_plots_wintr_corr[[res]][[InIdIh_names[InIdIh]]][[datatype_i]] <- p
        
        range_x_wintr_corr[[datatype_i]][[res]][1] <- min(density_null$x,
                                                          density_obs$x, 
                                                          range_x_wintr_corr[[datatype_i]][[res]][1])
        range_x_wintr_corr[[datatype_i]][[res]][2] <- max(density_null$x,
                                                          density_obs$x,
                                                          range_x_wintr_corr[[datatype_i]][[res]][2])
        
      }
      
    }
  }
}


range_x_wintr_corr_sel <- list()
for (res in responses){
  range_x_wintr_corr_sel[[res]] <-
    c(min(sapply(range_x_wintr_corr[!names(range_x_wintr_corr) %in% 
                                      c("all.plants", "no.poly", "pres.abs")],
                 function(x) x[[res]][1])),
      max(sapply(range_x_wintr_corr[!names(range_x_wintr_corr) %in% 
                                      c("all.plants", "no.poly", "pres.abs")],
                 function(x) x[[res]][2])))
}

f.plotadjust <- function(x){
  x + theme(plot.background = element_rect(fill = "transparent", colour = NA),
            panel.background = element_rect(fill = "transparent", colour = NA),
            axis.text.x = element_blank(),
            plot.margin = margin(0, 0, 0, 0)) +
    scale_x_continuous(limits = range_x_wintr_corr_sel[[res]])
}

wintr_corrplots <- list()
for (res in responses){
  for (var in c("ForMI", "Inonat", "Idwcut", "Iharv")){

    row1 <- plot_grid(plotlist = lapply(list(ggplot() + ggtitle(var) + 
                                               theme(plot.title = element_text(size = 25),
                                                     axis.line = element_blank(),
                                                     axis.ticks = element_blank(),
                                                     axis.text = element_blank()), 
                                             NULL, 
                                             posterior_plots_wintr_corr[[res]][[var]]$fw, 
                                             NULL, 
                                             NULL),
                                        f.plotadjust),
                      nrow = 1)
    row2 <- plot_grid(plotlist = lapply(list(NULL, 
                                             posterior_plots_wintr_corr[[res]][[var]]$fw.bm, 
                                             posterior_plots_wintr_corr[[res]][[var]]$fw.pw.bm,  
                                             posterior_plots_wintr_corr[[res]][[var]]$fw.pw,  
                                             posterior_plots_wintr_corr[[res]][[var]]$fw.mr),
                                        f.plotadjust),
                      nrow = 1)
    row3 <- plot_grid(plotlist = lapply(list(posterior_plots_wintr_corr[[res]][[var]]$bm,
                                             posterior_plots_wintr_corr[[res]][[var]]$fw.mr.bm,
                                             posterior_plots_wintr_corr[[res]][[var]]$fw.pw.mr.bm, 
                                             posterior_plots_wintr_corr[[res]][[var]]$fw.pw.mr,
                                             posterior_plots_wintr_corr[[res]][[var]]$pw),
                                        f.plotadjust),
                      nrow = 1)
    row4 <- plot_grid(plotlist = lapply(list(NULL,
                                             posterior_plots_wintr_corr[[res]][[var]]$mr.bm,
                                             posterior_plots_wintr_corr[[res]][[var]]$pw.mr.bm,
                                             posterior_plots_wintr_corr[[res]][[var]]$pw.mr,
                                             NULL),
                                        f.plotadjust),
                      nrow = 1)
    row5 <- plot_grid(plotlist = lapply(list(posterior_plots_wintr_corr[[res]][[var]]$raw,
                                             NULL,
                                             posterior_plots_wintr_corr[[res]][[var]]$mr,
                                             posterior_plots_wintr_corr[[res]][[var]]$pw.bm,
                                             NULL), 
                                        f.plotadjust),
                      nrow = 1)
    row6 <- plot_grid(NULL, NULL,
                      ggplot() + geom_blank() +
                        scale_x_continuous(limits = range_x_wintr_corr_sel[[res]]) +
                        theme(axis.line.x = element_blank(),
                              axis.ticks.x = element_blank(),
                              axis.text.x = element_text(angle = 90, size = 16,
                                                         hjust = 1, vjust = .5),
                              plot.background = element_rect(fill = "transparent", colour = NA),
                              panel.background = element_rect(fill = "transparent", colour = NA),
                              plot.margin = margin(0, 0, 0, 0)), 
                      NULL, NULL,
                      nrow = 1)
    
    wintr_corrplots[[res]][[var]] <- plot_grid(row1, row2, row3, row4, row5, row6, ncol = 1,
                                               rel_heights = c(1, 1, 1, 1, 1, .5))
  }
}

separator <- ggplot() + 
  annotate(geom = "segment", x = 1, xend = 1, y = -Inf, yend = Inf, lty = 3) +
  scale_y_continuous(limits = c(-1, 1)) +
  theme_nothing()

for (res in responses){
  p <- plot_grid(
    plot_grid(wintr_corrplots[[res]]$ForMI, separator, wintr_corrplots[[res]]$Inonat,
              nrow = 1, rel_widths = c(20, 1, 20)),
    
    ggplot() + 
      annotate(geom = "segment", x = Inf, xend = -Inf, y = 1, yend = 1, lty = 3) +
      scale_x_continuous(limits = c(-1, 1)) +
      theme_nothing(),
    
    plot_grid(wintr_corrplots[[res]]$Iharv, separator, wintr_corrplots[[res]]$Idwcut,
              nrow = 1, rel_widths = c(20, 1, 20)),
    
    ncol = 1, rel_heights = c(20, 1, 20)
  )

  plot_grid(plotlist = wintr_corrplots[[res]], nrow = 2, ncol = 2)

}


#  ACROSS BOTH ECOSYSTEMS -------------------------------- #####################
################################################################################.

# connectance and network size (analysis) ######################################
################################################################################.

# define function to get connectance and network size
f.nw.metrics.basic <- function(pl_i, corr_i){
  assign("int_global_wintr_target", 
         eval(parse(text  = paste0("int_global_wintr_", 
                                   gsub("\\.", "_", corr_i)))))
  assign("int_global_swnet_target", 
         eval(parse(text  = paste0("int_global_swnet_", 
                                   gsub("\\.", "_", corr_i)))))

  
  int_matrix_plot <- int_global_wintr_target %>% 
    bind_rows(int_global_swnet_target) %>% 
    filter(PlotID == pl_i) %>% 
    select(Host_Name_std, Herb_Name_std, interaction) %>% 
    spread(Herb_Name_std, interaction, fill = 0) %>% 
    column_to_rownames("Host_Name_std")
  
  con_i <- networklevel(int_matrix_plot, index = "connectance")
  ricp_i <- nrow(int_matrix_plot)
  rich_i <- ncol(int_matrix_plot)
  nwsize_i <- ricp_i * rich_i
  
  data.frame(PlotID = pl_i,
             con = con_i,
             ricp = ricp_i,
             rich = rich_i,
             nwsize = nwsize_i,
             row.names = NULL,
             stringsAsFactors = F)
}
  

# apply the function
nw.metrics.basic$raw <- mclapply(c(rownames(wintr_ct), rownames(swnet_ct)), 
                           f.nw.metrics.basic, "raw") %>% 
  do.call(rbind, .)
nw.metrics.basic$no.poly <- mclapply(c(rownames(wintr_ct), rownames(swnet_ct)), 
                                 f.nw.metrics.basic, "no.poly") %>% 
  do.call(rbind, .)
nw.metrics.basic$all.plants <- mclapply(c(rownames(wintr_ct), rownames(swnet_ct)), 
                                 f.nw.metrics.basic, "all.plants") %>% 
  do.call(rbind, .)
nw.metrics.basic$pres.abs <- mclapply(c(rownames(wintr_ct), rownames(swnet_ct)), 
                                        f.nw.metrics.basic, "pres.abs") %>% 
  do.call(rbind, .)


# stats ------------------------------------------------------------------------.

nw.metrics.basic <- lapply(nw.metrics.basic, function(x) x %>% 
                             mutate(nwsize_sqrt = sqrt(nwsize)) %>% 
                             left_join(lui.glob %>% bind_rows(formi), 
                                       by = "PlotID")
stats_b_basic <- list()

# forest models ----------------------------------------------------------------.

f.hiermod.basic.wintr <- function(res, corr_i){  
  data <- nw.metrics.basic[[corr_i]] %>% 
    filter(!is.na(ForMI_z))
  
  if (res %in% c("con", "nwsize_sqrt")){
    data$response <- (data[, res] - mean(data[, res]))/sd(data[, res])
    
    
    data.comb <- list(J = length(unique(data$Exploratory)),
                      N = nrow(data),
                      Exploratory = as.integer(as.factor(data$Exploratory)),
                      ForMI_z = as.vector(data$ForMI_z),
                      y = as.vector(data$response))
    
    
    fit.comb <- sampling(object = stan_formi_mod, data = data.comb,
                         chains = 4, warmup = 500, iter = 5500, control = list(adapt_delta = 0.99))
    
    data.split <- list(J = length(unique(data$Exploratory)),
                       N = nrow(data),
                       Exploratory = as.integer(as.factor(data$Exploratory)),
                       Inonat_z = as.vector(data$Inonat_z),
                       Idwcut_z = as.vector(data$Idwcut_z),
                       Iharv_z = as.vector(data$Iharv_z),
                       y = as.vector(data$response))
    
    fit.split <- sampling(object = stan_formisplit_mod, data = data.split,
                          chains = 4, warmup = 500, iter = 5500, control = list(adapt_delta = 0.99))
    
  } else {
    data$response <- data[, res]
    
    data.comb <- list(J = length(unique(data$Exploratory)),
                      N = nrow(data),
                      Exploratory = as.integer(as.factor(data$Exploratory)),
                      ForMI_z = as.vector(data$ForMI_z),
                      y = as.vector(data$response))
    
    
    fit.comb <- sampling(object = stan_formi_mod_poiss, data = data.comb,
                         chains = 4, warmup = 500, iter = 5500, control = list(adapt_delta = 0.99))
    
    data.split <- list(J = length(unique(data$Exploratory)),
                       N = nrow(data),
                       Exploratory = as.integer(as.factor(data$Exploratory)),
                       Inonat_z = as.vector(data$Inonat_z),
                       Idwcut_z = as.vector(data$Idwcut_z),
                       Iharv_z = as.vector(data$Iharv_z),
                       y = as.vector(data$response))
    
    fit.split <- sampling(object = stan_formisplit_mod_poiss, data = data.split,
                          chains = 4, warmup = 500, iter = 5500, control = list(adapt_delta = 0.99))
  }
  
  list(comb = extract(fit.comb)[c("a", "b")],
       split = extract(fit.split)[c("a", "b")])
  
}

stats_b_basic$raw[["forests"]] <- lapply(c("con", "ricp", "rich", "nwsize_sqrt"), 
                                     f.hiermod.basic.wintr, "raw")
names(stats_b_basic$raw[["forests"]]) <- c("con", "ricp", "rich", "nwsize_sqrt")

# grasslands models ------------------------------------------------------------.

f.hiermod.basic.swnet <- function(res, corr_i){  
  data <- nw.metrics.basic[[corr_i]] %>% 
    filter(!is.na(LUI_z))
  
  
  if (res %in% c("con", "nwsize_sqrt")){
    data$response <- (data[, res] - mean(data[, res]))/sd(data[, res])
    
    
    data.comb <- list(J = length(unique(data$Exploratory)),
                      N = nrow(data),
                      Exploratory = as.integer(as.factor(data$Exploratory)),
                      LUI_z = as.vector(data$LUI_z),
                      y = as.vector(data$response))
    
    
    fit.comb <- sampling(object = stan_lui_mod, data = data.comb,
                         chains = 4, warmup = 500, iter = 5500, control = list(adapt_delta = 0.99))
    
    data.split <- list(J = length(unique(data$Exploratory)),
                       N = nrow(data),
                       Exploratory = as.integer(as.factor(data$Exploratory)),
                       M_z = as.vector(data$M_z),
                       F_z = as.vector(data$F_z),
                       G_z = as.vector(data$G_z),
                       y = as.vector(data$response))
    
    fit.split <- sampling(object = stan_luisplit_mod, data = data.split,
                          chains = 4, warmup = 500, iter = 5500, control = list(adapt_delta = 0.99))
    
  } else {
    data$response <- data[, res]
    
    data.comb <- list(J = length(unique(data$Exploratory)),
                      N = nrow(data),
                      Exploratory = as.integer(as.factor(data$Exploratory)),
                      LUI_z = as.vector(data$LUI_z),
                      y = as.vector(data$response))
    
    
    fit.comb <- sampling(object = stan_lui_mod_poiss, data = data.comb,
                         chains = 4, warmup = 500, iter = 5500, control = list(adapt_delta = 0.99))
    
    data.split <- list(J = length(unique(data$Exploratory)),
                       N = nrow(data),
                       Exploratory = as.integer(as.factor(data$Exploratory)),
                       M_z = as.vector(data$M_z),
                       F_z = as.vector(data$F_z),
                       G_z = as.vector(data$G_z),
                       y = as.vector(data$response))
    
    fit.split <- sampling(object = stan_luisplit_mod_poiss, data = data.split,
                          chains = 4, warmup = 500, iter = 5500, control = list(adapt_delta = 0.99))
  }
  
  list(comb = extract(fit.comb)[c("a", "b")],
       split = extract(fit.split)[c("a", "b")])
  
}

stats_b_basic$raw[["grasslands"]] <- lapply(c("con", "ricp", "rich", "nwsize_sqrt"), 
                                        f.hiermod.basic.swnet, "raw")
names(stats_b_basic$raw[["grasslands"]]) <- c("con", "ricp", "rich", "nwsize_sqrt")


# feeding specialisation #######################################################
################################################################################.

spec_cat_f <- int_global_wintr_raw %>% 
  left_join(plant_agg_f %>% select(-c(PlotID, Cover)) %>% distinct(),
            by = "Host_Name_std") %>% 
  group_by(Herb_Name_std) %>% 
  summarise(n_spec = length(unique(Host_Species)),
            n_gen = length(unique(Host_Genus)),
            n_fam = length(unique(Host_Family))) %>% 
  ungroup() %>% 
  mutate(specialisation = ifelse(n_spec == 1 & n_gen == 1, "mono_spec",
                                 ifelse(n_gen == 1, "mono_gen",
                                        ifelse(n_fam == 1, "mono_fam",
                                               "poly"))))



spec_cat_g <- int_global_swnet_raw %>% 
  left_join(plant_agg_g %>% select(-c(PlotID, Cover)) %>% distinct(),
            by = "Host_Name_std") %>% 
  group_by(Herb_Name_std) %>% 
  summarise(n_spec = length(unique(Host_Species)),
            n_gen = length(unique(Host_Genus)),
            n_fam = length(unique(Host_Family))) %>% 
  ungroup() %>% 
  mutate(specialisation = ifelse(n_spec == 1 & n_gen == 1, "mono_spec",
                                 ifelse(n_gen == 1, "mono_gen",
                                        ifelse(n_fam == 1, "mono_fam",
                                               "poly"))))


feedspec.overall <- wintr_agg %>% 
  inner_join(spec_cat_f, by = "Herb_Name_std") %>% 
  mutate(system = "Forests") %>% 
  bind_rows(swnet_agg %>% 
              inner_join(spec_cat_g, by = "Herb_Name_std") %>% 
              mutate(system = "Grasslands")) %>% 
  group_by(system, PlotID, specialisation) %>% 
  summarise(NumberAdults = sum(NumberAdults),
            NumberSpecies = length(unique(Herb_Name_std))) %>% 
  ungroup() %>% 
  group_by(system, PlotID) %>% 
  mutate(spec_perc_abund = NumberAdults / sum(NumberAdults) * 100,
         spec_perc_rich = NumberSpecies / sum(NumberSpecies) * 100) %>% 
  ungroup() %>% 
  select(system, PlotID, specialisation, spec_perc_abund, spec_perc_rich) %>% 
  gather(type, value, -c(system, PlotID, specialisation)) %>% 
  spread(specialisation, value, fill = 0) %>% 
  gather(specialisation, value, -c(system, PlotID, type)) %>% 
  group_by(system, specialisation, type) %>% 
  summarise(value = mean(value)) %>% 
  ungroup() %>% 
  mutate(Herb_Group = "overall")


feedspec.group <- wintr_agg %>% 
  inner_join(spec_cat_f, by = "Herb_Name_std") %>% 
  mutate(system = "Forests") %>% 
  bind_rows(swnet_agg %>% 
              inner_join(spec_cat_g, by = "Herb_Name_std") %>% 
              mutate(system = "Grasslands")) %>% 
  group_by(system, PlotID, specialisation, Herb_Group) %>% 
  summarise(NumberAdults = sum(NumberAdults),
            NumberSpecies = length(unique(Herb_Name_std))) %>% 
  ungroup() %>% 
  group_by(system, PlotID) %>% 
  mutate(spec_perc_abund = NumberAdults / sum(NumberAdults) * 100,
         spec_perc_rich = NumberSpecies / sum(NumberSpecies) * 100) %>% 
  ungroup() %>% 
  select(system, PlotID, specialisation, spec_perc_abund, spec_perc_rich, Herb_Group) %>% 
  gather(type, value, -c(system, PlotID, specialisation, Herb_Group)) %>% 
  spread(specialisation, value, fill = 0) %>% 
  gather(specialisation, value, -c(system, PlotID, type, Herb_Group)) %>% 
  group_by(system, specialisation, type, Herb_Group) %>% 
  summarise(value = mean(value)) %>% 
  ungroup()

p.spec <- feedspec.overall %>% 
  bind_rows(feedspec.group) %>% 
  mutate(specialisation = factor(specialisation,
                                 levels = c("mono_spec", "mono_gen", "mono_fam",
                                            "poly"),
                                 labels = c("Species", "Genus", "Family",
                                            "Polyphagous")),
         Herb_Group = factor(Herb_Group, levels = c("overall", "Auchenorrhyncha", 
                                                    "Heteroptera", "Coleoptera", 
                                                    "Orthoptera")),
         type = factor(type, label = c("Abundance", "Richness"))) %>% 
  mutate(system = factor(system, levels = c("Grasslands", "Forests"))) %>% 
  ggplot(aes(x = Herb_Group, y = value, fill = specialisation)) +
  geom_col() +
  geom_vline(xintercept = 1.5, lty = 3) +
  scale_y_continuous(labels = c(paste0(seq(0,100,25), "%"))) +
  scale_fill_manual(values = c("#D2AB99", "#BDBEA9", "#8DB38B",
                               "#04724D"),
                    name = "Specialisation") +
  facet_grid(type ~ system) +
  xlab("Herbivore group") +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))


# taxonomic resolution of information ------------------------------------------.

res_info <- int.db %>% 
  
  group_by(Herb_Name_std) %>% 
  summarise(rank_high = ifelse(any(grepl("polyphagous", Host_Genus_Species)), 
                               "polyphagous", 
                               as.character(min(Host_Rank_red)))) %>% 
  mutate(rank_high = ifelse(rank_high %in% c("Subspecies"), "Species", rank_high),
         rank_high = ifelse(rank_high %in% c("Tribus"), "Family", rank_high),
         rank_high = ifelse(rank_high %in% c("Subdivision", "Cladus1", "Class", 
                                             "Order", "polyphagous"), 
                            "Higher Order", rank_high))


resrank.overall <-
  wintr_agg %>% 
  inner_join(res_info, by = "Herb_Name_std") %>% 
  mutate(system = "Forests") %>% 
  bind_rows(swnet_agg %>% 
              inner_join(res_info, by = "Herb_Name_std") %>% 
              mutate(system = "Grasslands")) %>% 
  group_by(system, PlotID, rank_high) %>% 
  summarise(NumberAdults = sum(NumberAdults),
            NumberSpecies = length(unique(Herb_Name_std))) %>% 
  ungroup() %>% 
  group_by(system, PlotID) %>% 
  mutate(spec_perc_abund = NumberAdults / sum(NumberAdults) * 100,
         spec_perc_rich = NumberSpecies / sum(NumberSpecies) * 100) %>% 
  ungroup() %>% 
  select(system, PlotID, rank_high, spec_perc_abund, spec_perc_rich) %>% 
  gather(type, value, -c(system, PlotID, rank_high)) %>% 
  spread(rank_high, value, fill = 0) %>% 
  gather(rank_high, value, -c(system, PlotID, type)) %>% 
  group_by(system, rank_high, type) %>% 
  summarise(value = mean(value)) %>% 
  ungroup() %>% 
  mutate(Herb_Group = "overall")


resrank.group <- wintr_agg %>% 
  inner_join(res_info, by = "Herb_Name_std") %>% 
  mutate(system = "Forests") %>% 
  bind_rows(swnet_agg %>% 
              inner_join(res_info, by = "Herb_Name_std") %>% 
              mutate(system = "Grasslands")) %>% 
  group_by(system, PlotID, rank_high, Herb_Group) %>% 
  summarise(NumberAdults = sum(NumberAdults),
            NumberSpecies = length(unique(Herb_Name_std))) %>% 
  ungroup() %>% 
  group_by(system, PlotID) %>% 
  mutate(spec_perc_abund = NumberAdults / sum(NumberAdults) * 100,
         spec_perc_rich = NumberSpecies / sum(NumberSpecies) * 100) %>% 
  ungroup() %>% 
  select(system, PlotID, rank_high, spec_perc_abund, spec_perc_rich, Herb_Group) %>% 
  gather(type, value, -c(system, PlotID, rank_high, Herb_Group)) %>% 
  spread(rank_high, value, fill = 0) %>% 
  gather(rank_high, value, -c(system, PlotID, type, Herb_Group)) %>% 
  group_by(system, rank_high, type, Herb_Group) %>% 
  summarise(value = mean(value)) %>% 
  ungroup()


p.res <- resrank.overall %>% 
  bind_rows(resrank.group) %>% 
  mutate(rank_high = factor(rank_high,
                            levels = c("Species", "Genus", "Family",
                                       "Higher Order")),
         Herb_Group = factor(Herb_Group, levels = c("overall", "Auchenorrhyncha", 
                                                    "Heteroptera", "Coleoptera", 
                                                    "Orthoptera")),
         type = factor(type, label = c("Abundance", "Richness"))) %>% 
  mutate(system = factor(system, levels = c("Grasslands", "Forests"))) %>% 
  ggplot(aes(x = Herb_Group, y = value, fill = rank_high)) +
  geom_col() +
  geom_vline(xintercept = 1.5, lty = 3) +
  # scale_x_discrete(labels = c("Abundance", "Richness")) +
  scale_y_continuous(labels = c(paste0(seq(0,100,25), "%"))) +
  scale_fill_manual(values = c("#CC7171", "#B5B48D", "#97D392",
                               "#1B6022"),
                    name = "Taxonomic resolution") +
  facet_grid(type ~ system) +
  xlab("Herbivore group") +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))


rel_widths <- c(7, 1, 7.8)

plot_grid(
plot_grid(ggplot_labeller(expression(bold("a")), 0, label.size, NA),
          empty,
          empty,
          ggplot_labeller(expression(bold("b")), 0, label.size, NA),
          empty,
          nrow = 1, rel_widths = c(.5, rel_widths[1] - .5, rel_widths[2],
                                   .5, rel_widths[3] - .5)),
plot_grid(p.spec, 
          
          ggplot() + 
            annotate(geom = "segment", x = 1, xend = 1, y = -Inf, yend = Inf, lty = 3) +
            scale_y_continuous(limits = c(-1, 1)) +
            theme_nothing(),
          
          p.res, nrow = 1, rel_widths = rel_widths),
ncol = 1, rel_heights = c(1, 20))


# plant community composition  #################################################
################################################################################.

plant_comp_f <- plant_agg_f %>% 
  filter(Cover > 0) %>% 
  mutate(Host_Group1 = ifelse(Host_Group1 == "clubmoss", "forb", Host_Group1)) %>% 
  group_by(PlotID, Host_Group1) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  select(PlotID, Host_Group1, n) %>% 
  spread(Host_Group1, n, fill = 0) %>% 
  gather(Host_Group1, n, -PlotID) %>% 
  full_join(plant_agg_f %>% 
              filter(Cover > 0) %>% 
              mutate(Host_Group1 = ifelse(Host_Group1 == "clubmoss", "forb", Host_Group1)) %>% 
              group_by(PlotID, Host_Group1) %>% 
              summarise(Cover = sum(Cover)) %>% 
              ungroup() %>% 
              select(PlotID, Host_Group1, Cover) %>% 
              spread(Host_Group1, Cover, fill = 0) %>% 
              gather(Host_Group1, Cover, -PlotID),
            by = c("PlotID", "Host_Group1")) %>% 
  group_by(Host_Group1) %>% 
  summarise(n = mean(n),
            Cover = mean(Cover))
  
plant_comp_g <- plant_agg_g %>% 
  filter(Cover > 0) %>% 
  mutate(Host_Group1 = ifelse(Host_Group1 == "clubmoss", "forb", Host_Group1)) %>% 
  group_by(PlotID, Host_Group1) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  select(PlotID, Host_Group1, n) %>% 
  spread(Host_Group1, n, fill = 0) %>% 
  gather(Host_Group1, n, -PlotID) %>% 
  full_join(plant_agg_g %>% 
              filter(Cover > 0) %>% 
              mutate(Host_Group1 = ifelse(Host_Group1 == "clubmoss", "forb", Host_Group1)) %>% 
              group_by(PlotID, Host_Group1) %>% 
              summarise(Cover = sum(Cover)) %>% 
              ungroup() %>% 
              select(PlotID, Host_Group1, Cover) %>% 
              spread(Host_Group1, Cover, fill = 0) %>% 
              gather(Host_Group1, Cover, -PlotID),
            by = c("PlotID", "Host_Group1")) %>% 
  group_by(Host_Group1) %>% 
  summarise(n = mean(n),
            Cover = mean(Cover))

colors.plantgroups <- c(fern = "#E1BC29",
                        forb = "#4D9DE0",
                        grass = "#A3ED5A",
                        shrub = "#E06655",
                        tree = "#A02424") 

labels.plantgroups <- c(fern = "Ferns",
                        forb = "Forbs",
                        grass = "Grasses",
                        shrub = "Shrubs",
                        tree = "Trees") 


plant_comp_f %>% 
  mutate(system = "Forests") %>% 
  bind_rows(plant_comp_g %>% 
              mutate(system = "Grasslands")) %>% 
  gather(var, value, -c(Host_Group1, system)) %>% 
  mutate(system = factor(system, levels = c("Grasslands", "Forests"))) %>% 
  ggplot(aes(x = system, y = value, fill = Host_Group1)) +
  geom_col() +
  scale_fill_manual(values = colors.plantgroups,
                    labels = labels.plantgroups,
                    name = "Plant group") +
  theme(axis.title = element_blank()) +
  facet_wrap(var ~ ., scales = "free_y",
             strip.position = "left", 
             labeller = as_labeller(c(n = "Mean number of species", 
                                      Cover = "Mean cover [%]") )) +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 11))

# overview of network metrics across systems  ##################################
################################################################################.

# the following code produces an overview for all
# network metrics for a set of land-use scenarios in the two ecosystems

# for "historic" reasons, basic network metrics (species richness, network size,
# connectance) were analyzed separately from other metrics (modularity, 
# nestedness, robustness). Thus, the code is overly complex.


png.g <- image_read("grasslands.png")
image_g <- image_ggplot(png.g)
png.f <- image_read("forests.png")
image_f <- image_ggplot(png.f)

corr_i <- "raw" # no sensitivity corrections

# ------------------------------------------------------------------------------.
# data for basic metrics -------------------------------------------------------.
# ------------------------------------------------------------------------------.

map.basic <- data.frame()
hdi.basic <- data.frame()
int.basic.map <- list()
int.basic.density <- list()

for (sys_i in names(stats_b_basic[[corr_i]])){
  for (res_i in names(stats_b_basic[[corr_i]][[sys_i]])){
    if (res_i %in% c("ricp", "rich")){
      func_inv <- exp # inverse function: exponential
    } else {
      func_inv <- function(x) { # inverse function: change back to original scale
        if (sys_i == "forests") data <- filter(nw.metrics.basic[[corr_i]], !is.na(ForMI))
        if (sys_i == "grasslands") data <- filter(nw.metrics.basic[[corr_i]], !is.na(LUI))
        x * sd(data[, res_i]) + mean(data[, res_i])
      } 
    }
    
    
    # Forest scenarios -----------------------------------------------------------.
    
    if (sys_i == "forests"){
      pred <- func_inv(
        (rowMeans(stats_b_basic[[corr_i]][[sys_i]][[res_i]]$split$a) + 
           stats_b_basic[[corr_i]][[sys_i]][[res_i]]$split$b[, 1] * max(formi$Inonat_z) +
           stats_b_basic[[corr_i]][[sys_i]][[res_i]]$split$b[, 2] * min(formi$Idwcut_z) +
           stats_b_basic[[corr_i]][[sys_i]][[res_i]]$split$b[, 3] * min(formi$Iharv_z))
      )
      map.basic <- data.frame(sys = sys_i, res = res_i, scenario = "Inonat_max",
                              MAP = as.numeric(map_estimate(pred)),
                              stringsAsFactors = F) %>% 
        bind_rows(map.basic, .)
      
      
      hdi.basic <- data.frame(sys = sys_i, res = res_i, scenario = "Inonat_max",
                              rbind(hdi(pred, ci = 1 - sqrt(0.05)),
                                    hdi(pred, ci = .95)),
                              hdi.type = c("core", "outer"),
                              stringsAsFactors = F) %>% 
        bind_rows(hdi.basic, .)
      # ------------------------------------------------------------------------.
      pred <- func_inv(
        (rowMeans(stats_b_basic[[corr_i]][[sys_i]][[res_i]]$split$a) + 
           stats_b_basic[[corr_i]][[sys_i]][[res_i]]$split$b[, 1] * min(formi$Inonat_z) +
           stats_b_basic[[corr_i]][[sys_i]][[res_i]]$split$b[, 2] * max(formi$Idwcut_z) +
           stats_b_basic[[corr_i]][[sys_i]][[res_i]]$split$b[, 3] * min(formi$Iharv_z))
      )
      map.basic <- data.frame(sys = sys_i, res = res_i, scenario = "Idwcut_max",
                              MAP = as.numeric(map_estimate(pred)),
                              stringsAsFactors = F) %>% 
        bind_rows(map.basic, .)
      
      
      hdi.basic <- data.frame(sys = sys_i, res = res_i, scenario = "Idwcut_max",
                              rbind(hdi(pred, ci = 1 - sqrt(0.05)),
                                    hdi(pred, ci = .95)),
                              hdi.type = c("core", "outer"),
                              stringsAsFactors = F) %>% 
        bind_rows(hdi.basic, .)
      # ------------------------------------------------------------------------.
      pred <- func_inv(
        (rowMeans(stats_b_basic[[corr_i]][[sys_i]][[res_i]]$split$a) + 
           stats_b_basic[[corr_i]][[sys_i]][[res_i]]$split$b[, 1] * min(formi$Inonat_z) +
           stats_b_basic[[corr_i]][[sys_i]][[res_i]]$split$b[, 2] * min(formi$Idwcut_z) +
           stats_b_basic[[corr_i]][[sys_i]][[res_i]]$split$b[, 3] * max(formi$Iharv_z))
      )
      map.basic <- data.frame(sys = sys_i, res = res_i, scenario = "Iharv_max",
                              MAP = as.numeric(map_estimate(pred)),
                              stringsAsFactors = F) %>% 
        bind_rows(map.basic, .)
      
      
      hdi.basic <- data.frame(sys = sys_i, res = res_i, scenario = "Iharv_max",
                              rbind(hdi(pred, ci = 1 - sqrt(0.05)),
                                    hdi(pred, ci = .95)),
                              hdi.type = c("core", "outer"),
                              stringsAsFactors = F) %>% 
        bind_rows(hdi.basic, .)
      # ------------------------------------------------------------------------.
      pred <- func_inv(
        (rowMeans(stats_b_basic[[corr_i]][[sys_i]][[res_i]]$comb$a) + 
           stats_b_basic[[corr_i]][[sys_i]][[res_i]]$comb$b * min(formi$ForMI_z))
      )
      map.basic <- data.frame(sys = sys_i, res = res_i, scenario = "ForMI_min",
                              MAP = as.numeric(map_estimate(pred)),
                              stringsAsFactors = F) %>% 
        bind_rows(map.basic, .)
      
      
      hdi.basic <- data.frame(sys = sys_i, res = res_i, scenario = "ForMI_min",
                              rbind(hdi(pred, ci = 1 - sqrt(0.05)),
                                    hdi(pred, ci = .95)),
                              hdi.type = c("core", "outer"),
                              stringsAsFactors = F) %>% 
        bind_rows(hdi.basic, .)
      # ------------------------------------------------------------------------.
      pred <- func_inv(
        (rowMeans(stats_b_basic[[corr_i]][[sys_i]][[res_i]]$comb$a) + 
           stats_b_basic[[corr_i]][[sys_i]][[res_i]]$comb$b * max(formi$ForMI_z))
      )
      map.basic <- data.frame(sys = sys_i, res = res_i, scenario = "ForMI_max",
                              MAP = as.numeric(map_estimate(pred)),
                              stringsAsFactors = F) %>% 
        bind_rows(map.basic, .)
      
      
      hdi.basic <- data.frame(sys = sys_i, res = res_i, scenario = "ForMI_max",
                              rbind(hdi(pred, ci = 1 - sqrt(0.05)),
                                    hdi(pred, ci = .95)),
                              hdi.type = c("core", "outer"),
                              stringsAsFactors = F) %>% 
        bind_rows(hdi.basic, .)
    }
    
    # ------------------------------------------------------------------------.
    # ------------------------------------------------------------------------.
    
    if (sys_i == "grasslands"){
      pred <- func_inv(
        (rowMeans(stats_b_basic[[corr_i]][[sys_i]][[res_i]]$split$a) + 
           stats_b_basic[[corr_i]][[sys_i]][[res_i]]$split$b[, 1] * max(lui.glob$M_z) +
           stats_b_basic[[corr_i]][[sys_i]][[res_i]]$split$b[, 2] * min(lui.glob$F_z) +
           stats_b_basic[[corr_i]][[sys_i]][[res_i]]$split$b[, 3] * min(lui.glob$G_z))
      )
      map.basic <- data.frame(sys = sys_i, res = res_i, scenario = "M_max",
                              MAP = as.numeric(map_estimate(pred)),
                              stringsAsFactors = F) %>% 
        bind_rows(map.basic, .)
      
      
      hdi.basic <- data.frame(sys = sys_i, res = res_i, scenario = "M_max",
                              rbind(hdi(pred, ci = 1 - sqrt(0.05)),
                                    hdi(pred, ci = .95)),
                              hdi.type = c("core", "outer"),
                              stringsAsFactors = F) %>% 
        bind_rows(hdi.basic, .)
      # ------------------------------------------------------------------------.
      pred <- func_inv(
        (rowMeans(stats_b_basic[[corr_i]][[sys_i]][[res_i]]$split$a) + 
           stats_b_basic[[corr_i]][[sys_i]][[res_i]]$split$b[, 1] * min(lui.glob$M_z) +
           stats_b_basic[[corr_i]][[sys_i]][[res_i]]$split$b[, 2] * max(lui.glob$F_z) +
           stats_b_basic[[corr_i]][[sys_i]][[res_i]]$split$b[, 3] * min(lui.glob$G_z))
      )
      map.basic <- data.frame(sys = sys_i, res = res_i, scenario = "F_max",
                              MAP = as.numeric(map_estimate(pred)),
                              stringsAsFactors = F) %>% 
        bind_rows(map.basic, .)
      
      
      hdi.basic <- data.frame(sys = sys_i, res = res_i, scenario = "F_max",
                              rbind(hdi(pred, ci = 1 - sqrt(0.05)),
                                    hdi(pred, ci = .95)),
                              hdi.type = c("core", "outer"),
                              stringsAsFactors = F) %>% 
        bind_rows(hdi.basic, .)
      # ------------------------------------------------------------------------.
      pred <- func_inv(
        (rowMeans(stats_b_basic[[corr_i]][[sys_i]][[res_i]]$split$a) + 
           stats_b_basic[[corr_i]][[sys_i]][[res_i]]$split$b[, 1] * min(lui.glob$M_z) +
           stats_b_basic[[corr_i]][[sys_i]][[res_i]]$split$b[, 2] * min(lui.glob$F_z) +
           stats_b_basic[[corr_i]][[sys_i]][[res_i]]$split$b[, 3] * max(lui.glob$G_z))
      )
      map.basic <- data.frame(sys = sys_i, res = res_i, scenario = "G_max",
                              MAP = as.numeric(map_estimate(pred)),
                              stringsAsFactors = F) %>% 
        bind_rows(map.basic, .)
      
      
      hdi.basic <- data.frame(sys = sys_i, res = res_i, scenario = "G_max",
                              rbind(hdi(pred, ci = 1 - sqrt(0.05)),
                                    hdi(pred, ci = .95)),
                              hdi.type = c("core", "outer"),
                              stringsAsFactors = F) %>% 
        bind_rows(hdi.basic, .)
      # ------------------------------------------------------------------------.
      pred <- func_inv(
        (rowMeans(stats_b_basic[[corr_i]][[sys_i]][[res_i]]$comb$a) + 
           stats_b_basic[[corr_i]][[sys_i]][[res_i]]$comb$b * min(lui.glob$LUI_z))
      )
      map.basic <- data.frame(sys = sys_i, res = res_i, scenario = "LUI_min",
                              MAP = as.numeric(map_estimate(pred)),
                              stringsAsFactors = F) %>% 
        bind_rows(map.basic, .)
      
      
      hdi.basic <- data.frame(sys = sys_i, res = res_i, scenario = "LUI_min",
                              rbind(hdi(pred, ci = 1 - sqrt(0.05)),
                                    hdi(pred, ci = .95)),
                              hdi.type = c("core", "outer"),
                              stringsAsFactors = F) %>% 
        bind_rows(hdi.basic, .)
      # ------------------------------------------------------------------------.
      pred <- func_inv(
        (rowMeans(stats_b_basic[[corr_i]][[sys_i]][[res_i]]$comb$a) + 
           stats_b_basic[[corr_i]][[sys_i]][[res_i]]$comb$b * max(lui.glob$LUI_z))
      )
      map.basic <- data.frame(sys = sys_i, res = res_i, scenario = "LUI_max",
                              MAP = as.numeric(map_estimate(pred)),
                              stringsAsFactors = F) %>% 
        bind_rows(map.basic, .)
      
      
      hdi.basic <- data.frame(sys = sys_i, res = res_i, scenario = "LUI_max",
                              rbind(hdi(pred, ci = 1 - sqrt(0.05)),
                                    hdi(pred, ci = .95)),
                              hdi.type = c("core", "outer"),
                              stringsAsFactors = F) %>% 
        bind_rows(hdi.basic, .)
      # ------------------------------------------------------------------------.
      # ------------------------------------------------------------------------.
      
    }
    # overall intercepts:
    pred <- func_inv(rowMeans(stats_b_basic[[corr_i]][[sys_i]][[res_i]]$comb$a))
    
    int.basic.map[[sys_i]] <- data.frame(sys = sys_i, res = res_i,
                                         MAP = as.numeric(map_estimate(pred)),
                                         stringsAsFactors = F) %>% 
      bind_rows(int.basic.map[[sys_i]], .)
    
    int.basic.density[[sys_i]] <- data.frame(sys = sys_i, res = res_i,
                                             density(pred),
                                             stringsAsFactors = F) %>% 
      bind_rows(int.basic.density[[sys_i]], .)
    
  }
}

# ------------------------------------------------------------------------------.
# data for modularity, nestedness, robustness ----------------------------------.
# ------------------------------------------------------------------------------.


scenarios.obs <- list()
int.obs <- list()
int.null <- list()


for (res in responses){
  
  # Forest scenarios -----------------------------------------------------------.
  
  scenarios.obs[[res]]$Inonat_max <- 
    c((rowMeans(stats_b_wintr[[corr_i]][[res]]$formi_split$a) + 
         stats_b_wintr[[corr_i]][[res]]$formi_split$b[, 1] * max(formi$Inonat_z) +
         stats_b_wintr[[corr_i]][[res]]$formi_split$b[, 2] * min(formi$Idwcut_z) +
         stats_b_wintr[[corr_i]][[res]]$formi_split$b[, 3] * min(formi$Iharv_z)) * 
        sd(nw.metrics.wintr[[corr_i]][[res]]) + mean(nw.metrics.wintr[[corr_i]][[res]]))
  
  scenarios.obs[[res]]$Idwcut_max <- 
    c((rowMeans(stats_b_wintr[[corr_i]][[res]]$formi_split$a) + 
         stats_b_wintr[[corr_i]][[res]]$formi_split$b[, 1] * min(formi$Inonat_z) +
         stats_b_wintr[[corr_i]][[res]]$formi_split$b[, 2] * max(formi$Idwcut_z) +
         stats_b_wintr[[corr_i]][[res]]$formi_split$b[, 3] * min(formi$Iharv_z)) * 
        sd(nw.metrics.wintr[[corr_i]][[res]]) + mean(nw.metrics.wintr[[corr_i]][[res]]))
  
  scenarios.obs[[res]]$Iharv_max <- 
    c((rowMeans(stats_b_wintr[[corr_i]][[res]]$formi_split$a) + 
         stats_b_wintr[[corr_i]][[res]]$formi_split$b[, 1] * min(formi$Inonat_z) +
         stats_b_wintr[[corr_i]][[res]]$formi_split$b[, 2] * min(formi$Idwcut_z) +
         stats_b_wintr[[corr_i]][[res]]$formi_split$b[, 3] * max(formi$Iharv_z)) * 
        sd(nw.metrics.wintr[[corr_i]][[res]]) + mean(nw.metrics.wintr[[corr_i]][[res]]))
  
  scenarios.obs[[res]]$ForMI_min <- 
    c((rowMeans(stats_b_wintr[[corr_i]][[res]]$formi_comb$a) + 
         stats_b_wintr[[corr_i]][[res]]$formi_comb$b * min(formi$ForMI_z)) * 
        sd(nw.metrics.wintr[[corr_i]][[res]]) + mean(nw.metrics.wintr[[corr_i]][[res]]))
  
  scenarios.obs[[res]]$ForMI_max <-
    c((rowMeans(stats_b_wintr[[corr_i]][[res]]$formi_comb$a) + 
         stats_b_wintr[[corr_i]][[res]]$formi_comb$b * max(formi$ForMI_z)) * 
        sd(nw.metrics.wintr[[corr_i]][[res]]) + mean(nw.metrics.wintr[[corr_i]][[res]]))
  
  # Grassland scenarios --------------------------------------------------------.
  
  scenarios.obs[[res]]$M_max <- 
    c((rowMeans(stats_b_swnet[[corr_i]][[res]]$lui_split$a) + 
         stats_b_swnet[[corr_i]][[res]]$lui_split$b[, 1] * max(lui.glob$M_z) +
         stats_b_swnet[[corr_i]][[res]]$lui_split$b[, 2] * min(lui.glob$F_z) +
         stats_b_swnet[[corr_i]][[res]]$lui_split$b[, 3] * min(lui.glob$G_z)) * 
        sd(nw.metrics.swnet[[corr_i]][[res]]) + mean(nw.metrics.swnet[[corr_i]][[res]]))
  
  scenarios.obs[[res]]$F_max <- 
    c((rowMeans(stats_b_swnet[[corr_i]][[res]]$lui_split$a) + 
         stats_b_swnet[[corr_i]][[res]]$lui_split$b[, 1] * min(lui.glob$M_z) +
         stats_b_swnet[[corr_i]][[res]]$lui_split$b[, 2] * max(lui.glob$F_z) +
         stats_b_swnet[[corr_i]][[res]]$lui_split$b[, 3] * min(lui.glob$G_z)) * 
        sd(nw.metrics.swnet[[corr_i]][[res]]) + mean(nw.metrics.swnet[[corr_i]][[res]]))
  
  scenarios.obs[[res]]$G_max <- 
    c(x = (rowMeans(stats_b_swnet[[corr_i]][[res]]$lui_split$a) + 
             stats_b_swnet[[corr_i]][[res]]$lui_split$b[, 1] * min(lui.glob$M_z) +
             stats_b_swnet[[corr_i]][[res]]$lui_split$b[, 2] * min(lui.glob$F_z) +
             stats_b_swnet[[corr_i]][[res]]$lui_split$b[, 3] * max(lui.glob$G_z)) * 
        sd(nw.metrics.swnet[[corr_i]][[res]]) + mean(nw.metrics.swnet[[corr_i]][[res]]))
  
  scenarios.obs[[res]]$LUI_min <- 
    c((rowMeans(stats_b_swnet[[corr_i]][[res]]$lui_comb$a) + 
         stats_b_swnet[[corr_i]][[res]]$lui_comb$b * min(lui.glob$LUI_z)) * 
        sd(nw.metrics.swnet[[corr_i]][[res]]) + mean(nw.metrics.swnet[[corr_i]][[res]]))
  
  scenarios.obs[[res]]$LUI_max <-
    c((rowMeans(stats_b_swnet[[corr_i]][[res]]$lui_comb$a) + 
         stats_b_swnet[[corr_i]][[res]]$lui_comb$b * max(lui.glob$LUI_z)) * 
        sd(nw.metrics.swnet[[corr_i]][[res]]) + mean(nw.metrics.swnet[[corr_i]][[res]]))
  
  
  # overall means --------------------------------------------------------------.
  
  int.obs[["forest"]][[res]] <- 
    rowMeans(stats_b_wintr[[corr_i]][[res]]$formi_comb$a) * 
    sd(nw.metrics.wintr[[corr_i]][[res]]) + 
    mean(nw.metrics.wintr[[corr_i]][[res]])
  
  
  int.obs[["grassland"]][[res]] <- 
    rowMeans(stats_b_swnet[[corr_i]][[res]]$lui_comb$a) * 
    sd(nw.metrics.swnet[[corr_i]][[res]]) + 
    mean(nw.metrics.swnet[[corr_i]][[res]])
  
}

# extract HDI from scenarios:
hdi.scenario.obs <- lapply(scenarios.obs, function(x) 
  lapply(x, function(x) 
    data.frame(rbind(hdi(x, ci = 1 - sqrt(0.05)),
                     hdi(x, ci = .95)),
               hdi.type = c("core", "outer"))))

# extract MAP estimate from scenarios:
est.scenario.obs <- lapply(scenarios.obs, function(x) 
  lapply(x, function(x) 
    data.frame(MAP = as.numeric(map_estimate(x)))))


# extract MAP estimate from observed intercept:
est.obs <- lapply(int.obs, function(x) lapply(x, function(x) data.frame(MAP = as.numeric(map_estimate(x))))) %>% 
  map(function(x) map(x, cbind.data.frame) %>% 
        enframe("response") %>% 
        unnest())

# extract density from observed intercept:
density.obs <- lapply(int.obs, function(x) lapply(x, function(x) data.frame(density(x)))) %>% 
  map(function(x) map(x, cbind.data.frame) %>% 
        enframe("response") %>% 
        unnest())


# set it all together ----------------------------------------------------------.

scenario.set <- data.frame(scenario = c("ForMI_min", "ForMI_max", 
                                        "Inonat_max", "Iharv_max", "Idwcut_max",
                                        "LUI_min", "LUI_max", 
                                        "M_max", "G_max", "F_max"),
                           aggr = c("comb", "comb", "split", "split", "split",
                                    "comb", "comb", "split", "split", "split"),
                           system = c(rep("forest", 5), rep("grassland", 5)),
                           stringsAsFactors = F)

d.hdi.scenario <- hdi.scenario.obs %>% 
  map(function(x) map(x, cbind.data.frame) %>% 
        enframe("scenario") %>% 
        unnest()) %>% 
  map(cbind.data.frame) %>% 
  enframe("response") %>% 
  unnest() %>% 
  left_join(scenario.set, by = "scenario") %>% 
  mutate(scenario = factor(scenario, levels = scenario.set$scenario))

d.est.scenario <- est.scenario.obs %>% 
  map(function(x) map(x, cbind.data.frame) %>% 
        enframe("scenario") %>% 
        unnest()) %>% 
  map(cbind.data.frame) %>% 
  enframe("response") %>% 
  unnest() %>% 
  left_join(scenario.set, by = "scenario") %>% 
  mutate(scenario = factor(scenario, levels = scenario.set$scenario))

# ------------------------------------------------------------------------------.
# combine the two sets of metrics ----------------------------------------------.
# ------------------------------------------------------------------------------.

map.scenarios <- d.est.scenario %>% 
  rename(res = response,
         sys = system) %>% 
  mutate(sys = paste0(sys, "s")) %>% 
  bind_rows(map.basic)

hdi.scenarios <- d.hdi.scenario %>% 
  rename(res = response,
         sys = system) %>% 
  mutate(sys = paste0(sys, "s")) %>% 
  bind_rows(hdi.basic)

int.map <- int.basic.map
int.density <- int.basic.density
for (i in c(1:2)){
  int.map[[i]] <- int.map[[i]] %>% 
    bind_rows(est.obs[[i]] %>% 
                rename(res = response))
  
  int.density[[i]] <- int.density[[i]] %>% 
    bind_rows(density.obs[[i]] %>% 
                rename(res = response))
}

# last edits -------------------------------------------------------------------.

labels.res <- c(ricp = expression(Plant~sp.~richness),
                rich = expression(Herbivore~sp.~richness),
                nwsize_sqrt = expression(sqrt(Network~size)),
                con = expression(Connectance),
                qmod = expression(Modularity),
                nested = expression(Nestedness),
                rob.A = expression(Robustness))

map.scenarios <- map.scenarios %>% 
  mutate(scenario = factor(scenario, levels = c("ForMI_min", "ForMI_max", 
                                                "Inonat_max", "Iharv_max", 
                                                "Idwcut_max",
                                                "LUI_min", "LUI_max", 
                                                "M_max", "G_max", "F_max")),
         res = factor(res, levels = names(labels.res), labels = labels.res))

hdi.scenarios <- hdi.scenarios %>% 
  mutate(scenario = factor(scenario, levels = c("ForMI_min", "ForMI_max", 
                                                "Inonat_max", "Iharv_max", 
                                                "Idwcut_max",
                                                "LUI_min", "LUI_max", 
                                                "M_max", "G_max", "F_max")),
         res = factor(res, levels = names(labels.res), labels = labels.res))

int.density <- lapply(int.density, function(df) df %>% 
                        group_by(res) %>% 
                        mutate(y = y / max(y), # scale to max 1
                               x_diff = median(diff(x))) %>% # add differences for later plotting (width of rects)
                        ungroup() %>% 
                        mutate(res = factor(res, levels = names(labels.res), 
                                            labels = labels.res)))

int.map <- lapply(int.map, function(df) df %>% 
                    mutate(res = factor(res, levels = names(labels.res), 
                                        labels = labels.res)))

# ------------------------------------------------------------------------------.
# plotting ---------------------------------------------------------------------.
# ------------------------------------------------------------------------------.

labels.x <-  c(ForMI_min = expression(ForMI[min]), 
               ForMI_max = expression(ForMI[max]),
               Inonat_max = expression(Inonat[max]), 
               Iharv_max = expression(Iharv[max]),
               Idwcut_max = expression(Idwcut[max]),
               LUI_min = expression(LUI[min]), 
               LUI_max = expression(LUI[max]),
               M_max = expression(Mowing[max]), 
               G_max = expression(Grazing[max]),
               F_max = expression(Fertilization[max]))

p <- map.scenarios %>% 
  mutate(scenario = factor(scenario, levels = rev(levels(scenario)))) %>% 
  ggplot() +
  geom_point(aes(x = MAP, y = scenario), size = 2, col = 1) +
  geom_rect(data = int.density$forests,
            aes(xmin = x+x_diff/2, xmax = x-x_diff/2,
                ymin = 5.5, ymax = Inf, alpha = y), fill = "#d95f02", col = NA) +
  geom_rect(data = int.density$grasslands,
            aes(xmin = x+x_diff/2, xmax = x-x_diff/2,
                ymin = -Inf, ymax = 5.5, alpha = y), fill = "#d95f02", col = NA) +
  geom_hline(yintercept = seq(1.5, 9.5, 1), lty = 3) +
  geom_hline(yintercept = 5.5, size = 1) +
  geom_segment(data = int.map$forests,
               aes(y = 5.5, yend = Inf,
                   x = MAP, xend = MAP),
               lty = 2, col = "#d95f02") +
  geom_segment(data = int.map$grasslands,
               aes(y = 0.5, yend = 5.5,
                   x = MAP, xend = MAP),
               lty = 2, col = "#d95f02") +
  geom_segment(data = hdi.scenarios,
               aes(y = scenario, yend = scenario, x = CI_low,
                   xend = CI_high, size = hdi.type),
               col = 1) +
  geom_point(aes(x = MAP, y = scenario), size = 2.55, col = "white") +
  geom_point(aes(x = MAP, y = scenario), size = 2, col = 1) +
  scale_size_manual(values = c(core = 1.5, outer = .5), guide = F) +
  scale_alpha_continuous(range = c(.1, .9), guide = F) +
  scale_y_discrete(labels = labels.x) +
  theme(strip.background = element_rect(fill = "grey80", colour = "grey80"),
        panel.spacing = unit(2, "lines"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  xlab("Predicted value") +
  ylab("Scenario") +
  facet_wrap(~ res, scales = "free_x", labeller = label_parsed, nrow = 2)

empty <- ggplot() + theme_nothing()

plot_grid(p,
          plot_grid(empty, image_f, image_g, empty,
                    image_f, image_g, empty, ncol = 1,
                    rel_heights = c(1, 5, 5, 2.5, 5, 5, 1.4)),
          nrow = 1, rel_widths = c(10, 1))















