# run the different structural equation models for different corrections / scenarios
# needs data sets nw.metrics.wintr, nw.metrics.swnet, nw.metrics.basic from main analyses

# additional packages
library(nlme)
library(piecewiseSEM)

################################################################################.
# PART 1: SEMs not including land-use intensity ----- ##########################
################################################################################.

# For both ecosystems, run SEM and summarize pathways to explain relations
# between network structure and robustness

# set some labels, colors 

pathways.cols <- c(con = "#ffffbf", qmod = "#d7191c", nested = "#2b83ba",
                   "con qmod" = "#fdae61", "con nested" = "#abdda4",
                   direct = "#a6611a")

labels.res <- c(qmod = "Modularity",
                nested = "Nestedness",
                rob.A = "Robustness")

pathways.labels <- c("Direct", "Connectance",
                     expression(Connectance %->% Modularity),
                     "Modularity",
                     expression(Connectance %->% Nestedness),
                     "Nestedness")

# ------------------------------------------------------------------------------.
# define an auxiliary function -------------------------------------------------
# ------------------------------------------------------------------------------.

# from a data frame containing the SEM coefficients, produce a summary

f.pathwaysummary <- function(d.sem.coefs){
  
  # set up data frame with all possible pathways
  pathways <- data.frame(Nr = 1, predictor = "nwsize_sqrt", intermediate1 = NA, intermediate2 = NA,
                         stringsAsFactors = F) %>% 
    add_row(Nr = 2, predictor = "nwsize_sqrt", intermediate1 = "con", intermediate2 = NA) %>% 
    add_row(Nr = 3, predictor = "nwsize_sqrt", intermediate1 = "con", intermediate2 = "qmod") %>% 
    add_row(Nr = 4, predictor = "nwsize_sqrt", intermediate1 = "con", intermediate2 = "nested") %>% 
    add_row(Nr = 5, predictor = "nwsize_sqrt", intermediate1 = "qmod", intermediate2 = NA) %>% 
    add_row(Nr = 6, predictor = "nwsize_sqrt", intermediate1 = "nested", intermediate2 = NA) %>% 
    add_row(Nr = 7, predictor = "con", intermediate1 = NA, intermediate2 = NA) %>% 
    add_row(Nr = 8, predictor = "con", intermediate1 = "qmod", intermediate2 = NA) %>% 
    add_row(Nr = 9, predictor = "con", intermediate1 = "nested", intermediate2 = NA) %>% 
    add_row(Nr = 10, predictor = "qmod", intermediate1 = NA, intermediate2 = NA) %>% 
    add_row(Nr = 11, predictor = "nested", intermediate1 = NA, intermediate2 = NA)
  
  # calculate pathways estiamtes
  for (i in 1:nrow(pathways)){
    
    predictor_i <- pathways[i, "predictor"]
    intermediates <- pathways[i, c("intermediate1", "intermediate2")]
    intermediates <- intermediates[!is.na(intermediates)]
    
    out <- 1
    if (length(intermediates) > 0){ # if i is an indirect pathway
      for (j in 1:length(intermediates)){
        out <- out * d.sem.coefs$estimate[d.sem.coefs$predictor == predictor_i & 
                                            d.sem.coefs$response == intermediates[j]]
        predictor_i <- intermediates[j]
      }
      
    } 
    out <- out * d.sem.coefs$estimate[d.sem.coefs$predictor == predictor_i & 
                                        d.sem.coefs$response == "rob.A"]
    pathways$estimate[i] <- out
  }
  
  # producte output
  pathways <- pathways %>% 
    mutate(intermediates = ifelse(is.na(intermediate1), "direct",
                                  ifelse(is.na(intermediate2),
                                         intermediate1,
                                         paste(intermediate1, intermediate2))),
           intermediates = factor(intermediates, 
                                  levels = c("direct", "con", "con qmod", "qmod", 
                                             "con nested", "nested")),
           predictor = factor(predictor, levels = c("nwsize_sqrt", "con", "qmod", "nested"),
                              labels = c("Network size", "Connectance", "Modularity", "Nestedness")))
  
  pathways
}

# ------------------------------------------------------------------------------.
# run SEMs and plot the summary ------------------------------------------------
# ------------------------------------------------------------------------------.

# summary data frames
fit.sem <- data.frame()
smry.tab <- data.frame()

plots.pathways <- list()

# loop through different sensitivity analyses
for (datatype.def in c("raw", "no.poly", "all.plants", "pres.abs")){
  
  # scale combined data
  data_comb <- nw.metrics.wintr[[datatype.def]] %>% 
    bind_rows(nw.metrics.swnet[[datatype.def]]) %>% 
    left_join(nw.metrics.basic[[datatype.def]], by = "PlotID") %>% 
    mutate_at(vars(nwsize_sqrt, con, qmod, nested), # NOT ROBUSTNESS
              ~ as.vector(scale(.)))
  
  # forests --------------------------------------------------------------------.

  data <- data_comb %>% 
    filter(Landuse == "F")
  
  modlist <- list(
    con = lme(con ~ nwsize_sqrt,
              random = ~ 1 | Exploratory,
              data),
    qmod = lme(qmod ~ con + nwsize_sqrt,
               random = ~ 1 | Exploratory,
               data),
    nested = lme(nested ~ con + nwsize_sqrt,
                 random = ~ 1 | Exploratory,
                 data),
    rob.A = lme(rob.A ~ con + nwsize_sqrt +nested + qmod,
                random = ~ 1 | Exploratory,
                data)
  )
  

  R2 <- rsquared(modlist) 
  R2 <- R2 %>% 
    rownames_to_column("response")
  
  fit.f <- sem.fit(modlist, data)

  sem.coefs.f <- sem.coefs(modlist, data = data, standardize = "none")  %>% 
    data.frame()
  
  smry.tab.f <- sem.coefs.f %>% 
    mutate(estimate = formatC(estimate, digits = 2),
           estimate = paste0(estimate, Var.6),
           predictor = factor(predictor, levels = c("nwsize_sqrt", "con", "qmod", "nested"))) %>%
    select(response, predictor, estimate) %>% 
    spread(predictor, estimate) %>% 
    left_join(R2 %>% 
                select(response, Marginal, Conditional) %>% 
                mutate_at(vars(Marginal, Conditional), 
                          ~formatC(., digits = 2, format = "f")), 
              by = "response")
  
  # grasslands -----------------------------------------------------------------.

  data <- data_comb %>% 
    filter(Landuse == "G")
  
  modlist <- list(
    con = lme(con ~ nwsize_sqrt,
              random = ~ 1 | Exploratory,
              data),
    qmod = lme(qmod ~ con + nwsize_sqrt,
               random = ~ 1 | Exploratory,
               data),
    nested = lme(nested ~ con + nwsize_sqrt,
                 random = ~ 1 | Exploratory,
                 data),
    rob.A = lme(rob.A ~ con + nwsize_sqrt +nested + qmod,
                random = ~ 1 | Exploratory,
                data)
  )
  
  R2 <- rsquared(modlist) 
  R2 <- R2 %>% 
    rownames_to_column("response")
  
  fit.g <- sem.fit(modlist, data)
  
  
  sem.coefs.g <- sem.coefs(modlist, standardize = "none", data = data)  %>% 
    data.frame()
  
  smry.tab.g <- sem.coefs.g %>% 
    mutate(estimate = formatC(estimate, digits = 2),
           estimate = paste0(estimate, Var.6),
           predictor = factor(predictor, levels = c("nwsize_sqrt", "con", "qmod", "nested"))) %>%
    select(response, predictor, estimate) %>% 
    spread(predictor, estimate) %>% 
    left_join(R2 %>% 
                select(response, Marginal, Conditional) %>% 
                mutate_at(vars(Marginal, Conditional), 
                          ~formatC(., digits = 2, format = "f")), 
              by = "response")
  
  
  # plotting -------------------------------------------------------------------.
  # save to a list
  
  pathways.f <- f.pathwaysummary(sem.coefs.f)
  
  pathways.g <- f.pathwaysummary(sem.coefs.g)
  
  pathways <- pathways.g %>% 
    mutate(system = "Grasslands") %>% 
    bind_rows(pathways.f %>% 
                mutate(system = "Forests")) %>% 
    mutate(system = factor(system, levels = c("Grasslands", "Forests")))
  

  pathways_sum <- pathways %>% 
    group_by(predictor, system) %>% 
    summarise(estimate = sum(estimate), .groups = "drop")
  

  
  plots.pathways[[datatype.def]] <-
    pathways %>% 
    mutate(predictor = factor(predictor, levels = rev(levels(predictor)))) %>% 
    ggplot(aes(x = predictor, y = estimate)) +
    geom_hline(yintercept = 0) +
    geom_col(aes(fill = intermediates), col = 1) +
    geom_segment(data = pathways_sum, aes(xend = predictor, yend = 0),
                 arrow = arrow(end = "first", angle = 90), size = 1.5) +
    scale_fill_manual(values = pathways.cols, labels = pathways.labels,
                      name = "Pathway") +
    xlab("Predictor") +
    ylab("Estimate") +
    facet_grid(~ system) +
    coord_flip() +
    theme(legend.text.align = 0)
  

  # summary table --------------------------------------------------------------

  smry.tab <- smry.tab %>% 
    bind_rows(smry.tab.f %>% 
                mutate(datatype = datatype.def,
                       system = "Forests")) %>% 
    bind_rows(smry.tab.g %>% 
                mutate(datatype = datatype.def,
                       system = "Grasslands"))
  
  fit.sem <- fit.sem %>% 
    bind_rows(fit.f$Fisher.C %>% 
                mutate(datatype = datatype.def,
                       system = "Forests")) %>% 
    bind_rows(fit.g$Fisher.C %>% 
                mutate(datatype = datatype.def,
                       system = "Grasslands"))
  
}


# summary output ---------------------------------------------------------------.

smry.tab %>% 
  mutate(system = factor(system, levels = c("Grasslands", "Forests"))) %>% 
  arrange(system)

fit.sem %>% 
  mutate(p_val_code = p_val_code(p.value),
         fisher.c = paste0(fisher.c, p_val_code))


# ------------------------------------------------------------------------------.
# Plot for main manuscript -----------------------------------------------------
# ------------------------------------------------------------------------------.

plots.pathways$raw +
  theme(strip.text = element_blank(),
        plot.margin = margin(2, 0, 0, 0, "cm"),
        axis.text = element_text(size = 12.5),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12.5))

# ------------------------------------------------------------------------------.
# Plot for supplement ----------------------------------------------------------
# ------------------------------------------------------------------------------.

plot_grid(
  plots.pathways$no.poly + theme(strip.text = element_blank(),
                                 legend.position = "none",
                                 axis.title.x = element_blank(),
                                 plot.margin = margin(2, 0, 0, 0, "cm"),
                                 axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)),
  plots.pathways$all.plants + theme(strip.text = element_blank(),
                                    legend.position = "none",
                                    axis.text.y = element_blank(),
                                    axis.title.y = element_blank(),
                                    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)),
  plots.pathways$pres.abs + theme(strip.text = element_blank(),
                                  axis.text.y = element_blank(),
                                  axis.title.y = element_blank(),
                                  axis.title.x = element_blank(),
                                  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)),
  nrow = 1, rel_widths = c(1, .7, 1.4), align = "h")


################################################################################.
# PART 2: SEMs including land-use intensity --------- ##########################
################################################################################.

# run SEMs for different sensitivity scenarios and produce output

for (datatype.def in c("raw", "no.poly", "all.plants", "pres.abs")){
  
  # ----------------------------------------------------------------------------.
  # forests --------------------------------------------------------------------
  # ----------------------------------------------------------------------------.

  data <- nw.metrics.wintr[[datatype.def]] %>% 
    left_join(nw.metrics.basic[[datatype.def]], by = "PlotID") %>% 
    mutate_at(vars(nwsize_sqrt, con, qmod, nested, rob.A),
              ~ as.vector(scale(.)))
  
  modlist <- list(
    mod_nwsize = lme(nwsize_sqrt ~ Inonat_z + Iharv_z + Idwcut_z,
               random = ~ 1 | Exploratory,
               data),
    mod_con = lme(con ~ Inonat_z + Iharv_z + Idwcut_z + nwsize_sqrt,
               random = ~ 1 | Exploratory,
               data),
    mod_qmod = lme(qmod ~ Inonat_z + Iharv_z + Idwcut_z + con + nwsize_sqrt,
               random = ~ 1 | Exploratory,
               data),
    mod_nested = lme(nested ~ Inonat_z + Iharv_z + Idwcut_z + con + nwsize_sqrt,
               random = ~ 1 | Exploratory,
               data),
    mod_robA = lme(rob.A ~ Inonat_z + Iharv_z + Idwcut_z + con + nwsize_sqrt +nested + qmod,
               random = ~ 1 | Exploratory,
               data)
  )

  (R2 <- rsquared(modlist)) 
  R2 <- R2 %>% 
    rownames_to_column("Mod")
  
  (fit <- sem.fit(modlist, data, corr.errors = c("Inonat_z ~~ Iharv_z", 
                                                 "Inonat_z ~~ Idwcut_z", 
                                                 "Idwcut_z ~~ Iharv_z")))
  
  sem.plot(modlist,
           coef.table = sem.coefs(modlist, data = data, 
                                  corr.errors = c("Inonat_z ~~ Iharv_z", 
                                                  "Inonat_z ~~ Idwcut_z", 
                                                  "Idwcut_z ~~ Iharv_z")))
  
  
  (sem.coefs.type <- sem.coefs(modlist, data = data, standardize = "none",
                               corr.errors = c("Inonat_z ~~ Iharv_z", 
                                               "Inonat_z ~~ Idwcut_z", 
                                               "Idwcut_z ~~ Iharv_z"))  %>% 
      data.frame() %>% 
      mutate(estimate10 = 10 * estimate))
  
  
   # p values for correlation are not correct. Do correctly:
  sem.coefs.type$p.value[sem.coefs.type$response == "~~ Inonat_z" &
                           sem.coefs.type$predictor == "~~ Idwcut_z"]  <- 
    cor.test(data$Inonat_z, data$Idwcut_z)$p.value
  
  sem.coefs.type$p.value[sem.coefs.type$response == "~~ Inonat_z" &
                           sem.coefs.type$predictor == "~~ Iharv_z"]  <- 
    cor.test(data$Inonat_z, data$Iharv_z)$p.value
  
  sem.coefs.type$p.value[sem.coefs.type$response == "~~ Idwcut_z" &
                           sem.coefs.type$predictor == "~~ Iharv_z"]  <- 
    cor.test(data$Idwcut_z, data$Iharv_z)$p.value
  
  sem.coefs.type$Var.6 <- p_val_code(sem.coefs.type$p.value)
  
  
  (sem.coefs.type <- sem.coefs.type %>% 
      mutate(predictor = factor(predictor, levels = c("Inonat_z", "Iharv_z", "Idwcut_z", "nwsize_sqrt", "con", "qmod", "nested"))) %>% 
      arrange(response, predictor))
  
  # effect pathways Inonat -----------------------------------------------------.
  
  # indirect basic
  (i.b.nonat <- sem.coefs.type$estimate[1] * sem.coefs.type$estimate[21] +
      sem.coefs.type$estimate[1] * sem.coefs.type$estimate[7] * sem.coefs.type$estimate[22] +
      sem.coefs.type$estimate[4] * sem.coefs.type$estimate[22])
  
  # indirect modularity
  (i.mod.nonat <- sem.coefs.type$estimate[1] * sem.coefs.type$estimate[11] * sem.coefs.type$estimate[23] +
      sem.coefs.type$estimate[1] * sem.coefs.type$estimate[7] *  sem.coefs.type$estimate[12] * sem.coefs.type$estimate[23] +
      sem.coefs.type$estimate[4] * sem.coefs.type$estimate[12] * sem.coefs.type$estimate[23] +
      sem.coefs.type$estimate[8] * sem.coefs.type$estimate[23])
  
  # indirect nestedness
  (i.nested.nonat <- sem.coefs.type$estimate[1] * sem.coefs.type$estimate[16] * sem.coefs.type$estimate[24] +
      sem.coefs.type$estimate[1] * sem.coefs.type$estimate[7] *  sem.coefs.type$estimate[17] * sem.coefs.type$estimate[24] +
      sem.coefs.type$estimate[4] * sem.coefs.type$estimate[17] * sem.coefs.type$estimate[24] +
      sem.coefs.type$estimate[13] * sem.coefs.type$estimate[24])
  
  # direct
  (d.nonat <- sem.coefs.type$estimate[18])
  
  (sum.nonat <- i.b.nonat + i.mod.nonat + i.nested.nonat + d.nonat)
  
  # effect pathways Iharv ------------------------------------------------------.
  
  # indirect basic
  (i.b.harv <- sem.coefs.type$estimate[2] * sem.coefs.type$estimate[21] +
      sem.coefs.type$estimate[2] * sem.coefs.type$estimate[7] * sem.coefs.type$estimate[22] +
      sem.coefs.type$estimate[5] * sem.coefs.type$estimate[22])
  
  # indirect modularity
  (i.mod.harv <- sem.coefs.type$estimate[2] * sem.coefs.type$estimate[11] * sem.coefs.type$estimate[23] +
      sem.coefs.type$estimate[2] * sem.coefs.type$estimate[7] *  sem.coefs.type$estimate[12] * sem.coefs.type$estimate[23] +
      sem.coefs.type$estimate[5] * sem.coefs.type$estimate[12] * sem.coefs.type$estimate[23] +
      sem.coefs.type$estimate[9] * sem.coefs.type$estimate[23])
  
  # indirect nestedness
  (i.nested.harv <- sem.coefs.type$estimate[2] * sem.coefs.type$estimate[16] * sem.coefs.type$estimate[24] +
      sem.coefs.type$estimate[2] * sem.coefs.type$estimate[7] *  sem.coefs.type$estimate[17] * sem.coefs.type$estimate[24] +
      sem.coefs.type$estimate[5] * sem.coefs.type$estimate[17] * sem.coefs.type$estimate[24] +
      sem.coefs.type$estimate[14] * sem.coefs.type$estimate[24])
  
  # direct
  (d.harv <- sem.coefs.type$estimate[19])
  
  (sum.harv <- i.b.harv + i.mod.harv + i.nested.harv + d.harv)
  
  # effect pathways Idwcut -----------------------------------------------------.
  
  # indirect basic
  (i.b.dwcut <- sem.coefs.type$estimate[3] * sem.coefs.type$estimate[21] +
      sem.coefs.type$estimate[3] * sem.coefs.type$estimate[7] * sem.coefs.type$estimate[22] +
      sem.coefs.type$estimate[6] * sem.coefs.type$estimate[22])
  
  # indirect modularity
  (i.mod.dwcut <- sem.coefs.type$estimate[3] * sem.coefs.type$estimate[11] * sem.coefs.type$estimate[23] +
      sem.coefs.type$estimate[3] * sem.coefs.type$estimate[7] *  sem.coefs.type$estimate[12] * sem.coefs.type$estimate[23] +
      sem.coefs.type$estimate[6] * sem.coefs.type$estimate[12] * sem.coefs.type$estimate[23] +
      sem.coefs.type$estimate[10] * sem.coefs.type$estimate[23])
  
  # indirect nestedness
  (i.nested.dwcut <- sem.coefs.type$estimate[3] * sem.coefs.type$estimate[16] * sem.coefs.type$estimate[24] +
      sem.coefs.type$estimate[3] * sem.coefs.type$estimate[7] *  sem.coefs.type$estimate[17] * sem.coefs.type$estimate[24] +
      sem.coefs.type$estimate[6] * sem.coefs.type$estimate[17] * sem.coefs.type$estimate[24] +
      sem.coefs.type$estimate[15] * sem.coefs.type$estimate[24])
  
  # direct
  (d.dwcut <- sem.coefs.type$estimate[20])
  
  (sum.dwcut <- i.b.dwcut + i.mod.dwcut + i.nested.dwcut + d.dwcut)
  
  
  smry <- data.frame(driver = "Inonat", i.basic = i.b.nonat, i.mod = i.mod.nonat, 
                     i.nested = i.nested.nonat, direct = d.nonat, total = sum.nonat,
                     stringsAsFactors = F) %>% 
    bind_rows(data.frame(driver = "Iharv", i.basic = i.b.harv, i.mod = i.mod.harv, 
                         i.nested = i.nested.harv, direct = d.harv, total = sum.harv,
                         stringsAsFactors = F)) %>% 
    bind_rows(data.frame(driver = "Idwcut", i.basic = i.b.dwcut, i.mod = i.mod.dwcut, 
                         i.nested = i.nested.dwcut, direct = d.dwcut, total = sum.dwcut,
                         stringsAsFactors = F)) %>% 
    mutate_at(vars(-driver), ~round(., 2))
  
  #-----------------------------------------------------------------------------.
  # grasslands -----------------------------------------------------------------
  #-----------------------------------------------------------------------------.
  
  data <- nw.metrics.swnet[[datatype.def]] %>% 
    left_join(nw.metrics.basic[[datatype.def]], by = "PlotID") %>% 
    mutate_at(vars(nwsize_sqrt, con, qmod, nested, rob.A),
              ~ as.vector(scale(.)))
  
  modlist <- list(
    mod_nwsize = lme(nwsize_sqrt ~ M_z + F_z + G_z,
               random = ~ 1 | Exploratory,
               data),
    mod_con = lme(con ~ M_z + F_z + G_z + nwsize_sqrt,
               random = ~ 1 | Exploratory,
               data),
    mod_qmod = lme(qmod ~ M_z + F_z + G_z + con + nwsize_sqrt,
               random = ~ 1 | Exploratory,
               data),
    mod_nested = lme(nested ~ M_z + F_z + G_z + con + nwsize_sqrt,
               random = ~ 1 | Exploratory,
               data),
    mod_robA = lme(rob.A ~ M_z + F_z + G_z + con + nwsize_sqrt +nested + qmod,
               random = ~ 1 | Exploratory,
               data)
  )

  (R2 <- rsquared(modlist)) 
  R2 <- R2 %>% 
    rownames_to_column("Mod")
  
  (fit <- sem.fit(modlist, data, corr.errors = c("M_z ~~ F_z", 
                                                 "M_z ~~ G_z", 
                                                 "F_z ~~ G_z")))
  sem.plot(modlist,
           coef.table = sem.coefs(modlist, data = data, 
                                  corr.errors = c("M_z ~~ F_z", 
                                                  "M_z ~~ G_z", 
                                                  "F_z ~~ G_z")))
  
  (sem.coefs.type <- sem.coefs(modlist, standardize = "none", data = data,
                               corr.errors = c("M_z ~~ F_z", 
                                               "M_z ~~ G_z", 
                                               "F_z ~~ G_z"))  %>% 
      data.frame() %>% 
      mutate(estimate10 = 10 * estimate))
  
  
  # p values for correlation are not correct. Do correctly:
  sem.coefs.type$p.value[sem.coefs.type$response == "~~ M_z" &
                           sem.coefs.type$predictor == "~~ F_z"]  <- 
      cor.test(data$M_std, data$F_std)$p.value
  
  sem.coefs.type$p.value[sem.coefs.type$response == "~~ M_z" &
                           sem.coefs.type$predictor == "~~ G_z"]  <- 
    cor.test(data$M_std, data$G_std)$p.value
  
  sem.coefs.type$p.value[sem.coefs.type$response == "~~ F_z" &
                           sem.coefs.type$predictor == "~~ G_z"]  <- 
    cor.test(data$F_std, data$G_std)$p.value
  
  sem.coefs.type$Var.6 <- p_val_code(sem.coefs.type$p.value)

      
      
  (sem.coefs.type <- sem.coefs.type %>% 
      mutate(predictor = factor(predictor, levels = c("M_z", "G_z", "F_z", "nwsize_sqrt", "con", "qmod", "nested"))) %>% 
      arrange(response, predictor))
  
  # effect pathways mowing -----------------------------------------------------.
  
  # indirect basic
  (i.b.mow <- sem.coefs.type$estimate[1] * sem.coefs.type$estimate[21] +
      sem.coefs.type$estimate[1] * sem.coefs.type$estimate[7] * sem.coefs.type$estimate[22] +
      sem.coefs.type$estimate[4] * sem.coefs.type$estimate[22])
  
  # indirect modularity
  (i.mod.mow <- sem.coefs.type$estimate[1] * sem.coefs.type$estimate[11] * sem.coefs.type$estimate[23] +
      sem.coefs.type$estimate[1] * sem.coefs.type$estimate[7] *  sem.coefs.type$estimate[12] * sem.coefs.type$estimate[23] +
      sem.coefs.type$estimate[4] * sem.coefs.type$estimate[12] * sem.coefs.type$estimate[23] +
      sem.coefs.type$estimate[8] * sem.coefs.type$estimate[23])
  
  # indirect nestedness
  (i.nested.mow <- sem.coefs.type$estimate[1] * sem.coefs.type$estimate[16] * sem.coefs.type$estimate[24] +
      sem.coefs.type$estimate[1] * sem.coefs.type$estimate[7] *  sem.coefs.type$estimate[17] * sem.coefs.type$estimate[24] +
      sem.coefs.type$estimate[4] * sem.coefs.type$estimate[17] * sem.coefs.type$estimate[24] +
      sem.coefs.type$estimate[13] * sem.coefs.type$estimate[24])
  
  # direct
  (d.mow <- sem.coefs.type$estimate[18])
  
  (sum.mow <- i.b.mow + i.mod.mow + i.nested.mow + d.mow)
  
  # effect pathways grazing ----------------------------------------------------.
  
  # indirect basic
  (i.b.graz <- sem.coefs.type$estimate[2] * sem.coefs.type$estimate[21] +
      sem.coefs.type$estimate[2] * sem.coefs.type$estimate[7] * sem.coefs.type$estimate[22] +
      sem.coefs.type$estimate[5] * sem.coefs.type$estimate[22])
  
  # indirect modularity
  (i.mod.graz <- sem.coefs.type$estimate[2] * sem.coefs.type$estimate[11] * sem.coefs.type$estimate[23] +
      sem.coefs.type$estimate[2] * sem.coefs.type$estimate[7] *  sem.coefs.type$estimate[12] * sem.coefs.type$estimate[23] +
      sem.coefs.type$estimate[5] * sem.coefs.type$estimate[12] * sem.coefs.type$estimate[23] +
      sem.coefs.type$estimate[9] * sem.coefs.type$estimate[23])
  
  # indirect nestedness
  (i.nested.graz <- sem.coefs.type$estimate[2] * sem.coefs.type$estimate[16] * sem.coefs.type$estimate[24] +
      sem.coefs.type$estimate[2] * sem.coefs.type$estimate[7] *  sem.coefs.type$estimate[17] * sem.coefs.type$estimate[24] +
      sem.coefs.type$estimate[5] * sem.coefs.type$estimate[17] * sem.coefs.type$estimate[24] +
      sem.coefs.type$estimate[14] * sem.coefs.type$estimate[24])
  
  # direct
  (d.graz <- sem.coefs.type$estimate[19])
  
  (sum.graz <- i.b.graz + i.mod.graz + i.nested.graz + d.graz)
  
  # effect pathways fertilisation ----------------------------------------------.
  
  # indirect basic
  (i.b.fert <- sem.coefs.type$estimate[3] * sem.coefs.type$estimate[21] +
      sem.coefs.type$estimate[3] * sem.coefs.type$estimate[7] * sem.coefs.type$estimate[22] +
      sem.coefs.type$estimate[6] * sem.coefs.type$estimate[22])
  
  # indirect modularity
  (i.mod.fert <- sem.coefs.type$estimate[3] * sem.coefs.type$estimate[11] * sem.coefs.type$estimate[23] +
      sem.coefs.type$estimate[3] * sem.coefs.type$estimate[7] *  sem.coefs.type$estimate[12] * sem.coefs.type$estimate[23] +
      sem.coefs.type$estimate[6] * sem.coefs.type$estimate[12] * sem.coefs.type$estimate[23] +
      sem.coefs.type$estimate[10] * sem.coefs.type$estimate[23])
  
  # indirect nestedness
  (i.nested.fert <- sem.coefs.type$estimate[3] * sem.coefs.type$estimate[16] * sem.coefs.type$estimate[24] +
      sem.coefs.type$estimate[3] * sem.coefs.type$estimate[7] *  sem.coefs.type$estimate[17] * sem.coefs.type$estimate[24] +
      sem.coefs.type$estimate[6] * sem.coefs.type$estimate[17] * sem.coefs.type$estimate[24] +
      sem.coefs.type$estimate[15] * sem.coefs.type$estimate[24])
  
  # direct
  (d.fert <- sem.coefs.type$estimate[20])
  
  (sum.fert <- i.b.fert + i.mod.fert + i.nested.fert + d.fert)
  
  
  smry <- data.frame(driver = "Mowing", i.basic = i.b.mow, i.mod = i.mod.mow, 
                     i.nested = i.nested.mow, direct = d.mow, total = sum.mow,
                     stringsAsFactors = F) %>% 
    bind_rows(data.frame(driver = "Grazing", i.basic = i.b.graz, i.mod = i.mod.graz, 
                         i.nested = i.nested.graz, direct = d.graz, total = sum.graz,
                         stringsAsFactors = F)) %>% 
    bind_rows(data.frame(driver = "Fertilisation", i.basic = i.b.fert, i.mod = i.mod.fert, 
                         i.nested = i.nested.fert, direct = d.fert, total = sum.fert,
                         stringsAsFactors = F)) %>% 
    mutate_at(vars(-driver), ~round(., 2))
  
  # ----------------------------------------------------------------------------.
  # export ---------------------------------------------------------------------
  # ----------------------------------------------------------------------------.
  
  data.frame(
    pathways = c("indirect basic", "indirect qmod", "indirect nested", "direct", "total"),
    NoNat = c(i.b.nonat, i.mod.nonat, i.nested.nonat, d.nonat, sum.nonat),
    Harv = c(i.b.harv, i.mod.harv, i.nested.harv, d.harv, sum.harv),
    DWcut = c(i.b.dwcut, i.mod.dwcut, i.nested.dwcut, d.dwcut, sum.dwcut),
    Mowing = c(i.b.mow, i.mod.mow, i.nested.mow, d.mow, sum.mow),
    Grazing = c(i.b.graz, i.mod.graz, i.nested.graz, d.graz, sum.graz),
    Fertilisation = c(i.b.fert, i.mod.fert, i.nested.fert, d.fert, sum.fert)
  ) %>% column_to_rownames("pathways") %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column("Driver") %>% 
    write.table("targetfile")
  
}

