
# Convergence Checks for all real data examples ---------------------------
library(cmdstanr)
library(tidyverse)




convergence_all <- NULL


# Shorebirds --------------------------------------------------------------


sp_data <- "Shorebird"
model <- "gamye"
model_variant <- "spatial"
fit <- readRDS(paste0(paste("output/fit",sp_data,model_variant,model,
                            sep = "_"),".rds"))
conv_tmp <- fit$summary() %>% 
  mutate(species = sp_data,
         model = model,
         model_variant = model_variant)
convergence_all <- bind_rows(convergence_all,conv_tmp)

model_variant <- "hier"
fit <- readRDS(paste0(paste("output/fit",sp_data,model_variant,model,
                            sep = "_"),".rds"))
conv_tmp <- fit$summary()%>% 
  mutate(species = sp_data,
         model = model,
         model_variant = model_variant)
convergence_all <- bind_rows(convergence_all,conv_tmp)


# CBC --------------------------------------------------------------


sp_data <- "CBC"

model <- "gamye"
model_variant <- "spatial"
fit <- readRDS(paste0(paste("output/fit",sp_data,model_variant,model,
                            sep = "_"),".rds"))
conv_tmp <- fit$summary() %>% 
  mutate(species = sp_data,
         model = model,
         model_variant = model_variant)
convergence_all <- bind_rows(convergence_all,conv_tmp)

model_variant <- "hier"
fit <- readRDS(paste0(paste("output/fit",sp_data,model_variant,model,
                            sep = "_"),".rds"))
conv_tmp <- fit$summary()%>% 
  mutate(species = sp_data,
         model = model,
         model_variant = model_variant)
convergence_all <- bind_rows(convergence_all,conv_tmp)




model <- "first_diff"
model_variant <- "spatial"
fit <- readRDS(paste0(paste("output/fit",sp_data,model_variant,model,
                            sep = "_"),".rds"))
conv_tmp <- fit$summary() %>% 
  mutate(species = sp_data,
         model = model,
         model_variant = model_variant)
convergence_all <- bind_rows(convergence_all,conv_tmp)

model_variant <- "hier"
fit <- readRDS(paste0(paste("output/fit",sp_data,model_variant,model,
                            sep = "_"),".rds"))
conv_tmp <- fit$summary()%>% 
  mutate(species = sp_data,
         model = model,
         model_variant = model_variant)
convergence_all <- bind_rows(convergence_all,conv_tmp)







# BBS --------------------------------------------------------------


sp_data <- "BBS"

model <- "gamye"
model_variant <- "spatial"
conv_tmp <- readRDS(paste0("output/",paste(model,model_variant,
                            sep = "_"),"_param_summary.rds"))

conv_tmp <- conv_tmp %>% 
  mutate(species = sp_data,
         model = model,
         model_variant = model_variant)
convergence_all <- bind_rows(convergence_all,conv_tmp)

model_variant <- "hier"
conv_tmp <- readRDS(paste0("output/",paste(model,model_variant,
                                           sep = "_"),"_param_summary.rds"))

conv_tmp <- conv_tmp %>% 
  mutate(species = sp_data,
         model = model,
         model_variant = model_variant)
convergence_all <- bind_rows(convergence_all,conv_tmp)





model <- "first_diff"
model_variant <- "spatial"
conv_tmp <- readRDS(paste0("output/",paste(model,model_variant,
                                           sep = "_"),"_param_summary.rds"))

conv_tmp <- conv_tmp %>% 
  mutate(species = sp_data,
         model = model,
         model_variant = model_variant)
convergence_all <- bind_rows(convergence_all,conv_tmp)

model_variant <- "hier"
conv_tmp <- readRDS(paste0("output/",paste(model,model_variant,
                                           sep = "_"),"_param_summary.rds"))

conv_tmp <- conv_tmp %>% 
  mutate(species = sp_data,
         model = model,
         model_variant = model_variant)
convergence_all <- bind_rows(convergence_all,conv_tmp)


convergence_all <- convergence_all %>% 
  mutate(variable_type = stringr::str_extract(variable, "^\\w+"))



saveRDS(convergence_all,
        file = "output/convergence_all_real_data.rds")



ess_fail <- convergence_all %>% 
  filter(ess_bulk < 400)
rhat_fail <- convergence_all %>% 
  filter(rhat > 1.02)

convergence_summary <- convergence_all %>% 
  filter(is.finite(rhat),
         is.finite(ess_bulk)) %>% 
  group_by(species,model,model_variant,variable_type) %>% 
  summarise(max_rhat = max(rhat,na.rm = TRUE),
            min_ess_bulk = min(ess_bulk,na.rm = TRUE),
            .groups = "keep")
