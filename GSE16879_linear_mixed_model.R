# load libraries
library(stringr)
library(lmerTest)
library(tidyverse)
library(purrr)
library(easystats)

# load dataset
data <- read.csv("AntiTNF_expr_pheno.csv")

# filter out controls
data %>%
  filter(str_detect(title, "UC|CD")) ->
data_filter_nocontrols


# Create new column for time point
data_filter_nocontrols %>%
  separate(title, c("ID", "Relationship to IFX"), sep="_") ->
data_filter_nocontrols 


# Rename Response and Time Point Columns
data_filter_nocontrols %>%
  rename(Response=response.to.infliximab.ch1) %>%
  rename(Time=before.or.after.first.infliximab.treatment.ch1) ->
data_filter_nocontrols


# Table of Responders by Time point (CDcR2 has no after IFX biopsy)
table(data_filter_nocontrols$Response, data_filter_nocontrols$Time)


##Model
model<-lmer(X1007_s_at ~ (1|ID) + Response + Time + Response:Time, data = data_filter_nocontrols)
anova(model)

## Create vector of transcript only information for Model 
# data_filter_nocontrols %>%
#   select(X1007_s_at:AFFX.TrpnX.M_at) %>%
#   names() ->
# transcripts


## Run model on each transcript
data_filter_nocontrols %>%
  dplyr::select(X1007_s_at:X121_at) %>%
  map(~lmer(.x ~ (1|ID) + Response + Time + Response:Time, data=data_filter_nocontrols)) ->
models_list


# summary
models_list %>% 
  map(summary, .id = ".x") %>% 
  map_dfr(insight::get_parameters, .id = "transcript") %>% 
  janitor::clean_names() %>% 
  mutate(parameter = as.numeric(parameter)) %>% 
  filter(parameter <5) %>% 
  select(transcript: estimate_1, estimate_5) %>% 
  rename(estimate = estimate_1) %>% 
  rename(p_value = estimate_5) %>% 
  mutate(parameter = case_when(
    parameter == 1 ~ "intercept",
    parameter == 2 ~ "response_yes",
    parameter == 3 ~ "time_beforeIFX",
    parameter == 4 ~ "time*response"
  ))



## Other loop example
#Assign outcome variables (i.e. transcripts)
# out_start=9
# out_end= 54683
# out_nvar=out_end-out_start+1
# 
# out_variable=rep(NA, out_nvar)
# out_beta=rep(NA, out_nvar)
# out_se = rep(NA, out_nvar)
# out_pvalue=rep(NA, out_nvar)

# Loop to run model on every transcript
for (i in out_start:out_end){
  outcome = colnames(data_filter_nocontrols)[i]
    model <- lmer(get(outcome) ~ (1|ID) + Response*Time, data = data_filter_nocontrols)}
    
# Create dataframe with results
outcome = data.frame(out_variable, out_beta, out_se, out_pvalue)