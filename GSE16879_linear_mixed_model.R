# load libraries
library(stringr)
library(lmerTest)
library(tidyverse)
library(purrr)

# load dataset
data <- read.csv("AntiTNF_expr_pheno.csv", check.names = FALSE)


# filter out controls
data %>%
  janitor::clean_names() %>% 
  filter(str_detect(data$ID, "UC|CD")) ->
data_filter_nocontrols


# Create new column for time point
data_filter_nocontrols %>%
  separate(data$ID, c("ID", "Relationship to IFX"), sep="_") ->
data_filter_nocontrols 


## Make Predictors Factor Variables
data_filter_nocontrols %>%
  mutate(Response = as.factor(Response)) %>%
  mutate(Time = as.factor(Time))  ->
data_filter_nocontrols


# Table of Responders by Time point (CDcR2 has no after IFX biopsy)
table(data_filter_nocontrols$Response, data_filter_nocontrols$Time)


##Model
model<-lmer(x1007_s_at ~ (1|ID) + Disease + Tissue + Response + Time + Response:Time, data = data_filter_nocontrols)
summary(model)


## Run model on each transcript (affx_trpn_x_m_at when ready to run whole dataset)
data_filter_nocontrols %>%
  dplyr::select(x1007_s_at:x1552379_at) %>%
  map(~lmer(.x ~ (1|ID) + Response + Time + Response:Time, data=data_filter_nocontrols)) %>%
  map_dfr(~ broom::tidy(.), .id = 'ID') ->
models_dtf


## Add FDR (False Discovery Rate) column
models_dtf %>%
  mutate(p.fdr = p.adjust(p.value, method="fdr")) ->
models_dtf


## Create separate data frame for interaction terms only
models_dtf %>%
  select(ID, term, estimate, p.value, p.fdr) %>%
  filter(term=='ResponseYes:TimeBefore first infliximab treatment') %>%
  select(-term) ->
interactions_only_dtf


## Remove "x" prefix assigned to numeric transcripts
sub("x", "", interactions_only_dtf$ID) ->
  interactions_only_dtf$ID

## Merge with transcript information
transcript_info <- read.csv("Transcript Information.csv")

transcript_info %>%
  select('ID', 'Gene.Title', 'Gene.Symbol')->
transcript_info

left_join(interactions_only_dtf, transcript_info, by= "ID") ->
Final_Results_Interactions_Transcripts


## Arrange by lowest p-value
Final_Results_Interactions_Transcripts %>%
  arrange(p.fdr) ->
Final_Results_Interactions_Transcripts

  
# Write output to file
write.csv(Final_Results_Interactions_Transcripts, file="Final_Results_Time_Response_Interactions_Transcripts.csv")
