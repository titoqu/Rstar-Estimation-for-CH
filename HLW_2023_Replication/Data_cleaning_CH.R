library(tidyverse)

df<-readxl::read_xlsx("inputData/Data_prof.xlsx")
names(df)<-c("time","GDP","CPI_core","CPI_2000", "interest")

# create covid index averaging by cantons weights in CH economy
df_cantonal<-read.csv("inputData/GDP_per_canton.csv")
names(df_cantonal)<-c("canton","GDP")
df_cantonal<- df_cantonal %>% mutate(weight = GDP/sum(GDP))
sum(df_cantonal$weight) #sanity check
df_cantonal<- df_cantonal %>% select(-GDP) %>% 
  pivot_wider(names_from = canton, values_from = weight) #wide format

df_covid<-read.csv("inputData/Data_prof_COVID.csv")
names(df_covid) <- gsub("ch\\.kof\\.stringency\\.(\\w{2})\\.stringency_plus",
                        "\\U\\1", names(df_covid), perl = TRUE) # clean canton names
df_covid<-df_covid %>% select(-"CH") # in theory indicates the federal index, but not precise

#average out using weights
df_covid <- df_covid %>%
  mutate(weighted_sum = rowSums(across(AG:ZH) * as.numeric(df_cantonal[1, ])))



