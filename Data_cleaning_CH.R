library(tidyverse)
library(tsbox)
library(openxlsx)

# create covid index averaging by cantons weights in CH economy --------------------------

df_cantonal<-read.csv("inputData/GDP_per_canton.csv")
names(df_cantonal)<-c("canton","GDP")
df_cantonal<- df_cantonal %>% mutate(weight = GDP/sum(GDP))
sum(df_cantonal$weight) #sanity check
df_cantonal<- df_cantonal %>% select(-GDP) %>% 
  pivot_wider(names_from = canton, values_from = weight) #wide format

df_covid<-read.csv("inputData/raw_covid_CH.csv")
names(df_covid) <- gsub("ch\\.kof\\.stringency\\.(\\w{2})\\.stringency_plus",
                        "\\U\\1", names(df_covid), perl = TRUE) # clean canton names

df_covid <- df_covid %>% select(-"CH") %>%
  mutate(covid.ind = rowSums(across(AG:ZH) * as.numeric(df_cantonal[1, ]))) %>% 
  select(c(date,covid.ind)) %>% ts_ts() %>% ts_frequency("quarter") %>% ts_df()

#linearly interpolate from Q1 2022 (last value) to Q4 2024 
last_positive <- df_covid %>% filter(value > 0) %>% slice_tail(n = 1)
start_time <- last_positive$time
end_time <- as.Date("2024-10-01")
interp_dates <- seq(start_time, end_time, by = "quarter")

interp_values <- approx(
  x = c(start_time, end_time),
  y = c(last_positive$value, 0), # Interpolating from last positive value to 0
  xout = interp_dates
)

interp_df <- data.frame(time = interp_values$x, value = interp_values$y)
df_covid <- bind_rows(df_covid %>% filter(time < start_time), interp_df)

# create final dataset --------------------------------------------------------------

df_data<-readxl::read_xlsx("inputData/raw_CH.xlsx")
names(df_data)<-c("date","GDP","CPI_core","CPI_2020", "interest")
df_data$date<- df_data$date %>% as.Date()

# df_data$inflation<- select(df_data, c(date, CPI_2020)) %>% ts_ts() %>% ts_pcy() #prof. method (inaccurate?)

df_data<- df_data %>%
  mutate(
    gdp.log = log(GDP, base = exp(1)),
    inflation = (CPI_2020 - lag(CPI_2020)) / lag(CPI_2020),
    inflation = (1 + inflation)^4 - 1,  # Annualize inflation
    inflation.expectations = (lag(inflation, 1) + lag(inflation, 2) + lag(inflation, 3) + lag(inflation, 4)) / 4
  )



df_final <- df_data %>%  #include covid.index
  full_join(df_covid, by = c("date" = "time")) %>% 
  rename(covid.ind=value) %>% 
  mutate(covid.ind = ifelse(is.na(covid.ind), 0, covid.ind)) %>% 
  select(c("date","gdp.log","inflation","inflation.expectations","interest","covid.ind")) %>% 
  na.omit()


wb <- loadWorkbook("inputData/Input.xlsx")
writeData(wb, sheet = "CH input data", x = df_final)
saveWorkbook(wb, "inputData/Input.xlsx", overwrite = TRUE)

