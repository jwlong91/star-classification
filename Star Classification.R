# load libraries and set options
library(tidyverse)
library(GGally)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# load data
df <- read_csv('Stars.csv')

# check for na's
any(is.na(df))

# rename fields, remove zero indexing on star type, and create log temp
df <- df %>%
  rename(temp = `Temperature (K)`, lum = `Luminosity(L/Lo)`,
         rad = `Radius(R/Ro)`, mag = `Absolute magnitude(Mv)`,
         type = `Star type`, color = `Star color`,
         class = `Spectral Class`) %>%
  mutate(type = type+1,
         ln_temp = log(temp),
         star_type_title = as.factor(case_when(
           type == 1 ~ '1. Brown Dwarf',
           type == 2 ~ '2. Red Dwarf',
           type == 3 ~ '3. White Dwarf',
           type == 4 ~ '4. Main Sequence',
           type == 5 ~ '5. Supergiant',
           type == 6 ~ '6. Hypergiant'
         )))

# correlations between variables
ggcorr(df)

# investigate relationship between log temp and magnitude
df %>%
  ggplot(aes(x = log(temp), y = mag, color = star_type_title)) +
  geom_point() + 
  labs(title = 'Star Types by Temperature and Magnitude',
       x = 'Log Temperature (K)',
       y = 'Absolute Magnitude (Mv)',
       color = 'Star Type')

# shuffle data frame
rows <- sample(nrow(df))
df <- df[rows,]

# split to train and test data frames
dfTrain <- df[1:floor(nrow(df)*0.8),]
dfTest <- df[(floor(nrow(df)*0.8)+1):nrow(df),]

# set up data for stan
# beta1 = log temperature
# beta2 = magnitude
stars_dat <- list(N = nrow(dfTrain),
                  K = 6,
                  D = 2,
                  x = dfTrain %>% select(ln_temp, mag) %>% as.matrix(ncol = 2),
                  y = dfTrain$type,
                  N_new = nrow(dfTest),
                  x_new = dfTest %>% select(ln_temp, mag) %>% as.matrix(ncol = 2))

# set up stan model
m1 <- 
  "data {
int<lower=2> K; // number of categories
int<lower=1> N; // number of observations
int<lower=1> D; // number of predictors
int<lower=1,upper=K> y[N];
matrix[N,D] x;

int<lower=0> N_new;
matrix[N_new,D] x_new;
}

parameters {
matrix[D,K] beta;  // variable effects

}

model {
matrix[N,K] x_beta = x * beta;

to_vector(beta) ~ normal(0,5);

for (n in 1:N)
  y[n] ~ categorical_logit(x_beta[n]');
}

generated quantities {
matrix[N_new,K] y_new;

matrix[N_new,K] x_beta_new_exp = exp(x_new * beta);

vector[N_new] sum_eta;

for (n in 1:N_new)
  sum_eta[n] = sum(x_beta_new_exp[n]);

for (n in 1:N_new)
  for (k in 1:K)
    y_new[n,k] = x_beta_new_exp[n,k]/sum_eta[n];

}

"

# fit model
fit1 <- stan(model_code = m1, data = stars_dat, iter = 1000, chains = 2, 
             seed = 42, control = list(max_treedepth = 15))

# print model results
print(fit1)

# create traceplots of betas to check convergence
traceplot(fit1, pars = c('beta'))

# check coefficient values
rstan::stan_plot(fit1, pars = c('beta'))
stan_hist(fit1, pars = c('beta'))

# extract samples
dfSamples <- as.data.frame(fit1)

# manipulate data to get to tidy samples
dfPredictions <- dfSamples[,13:((nrow(dfTest)*6)+12)] %>%
  mutate(fake_pivot = 1) %>%
  pivot_longer(cols = -fake_pivot) %>%
  mutate(rowid = as.numeric(str_sub(name, 7, (str_length(name)-3))),
         sample_type = str_sub(name, (str_length(name)-1), (str_length(name)-1)),
         sample_type_title = str_c('Type', sample_type, sep = ' ')) %>%
  select(rowid, sample_type, sample_type_title, value)

# join in actual type
dfAccuracy <- dfTest %>%
  rowid_to_column() %>%
  select(rowid, type, ln_temp, mag) %>%
  inner_join(dfPredictions, by = 'rowid') %>%
  mutate(accurate = case_when(
    type == sample_type ~ TRUE,
    TRUE ~ FALSE
  ))

# graph distributions
dfAccuracy %>%
  filter(rowid <= 12) %>%
  ggplot(aes(x = value, fill = accurate)) +
  geom_histogram() +
  scale_y_continuous(breaks = c(0,500,1000)) +
  facet_grid(rowid~sample_type_title) +
  labs(title = 'Predicted Star Types',
       subtitle = 'By Test Star',
       x = 'Likelihood of Type',
       y = 'Count of Samples',
       fill = 'Accurate',
       caption = 'Test Stars 1-12')

dfAccuracy %>%
  filter(rowid > 12 & rowid <= 24) %>%
  ggplot(aes(x = value, fill = accurate)) +
  geom_histogram() +
  scale_y_continuous(breaks = c(0,500,1000)) +
  facet_grid(rowid~sample_type_title) +
  labs(title = 'Predicted Star Types',
       subtitle = 'By Test Star',
       x = 'Likelihood of Type',
       y = 'Count of Samples',
       fill = 'Accurate',
       caption = 'Test Stars 13-24')

dfAccuracy %>%
  filter(rowid > 25 & rowid <= 36) %>%
  ggplot(aes(x = value, fill = accurate)) +
  geom_histogram() +
  scale_y_continuous(breaks = c(0,500,1000)) +
  facet_grid(rowid~sample_type_title) +
  labs(title = 'Predicted Star Types',
       subtitle = 'By Test Star',
       x = 'Likelihood of Type',
       y = 'Count of Samples',
       fill = 'Accurate',
       caption = 'Test Stars 25-36')

dfAccuracy %>%
  filter(rowid > 36) %>%
  ggplot(aes(x = value, fill = accurate)) +
  geom_histogram() +
  scale_y_continuous(breaks = c(0,500,1000)) +
  facet_grid(rowid~sample_type_title) +
  labs(title = 'Predicted Star Types',
       subtitle = 'By Test Star',
       x = 'Likelihood of Type',
       y = 'Count of Samples',
       fill = 'Accurate',
       caption = 'Test Stars 37-48')

# create summary accuracy
dfAccuracySummary <- dfAccuracy %>%
  group_by(rowid, type, ln_temp, mag, sample_type, accurate) %>%
  summarize(mean_value = mean(value)) %>%
  ungroup() %>%
  arrange(rowid, desc(mean_value)) %>%
  group_by(rowid) %>%
  slice(1)

# calculate percent accurate
dfAccuracySummary %>%
  group_by(accurate) %>%
  summarize(count = n())

# calculate percent accurate by true type
dfAccuracySummary %>%
  group_by(type, accurate) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  group_by(type) %>%
  mutate(total = sum(count),
         pct = count/total,
         star_type_title = as.factor(case_when(
           type == 1 ~ '1. Brown Dwarf',
           type == 2 ~ '2. Red Dwarf',
           type == 3 ~ '3. White Dwarf',
           type == 4 ~ '4. Main Sequence',
           type == 5 ~ '5. Supergiant',
           type == 6 ~ '6. Hypergiant'
         ))) %>%
  filter(accurate == TRUE) %>%
  ggplot(aes(x = star_type_title, y = pct, fill = star_type_title)) +
  geom_col() +
  geom_text(aes(label = round(pct, 2)), vjust = -0.25) +
  labs(title = 'Percent Accurate by Actual Star Type',
       x = 'Star Type',
       y = 'Percent Accurately Classified') +
  theme(legend.position = 'none')

# plot on original scale
dfAccuracySummary %>%
  ggplot(aes(x = ln_temp, y = mag, color = accurate)) +
  geom_point() + 
  labs(title = 'Assessing Accuracy of Star Type Classification',
       x = 'Log Temperature (K)',
       y = 'Absolute Magnitude (Mv)',
       color = 'Accurate')


