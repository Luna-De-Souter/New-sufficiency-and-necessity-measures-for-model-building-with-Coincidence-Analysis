### This script creates a plot that illustrates how PA-consistency converges 
### to 0.5 with increasing number of cases when antecedent and outcome are 
### independent, regardless of the prevalence of the outcome. The script also 
### creates a plot that illustrates how AAC-consistency converges to 0.5 with 
### increasing number of cases when antecedent and outcome are independent, 
### regardless of the relative frequency of the antecedent.


# set seed to ensure replicability
set.seed(552)

# load libraries
library(ggplot2)

# define a function that calculates PA-consistency
calculate_PA_con <- function(Phi, Y){
  sum(Phi*Y)/(sum(Phi*Y) + sum(Y)/sum(1-Y)*sum(Phi*(1-Y)))
}

# set the number of cases to converge to
n <- 10000
# set the number of cases per interval
interval <- 100
# calculate the number of intervals
num_intervals <- n/interval

## calculate PA-consistencies for prevalence of 0.2
# generate antecedent values for all cases
Phi <- rbinom(n, 1, 0.4)
# generate outcome values for all cases, with an expected prevalence of 0.2
Y <- rbinom(n, 1, 0.2)
# generate a list containing the PA-consistencies of the first 100, 200, ..., n cases
PA_cons_02 <- rep(NA, num_intervals)
for (i in 1:num_intervals){
  end <- i*interval
  PA_cons_02[i] <- calculate_PA_con(Phi[1:end], Y[1:end])
}

## analogously calculate PA-consistencies for prevalence of 0.5
Phi <- rbinom(n, 1, 0.4)
Y <- rbinom(n, 1, 0.5)
PA_cons_05 <- rep(NA, num_intervals)
for (i in 1:num_intervals){
  end <- i*interval
  PA_cons_05[i] <- calculate_PA_con(Phi[1:end], Y[1:end])
}

## analogously calculate PA-consistencies for prevalence of 0.8
Phi <- rbinom(n, 1, 0.4)
Y <- rbinom(n, 1, 0.8)
PA_cons_08 <- rep(NA, num_intervals)
for (i in 1:num_intervals){
  end <- i*interval
  PA_cons_08[i] <- calculate_PA_con(Phi[1:end], Y[1:end])
}

#create a convergence plot of these three series
df <- data.frame(interval=1:num_intervals, PA_con02=PA_cons_02, PA_con05=PA_cons_05, PA_con08=PA_cons_08)
ggplot(df, aes(x=interval)) +
  geom_line(aes(y=PA_con02, color='Prevalence 0.2')) +
  geom_line(aes(y=PA_con05, color='Prevalence 0.5')) +
  geom_line(aes(y=PA_con08, color='Prevalence 0.8')) +
  geom_point(aes(y=PA_con02, color='Prevalence 0.2'), size=0.5) +
  geom_point(aes(y=PA_con05, color='Prevalence 0.5'), size=0.5) +
  geom_point(aes(y=PA_con08, color='Prevalence 0.8'), size=0.5) +
  labs(x='Number of cases (in hundreds)', y='PA-consistency', title='Convergence of PA-consistency') +
  scale_y_continuous(breaks=seq(0, 1, by=0.1), limits=c(0.3,0.7)) +
  scale_x_continuous(breaks=seq(0, 100, by=20)) +
  scale_color_manual(values=c('Prevalence 0.2'='blue', 'Prevalence 0.5'='red', 'Prevalence 0.8'='green')) +
  theme_minimal() +
  theme(legend.title=element_blank(), legend.position=c(0.8, 0.2))




# define a function that calculates AAC-consistency
calculate_AAC_con <- function(Phi, Y){
  sum((1-Phi)*(1-Y))/(sum((1-Phi)*(1-Y)) + sum(1-Phi)/sum(Phi)*sum(Phi*(1-Y)))
}

# set the number of cases to converge to
n <- 10000
# set the number of cases per interval
interval <- 100
# calculate the number of intervals
num_intervals <- n/interval

## calculate AAC-consistencies for prevalence of 0.2
# generate antecedent values for all cases
Phi <- rbinom(n, 1, 0.2)
# generate outcome values for all cases, with an expected prevalence of 0.2
Y <- rbinom(n, 1, 0.3)
# generate a list containing the AAC-consistencies of the first 100, 200, ..., n cases
AAC_cons_02 <- rep(NA, num_intervals)
for (i in 1:num_intervals){
  end <- i*interval
  AAC_cons_02[i] <- calculate_AAC_con(Phi[1:end], Y[1:end])
}

## analogously calculate AAC-consistencies for prevalence of 0.5
Phi <- rbinom(n, 1, 0.5)
Y <- rbinom(n, 1, 0.3)
AAC_cons_05 <- rep(NA, num_intervals)
for (i in 1:num_intervals){
  end <- i*interval
  AAC_cons_05[i] <- calculate_AAC_con(Phi[1:end], Y[1:end])
}

## analogously calculate AAC-consistencies for prevalence of 0.8
Phi <- rbinom(n, 1, 0.8)
Y <- rbinom(n, 1, 0.3)
AAC_cons_08 <- rep(NA, num_intervals)
for (i in 1:num_intervals){
  end <- i*interval
  AAC_cons_08[i] <- calculate_AAC_con(Phi[1:end], Y[1:end])
}

#create a convergence plot of these three series
df <- data.frame(interval=1:num_intervals, AAC_con02=AAC_cons_02, AAC_con05=AAC_cons_05, AAC_con08=AAC_cons_08)
ggplot(df, aes(x=interval)) +
  geom_line(aes(y=AAC_con02, color='Antecedent relative frequency 0.2')) +
  geom_line(aes(y=AAC_con05, color='Antecedent relative frequency 0.5')) +
  geom_line(aes(y=AAC_con08, color='Antecedent relative frequency 0.8')) +
  geom_point(aes(y=AAC_con02, color='Antecedent relative frequency 0.2'), size=0.5) +
  geom_point(aes(y=AAC_con05, color='Antecedent relative frequency 0.5'), size=0.5) +
  geom_point(aes(y=AAC_con08, color='Antecedent relative frequency 0.8'), size=0.5) +
  labs(x='Number of cases (in hundreds)', y='AAC-consistency', title='Convergence of AAC-consistency') +
  scale_y_continuous(breaks=seq(0, 1, by=0.1), limits=c(0.3,0.7)) +
  scale_x_continuous(breaks=seq(0, 100, by=20)) +
  scale_color_manual(values=c('Antecedent relative frequency 0.2'='blue', 'Antecedent relative frequency 0.5'='red', 'Antecedent relative frequency 0.8'='green')) +
  theme_minimal() +
  theme(legend.title=element_blank(), legend.position=c(0.8, 0.2))


