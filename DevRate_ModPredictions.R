

library(tidyverse)
library(lme4)

###
# Create an example data frame ====
###

EXAMPLEDATA <- expand.grid(Pond = LETTERS[1:10], Clutch = 1:10, ID = 1:6, day = seq(1, 13, by = 2)) %>%
  mutate(stage = 1 + log(day)*10) 

head(EXAMPLEDATA) # This is what your input data should look like. "Day 1" is the stage at time of stocking

# Here we generate some structured variance and errors, Then we use these to make more realistic stage estimates to fit in the model
INPUT <- EXAMPLEDATA %>%
  left_join(data.frame(Pond = unique(EXAMPLEDATA$Pond), B_pond = rnorm(length(unique(EXAMPLEDATA$Pond)), 1, 0.2), b_pond = rnorm(length(unique(EXAMPLEDATA$Pond)), 0, 1)),
            by = "Pond") %>%
  left_join(expand.grid(Pond = unique(EXAMPLEDATA$Pond), Clutch = unique(EXAMPLEDATA$Clutch)) %>% mutate(B_clutch = rnorm(nrow(.), 1, 0.1), b_clutch = rnorm(nrow(.), 0, 0.2)),
            by = c("Pond", "Clutch")) %>%
  left_join(expand.grid(Pond = unique(EXAMPLEDATA$Pond), Clutch = unique(EXAMPLEDATA$Clutch), ID = unique(EXAMPLEDATA$ID)) %>% mutate(B_ID = rnorm(nrow(.), 1, 0.01), b_ID = 0),
            by = c("Pond", "Clutch", "ID")) %>%
  mutate(error = rnorm(nrow(.))) %>%
  mutate(Stage = 7*log(day)*B_pond*B_clutch*B_ID + 5 + b_pond + b_clutch + b_ID + error)

INPUT %>%
  ggplot(aes(x = day, y = Stage, group = interaction(Pond, Clutch, ID))) +
  geom_line() +
  theme_bw() # Each embryo has it's own growth trajectory



###
# Fit the developmental rate model ====
###

# First, we transform the days to the natural log of days
INPUT <- INPUT %>% mutate(ln.day = log(day))

mod_dev <- lmer(Stage ~ ln.day + (0 + ln.day | Pond:Clutch:ID) + (0 + ln.day | Pond:Clutch) + (0 + ln.day | Pond) + (1 | Pond:Clutch) + (1 | Pond), data = INPUT) # This model assumes that all embryos in the same clutch were fertilized at the same time (i.e. no random intercept for ID). Otherwise it assumes correlated errors for individuals within clutch and within ponds in both fertilizations time and developmental rate.

summary(mod_dev) # The main effect for ln.day should be close to 7 and the intercept should be close to 5, which is what we specified in the example dataset before including the random effects



###
# Estimate embryonic period ====
###

# The efficient way to do this would be to back solve the formula so that Day is on the left hand side--then we could simple put in stage of 1 or 20 and solve for the day for each individual. 
# BUT, that is beyond my brainpower right now to figure out how to do that with the random effect model.
# So, my inefficient solution is to just make a massive prediction data frame and then find the Day (within 1/10th of a Day) when the embryo is predicted to be at stage 1 and 20.

EP.calc <- function(dat, lmer.mod){
new.dat1 <- dat %>% group_by(Pond, Clutch, ID) %>% tally() %>% select(-n) %>% mutate(Name = paste(Pond, Clutch, ID, sep = "_")) #Get unique IDs
new.dat2 <- new.dat1 %>%
  left_join(expand.grid(Name = new.dat1$Name, ln.day = log(seq(0.01, max(dat$day) + 4, by = 0.01))),
            by = "Name") # This creates a massive prediction dataframe
  
new.dat2$fit <- predict(lmer.mod, newdata = new.dat2, re.form = NULL) # This predicts the stage for all days. It takes quite a while...

new.dat2 %>% mutate(fit = abs(fit - 1)) %>% filter(Pond == "A" & Clutch == 1) %>% ggplot(aes(x = exp(ln.day), y = fit)) + geom_point()

EP <- new.dat2 %>%
  mutate(fit.stage = abs(fit -1)) %>% # Find the day at stage 1. This makes the closest fit value to stage 1 the absolute minimum
  group_by(Pond, Clutch, ID) %>%
  filter(fit.stage == min(fit.stage)) %>%
  mutate(Day.GS01 = exp(ln.day)) %>%
  select(-fit, -fit.stage, -ln.day, -Name) %>%
  left_join(new.dat2 %>%
              mutate(fit.stage = abs(fit -20)) %>% # Find the day at stage 20. This makes the closest fit value to stage 20 the absolute minimum
              group_by(Pond, Clutch, ID) %>%
              filter(fit.stage == min(fit.stage)) %>%
              mutate(Day.GS20 = exp(ln.day)) %>%
              select(-fit, -fit.stage, -ln.day, -Name),
            by = c("Pond", "Clutch", "ID")) %>%
  mutate(EP = Day.GS20 - Day.GS01) %>% # Estimate the embryonic period as the time between the day predicted for GS 20 and GS1
  select(-Day.GS20, -Day.GS01)

EP
} # This function estimates the embryonic period from stage 1 to 20 in days.

T0 <- Sys.time()
OUTPUT <- EP.calc(INPUT, mod_dev)
T1 <- Sys.time()
T1 - T0 # Time this function.