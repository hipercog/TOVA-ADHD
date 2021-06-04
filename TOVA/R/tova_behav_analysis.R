#Load data
tova_parsed <- read.xlsx("/Users/jpalomak/Downloads/tova_parsed__2013-08-22_12-15-08.xlsx", 1)

#Wrangle data
tova_wrangled <- tova_parsed %>% select(SESNUM, PatientHealthy, NAME, GroupName, GENDER,
                                        CORRSPQ1, CORRSPQ2, CORRSPQ3, CORRSPQ4,
                                        COMERRQ1, COMERRQ2, COMERRQ3, COMERRQ4,
                                        OMERRQ1, OMERRQ2, OMERRQ3, OMERRQ4) %>% 
  filter(SESNUM == 1) %>% rename(QA1 = CORRSPQ1, QA2 = CORRSPQ2, QA3 = CORRSPQ3, QA4 = CORRSPQ4, 
                                 QB1 = COMERRQ1, QB2 = COMERRQ2, QB3 = COMERRQ3, QB4 = COMERRQ4,
                                 QC1 = OMERRQ1, QC2 = OMERRQ2, QC3 = OMERRQ3, QC4 = OMERRQ4)

         
tova_wrangled <- tova_wrangled %>% gather(var, value, QA1:QC4) %>%
  separate(var, into = c("var", "dv", "quarter"), 1:2) %>%
  unite(var, var, quarter, sep="") %>%
  spread(dv, value) %>% rename(quarter = var, corrsp = A, comerr = B, omerr = C) %>%
  mutate(prop_err = (comerr)/(corrsp+comerr),
         total_presses = corrsp+comerr,
         corrsp_fixed = log(max(corrsp+1) - corrsp)) #fixed skewed distribution (comerr and total_presses are helplessly skewed.)

tova_wrangled$NAME <- as.factor(tova_wrangled$NAME)
tova_wrangled$quarter <- as.factor(tova_wrangled$quarter)
tova_wrangled$GroupName <- as.factor(tova_wrangled$GroupName)
tova_wrangled$PatientHealthy <- as.factor(tova_wrangled$PatientHealthy)

#Visualize total number of presses between groups and test quarters
total_presses <- ggplot(tova_wrangled, aes(quarter, total_presses, color=PatientHealthy)) +
  geom_boxplot(position=position_dodge(.9)) +
  #geom_jitter(position = position_jitterdodge(jitter.width = .15, jitter.height = .05, dodge.width=.9), alpha = 0.30) + 
  theme_bw() + labs(title="Total number of presses", x=NULL, y=NULL) + theme(legend.position = "bottom") +
  geom_vline(xintercept=1.5, linetype=2, alpha=.5) +
  geom_vline(xintercept=2.5, linetype=2, alpha=.5) +
  geom_vline(xintercept=3.5, linetype=2, alpha=.5)

#Visualize number of commission errors between groups and test quarters
commission_errors <- ggplot(tova_wrangled, aes(quarter, comerr, color=PatientHealthy)) + 
  geom_boxplot(position=position_dodge(.9)) +
  #geom_jitter(position = position_jitterdodge(jitter.width = .15, jitter.height = .05, dodge.width=.9), alpha = 0.30) + 
  theme_bw() + labs(title="Commission errors", x=NULL, y=NULL) + theme(legend.position = "bottom") +
  geom_vline(xintercept=1.5, linetype=2, alpha=.5) +
  geom_vline(xintercept=2.5, linetype=2, alpha=.5) +
  geom_vline(xintercept=3.5, linetype=2, alpha=.5)

#Visualize number of correct responses between groups and test quarters
correct_responses <- ggplot(tova_wrangled, aes(quarter, corrsp, color=PatientHealthy))  +
  geom_boxplot(position=position_dodge(.9)) +
  #geom_jitter(position = position_jitterdodge(jitter.width = .15, jitter.height = 0.5, dodge.width=.9), alpha = 0.30) + 
  theme_bw() + labs(title="Correct responses", x=NULL, y=NULL) + theme(legend.position = "bottom") +
  geom_vline(xintercept=1.5, linetype=2, alpha=.5) +
  geom_vline(xintercept=2.5, linetype=2, alpha=.5) +
  geom_vline(xintercept=3.5, linetype=2, alpha=.5)

#Visualize proportion of commission errors between groups and test quarters
proportion_commission <- ggplot(tova_wrangled, aes(quarter, prop_err, color=PatientHealthy)) +
  geom_boxplot(position=position_dodge(.9)) +
  #geom_jitter(position = position_jitterdodge(jitter.width = .15, jitter.height = 0, dodge.width=.9), alpha = 0.30) + 
  theme_bw() + labs(title="Proportion of commission errors", x=NULL, y=NULL) + theme(legend.position = "bottom") +
  geom_vline(xintercept=1.5, linetype=2, alpha=.5) +
  geom_vline(xintercept=2.5, linetype=2, alpha=.5) +
  geom_vline(xintercept=3.5, linetype=2, alpha=.5)

ggarrange(total_presses, commission_errors, correct_responses, proportion_commission, 
          ncol=2, nrow=2, common.legend = TRUE, legend="bottom")


#ADHD vs control for correct responses (fixed distribution)
tova_model1 <- lmer(corrsp_fixed ~ PatientHealthy*quarter + (1|NAME), data=tova_wrangled)

tova_effects <- effect(c("PatientHealthy*quarter"), tova_model1)
tova_effects <- as.data.frame(tova_effects)

ggplot(tova_effects, aes(PatientHealthy, fit, fill=quarter)) + geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=lower, ymax=upper), position=position_dodge(.9), width=.2) +
  theme_bw() + labs(y="log(max(corrsp+1) - corrsp)", title="Model prediction") + theme(legend.position="bottom")


