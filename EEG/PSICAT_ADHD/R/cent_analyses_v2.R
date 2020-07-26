#Näitä ennen pitää ajaa exgauss.rmd -stä Tuiskin alustavat aineiston wranglaukset

sbj_params_ml$all_groups <- interaction(sbj_params_ml$congruency,sbj_params_ml$shape,sbj_params_ml$group)
levels(sbj_params_ml$all_groups)
relevel(sbj_params_ml$all_groups, ref = "Congruent.shape.ctrl")
#sbj_params_ml$all_groups = factor(sbj_params_ml$all_groups,levels(sbj_params_ml$all_groups)[c(1:8)])
mtau.log_new <- lmer(log(tau) ~ all_groups + (1|Subject), data=sbj_params_ml)
summary(mtau.log_new)

#Collect effects for plotting, simple model
ef_cent_lmm <- effect("group:shape:congruency", mtau.log)
x_cent_lmm <- as.data.frame(ef_cent_lmm)
x_cent_lmm

ggplot(x_cent_lmm, aes(x=congruency, y=fit, fill=shape)) + 
  geom_bar(position=position_dodge(), stat="identity", colour="black", size=.2) +
  geom_errorbar(aes(ymin=lower, ymax=upper),width=.2, position=position_dodge(.9)) + 
  xlab("Congruency") + ylab("Rtime") + theme_set(theme_gray(base_size = 12)) + facet_wrap("group") + coord_cartesian(ylim = c(3.7, 5)) +
  scale_fill_hue(name="Shape")

cbind(new_data, rtime = predict(mtau.log, newdata = new_data, re.form=NA)) %>%
  ggplot(aes(congruency, exp(rtime), linetype=shape)) + geom_point() + 
  geom_line(aes(group=shape),size=1) +
  facet_wrap(~group) + theme_bw() + ggtitle("tau values (log model)")

ggplot(sbj_params_ml) + aes(all_groups, log(tau)) + geom_point() + theme_bw() + ggtitle("tau values") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#difflsmeans(mtau.log, c("congruency", "shape", "group"))

#Lasketaan kolmisuuntainen interaktio kontrastien avulla:

#"MAIN EFFECT" of congruency within ADHD
contrast_congr_adhd <- rbind("contrast" = c(0,0,0,0,1,-1,1,-1))
summary(glht(mtau.log_new, linfct = mcp(all_groups = contrast_congr_adhd)))

#"MAIN EFFECT" of shape within ADHD
contrast_shape_adhd <- rbind("contrast" = c(0,0,0,0,1,1,-1,-1))
summary(glht(mtau.log_new, linfct = mcp(all_groups = contrast_shape_adhd)))

#"MAIN EFFECT" of congruency within CTRL
contrast_congr_ctrl <- rbind("contrast" = c(1,-1,1,-1,0,0,0,0))
summary(glht(mtau.log_new, linfct = mcp(all_groups = contrast_congr_ctrl)))

#"MAIN EFFECT" of shape within CTRL
contrast_shape_ctrl <- rbind("contrast" = c(1,1,-1,-1,0,0,0,0))
summary(glht(mtau.log_new, linfct = mcp(all_groups = contrast_shape_ctrl)))

#interaction WITHIN ADHD
comparison_adhd <- (contrast_congr_adhd*contrast_shape_adhd)
summary(glht(mtau.log_new, linfct = mcp(all_groups = comparison_adhd)))

#interaction WITHIN CTRL
comparison_ctrl <- (contrast_congr_ctrl*contrast_shape_ctrl)
summary(glht(mtau.log_new, linfct = mcp(all_groups = comparison_ctrl)))

#Kolmisuuntainen interaktio on sama asia, kuin WITHIN ADHD ja WITHIN CTRL -interaktioiden erotus:
#difference between interactions adhd vs control (equals three-way interaction group*shape*congruency)
comparison_final <- (comparison_adhd - comparison_ctrl)
summary(glht(mtau.log_new, linfct = mcp(all_groups = comparison_final)))

#Interaktion tulkinta: IVn yhteys DVhen erilainen RIIPPUEN jonkin toisen IVn tasosta. Kolmisuuntaisen interaktion tulkinta:
#Interaktio käyttäytyy eri tavalla riippuen jonkin kolmannen IVn tasosta.

#Seuraavaksi, Benin spesifi tutkimuskysymys oli: Onko congruency -faktorin efekti suurempi kuin shape -faktorin efekti ADHD- CTRL- ryhmien sisällä,
#ja, erityisesti: Onko näiden EFEKTIEROJEN välillä eroja adhd- ja ctrl -ryhmien välillä. Tämä selvitetään seuraavilla kontrasteilla:

difference1 <- contrast_congr_ctrl - contrast_shape_ctrl #congruency- ja shape -faktorien efektierot WITHIN ctrl
difference2 <- contrast_congr_adhd - contrast_shape_adhd #congruency- ja shape -faktorien efektierot WITHIN adhd
difference3 <- difference1 - difference2 #näiden välinen ero
summary(glht(mtau.log_new, linfct = mcp(all_groups = difference3)))




#------ alla unsorted kamaa

#contrast for incongruent_shape_adhd vs incongruent_nonshape_adhd:
contrast1 <- rbind("contrast" = c(0,0,0,0,0,1,0,-1))
summary(glht(mtau.log_new, linfct = mcp(all_groups = contrast1)))

#contrast for incongruent_shape_ctr vs incongruent_nonshape_ctr
contrast2 <- rbind("contrast" = c(0,1,0,-1,0,0,0,0))
summary(glht(mtau.log_new, linfct = mcp(all_groups = contrast2)))

#contrast difference
(contrast1minus2 <- contrast1 - contrast2)
summary(glht(mtau.log_new, linfct = mcp(all_groups = contrast1minus2)))

#Syy miksi ylla ei saada 3-way interaktiona merkitsevaksi, on se, etta siina on oikeastaan laskettu eri interaktio!

contrast_congr_shape_adhd <- rbind("contrast" = c(0,0,0,0,1,-1,0,0))
summary(glht(mtau.log_new, linfct = mcp(all_groups = contrast_congr_shape_adhd)))

contrast_congr_shape_ctrl <- rbind("contrast" = c(1,-1,0,0,0,0,0,0))
summary(glht(mtau.log_new, linfct = mcp(all_groups = contrast_congr_shape_ctrl)))

comparison_shape_final <- (contrast_congr_shape_adhd - contrast_congr_shape_ctrl)
summary(glht(mtau.log_new, linfct = mcp(all_groups = comparison_shape_final)))


#alla oleva tuottaa about samat tulokset kuin lmer-malli
# contrast_main_congr <- rbind("contrast" = c(1,-1,1,-1,1,-1,1,-1))
# summary(glht(mtau.log_new, linfct = mcp(all_groups = contrast_main_congr)))
# 
# contrast_main_shape <- rbind("contrast" = c(1,1,-1,-1,1,1,-1,-1))
# summary(glht(mtau.log_new, linfct = mcp(all_groups = contrast_main_shape)))

# contrast_main_group <- rbind("contrast" = c(1,1,1,1,-1,-1,-1,-1))
# summary(glht(mtau.log_new, linfct = mcp(all_groups = contrast_main_group)))
# 
# summary(glht(mtau.log_new, linfct = mcp(all_groups = contrast_main_shape*contrast_main_congr*contrast_main_group)))


#So:
summary(glht(mtau.log_new, linfct = mcp(all_groups = contrast_main_congr*contrast_main_group)))

#equals:
summary(glht(mtau.log_new, linfct = mcp(all_groups = contrast_congr_ctrl-contrast_congr_adhd)))



  
