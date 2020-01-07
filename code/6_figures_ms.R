# 6: Figures for manuscript
# Comparison of SDM approaches
# Tim Szewczyk

# This script plots the output from `3_fitModels.R` and `4_analysis.R`.

########
## Make appendixes
########

rmarkdown::render("ms/supp/App_1.Rmd", 
                  output_file="App_1.pdf",
                  knit_root_dir="../")

rmarkdown::render("ms/supp/App_2.Rmd", 
                  output_file="App_2.pdf",
                  knit_root_dir="../")



########
## Setup
########

# load workspace
pkgs <- c("tidyverse", "magrittr", "stringr", "here", "viridis",
          "grid", "gtable", "ggnewscale")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "fn", full.names=T), source)
theme_set(theme_bw())
sp.names <- c("shrub", "biennial")
SDM_col <- c(MaxEnt="#fd8d3c", IPM="#08306b", 
             'CA[i]'="#2171b5", 'CA[p]'="#6baed6")

# TSS
TSS.df <- rbind(read.csv("out/sp1/out_TSS.csv"),
                read.csv("out/sp2/out_TSS.csv")) %>% 
  rbind(., expand.grid(TSS=NA, sensitivity=NA, specificity=NA,
                       Boundary=unique(.$Boundary),
                       SDM="MxE", issue=c("noSB", "underDisp", "overDisp"),
                       sp=unique(.$sp))) 
TSS.df$issue <- factor(TSS.df$issue, 
                       levels=c("none", "noise", "sampBias", "nonEq", 
                                "wrongCov", "noSB", "underDisp", "overDisp"),
                       labels=c("Ideal", "Measurement error", "Sampling bias", 
                                "Non equilibrium", "Incorrect covariates", 
                                "No Seedbank", "Under dispersal",
                                "Over dispersal"))
TSS.df$sp <- factor(TSS.df$sp, levels=c("sp1", "sp2"), labels=sp.names)
TSS.df$Boundary <- paste0("True~range:~", TSS.df$Boundary)
TSS.df$SDM <- lvls_revalue(TSS.df$SDM, c("CA[p]", "CA[i]", "IPM", "MaxEnt"))
TSS.ci <- TSS.df %>% group_by(sp, SDM, issue, Boundary) %>% 
  summarise(med=median(TSS, na.rm=T), 
            mn=mean(TSS, na.rm=T),
            sd=sd(TSS, na.rm=T),
            q025=quantile(TSS, probs=0.025, na.rm=T), 
            q25=quantile(TSS, probs=0.25, na.rm=T), 
            q75=quantile(TSS, probs=0.75, na.rm=T), 
            q975=quantile(TSS, probs=0.975, na.rm=T)) %>%
  arrange(Boundary, issue, sp, desc(med)) %>%
  ungroup %>% group_by(Boundary, issue, sp) %>%
  mutate(rank_SDM=row_number())  %>%
  full_join(., filter(., issue=="Ideal"), by=c("sp", "SDM", "Boundary"), 
            suffix=c("", ".none"))

# Sensitivity & specificity
senspe.ci <- TSS.df %>% group_by(sp, SDM, issue, Boundary) %>% 
  summarise(med.sen=median(sensitivity, na.rm=T), 
            mn.sen=mean(sensitivity, na.rm=T),
            q025.sen=quantile(sensitivity, probs=0.025, na.rm=T), 
            q25.sen=quantile(sensitivity, probs=0.25, na.rm=T), 
            q75.sen=quantile(sensitivity, probs=0.75, na.rm=T), 
            q975.sen=quantile(sensitivity, probs=0.975, na.rm=T),
            med.spe=median(specificity, na.rm=T), 
            mn.spe=mean(specificity, na.rm=T),
            q025.spe=quantile(specificity, probs=0.025, na.rm=T), 
            q25.spe=quantile(specificity, probs=0.25, na.rm=T), 
            q75.spe=quantile(specificity, probs=0.75, na.rm=T), 
            q975.spe=quantile(specificity, probs=0.975, na.rm=T)) %>%
  arrange(Boundary, issue, sp, desc(med.sen)) %>%
  full_join(., filter(., issue=="Ideal"), by=c("sp", "SDM", "Boundary"), 
            suffix=c("", ".none"))

# Log likelihood
LL.df <- rbind(read.csv("out/sp1/out_LL.csv"),
                read.csv("out/sp2/out_LL.csv")) %>% 
  rbind(., expand.grid(LogLik=NA, Boundary=unique(.$Boundary),
                       SDM="MxE", issue=c("noSB", "underDisp", "overDisp"),
                       sp=unique(.$sp))) 
LL.df$issue <- factor(LL.df$issue, 
                       levels=c("none", "noise", "sampBias", "nonEq", 
                                "wrongCov", "noSB", "underDisp", "overDisp"),
                       labels=c("Ideal", "Measurement error", "Sampling bias", 
                                "Non equilibrium", "Incorrect covariates", 
                                "No Seedbank", "Under dispersal",
                                "Over dispersal"))
LL.df$sp <- factor(LL.df$sp, levels=c("sp1", "sp2"), labels=sp.names)
LL.df$Boundary <- paste0("True~range:~", LL.df$Boundary)
LL.df$SDM <- lvls_revalue(LL.df$SDM, c("CA[p]", "CA[i]", "IPM", "MaxEnt"))

# Predictions
P.df <- rbind(read.csv("out/sp1/out_P.csv"),
                read.csv("out/sp2/out_P.csv"))
P.df$CI_95 <- NA
P.df$CI_95[P.df$SDM=="IPM"] <- P.df$lambda.975.f[P.df$SDM=="IPM"] - 
  P.df$lambda.025.f[P.df$SDM=="IPM"]
P.df$CI_95[grepl("CA", P.df$SDM)] <- P.df$Surv.S.975.f[grepl("CA", P.df$SDM)] - 
  P.df$Surv.S.025.f[grepl("CA", P.df$SDM)]
P.df$S_in_CI <- NA
P.df$S_in_CI[P.df$SDM!="MxE"] <- FALSE
P.df$S_in_CI[grepl("CA", P.df$SDM) & P.df$Surv.S.975.f >= P.df$Surv.S & 
               P.df$Surv.S.025.f <= P.df$Surv.S] <- TRUE
P.df$S_in_CI[P.df$SDM=="IPM" & P.df$lambda.975.f >= P.df$lambda & 
               P.df$lambda.025.f <= P.df$lambda] <- TRUE
P.df$issue <- factor(P.df$issue, 
                     levels=c("none", "noise", "sampBias", "nonEq", 
                              "wrongCov", "noSB", "underDisp", "overDisp"),
                     labels=c("Ideal", "Measurement error", "Sampling bias", 
                              "Non equilibrium", "Incorrect covariates", 
                              "No Seedbank", "Under dispersal",
                              "Over dispersal"))
P.df$sp <- factor(P.df$sp, levels=c("sp1", "sp2"), labels=sp.names)
P.df$SDM <- lvls_revalue(P.df$SDM, c("CA[p]", "CA[i]", "IPM", "MaxEnt"))
P.sum <- P.df %>% group_by(sp, SDM, issue) %>%
  summarise(propDisputedCells=sum(prP>0 & prP<1)/n(),
            propSinCI=sum(S_in_CI)/n(),
            mean95CI=mean(CI_95, na.rm=T),
            median95CI=median(CI_95, na.rm=T))
P.sum <- rbind(ungroup(P.sum), 
               expand.grid(sp=unique(P.sum$sp), SDM="MaxEnt", 
                           issue=c("No Seedbank", 
                                   "Under dispersal", "Over dispersal"),
                           propDisputedCells=NA, propSinCI=NA,
                           mean95CI=NA, median95CI=NA))






########
## TSS: summaries
########
# all TSS values
summary(TSS.df$TSS)

# all TSS medians
summary(TSS.ci$med)

# SDM rank within scenarios: scenarios with all SDMs only
SDM.ranks <- TSS.ci %>% ungroup %>% 
  filter(!issue %in% c("No Seedbank", "Under dispersal", "Over dispersal")) %>% 
  group_by(SDM, Boundary, sp) %>% 
  summarise(min_rank=min(rank_SDM), mean_rank=mean(rank_SDM), 
            med_rank=median(rank_SDM), max_rank=max(rank_SDM)) %>%
  arrange(Boundary, mean_rank)

# Scenario ranks within SDMs
scenario.ranks <- TSS.ci %>% ungroup %>% 
  # filter(SDM!="MaxEnt") %>%
  # filter(!issue %in% c("No Seedbank", "Under dispersal", "Over dispersal")) %>% 
  arrange(Boundary, SDM, sp, desc(med)) %>%
  group_by(Boundary, SDM, sp) %>%
  mutate(rank_issue=row_number()) %>% 
  ungroup %>% group_by(issue, Boundary) %>% 
  summarise(min_rank=min(rank_issue), mean_rank=mean(rank_issue), 
            med_rank=median(rank_issue), max_rank=max(rank_issue)) %>%
  arrange(Boundary, mean_rank)

# SDM rank within scenarios: scenarios with all SDMs only
senspe.ranks <- senspe.ci %>% ungroup %>% 
  filter(!issue %in% c("No Seedbank", "Under dispersal", "Over dispersal")) %>% 
  arrange(Boundary, issue, sp, desc(med.sen)) %>%
  ungroup %>% group_by(Boundary, issue, sp) %>%
  mutate(rank_sen=row_number()) %>%
  arrange(Boundary, issue, sp, desc(med.spe)) %>%
  ungroup %>% group_by(Boundary, issue, sp) %>%
  mutate(rank_spe=row_number())
# Sensitivity ranks
senspe.ranks %>% ungroup %>% group_by(SDM, Boundary) %>% 
  summarise(min_rank=min(rank_sen), mean_rank=mean(rank_sen), 
            med_rank=median(rank_sen), max_rank=max(rank_sen)) %>%
  arrange(Boundary, mean_rank)
# Specificity ranks
senspe.ranks %>% ungroup %>% group_by(SDM, Boundary) %>% 
  summarise(min_rank=min(rank_spe), mean_rank=mean(rank_spe), 
            med_rank=median(rank_spe), max_rank=max(rank_spe)) %>%
  arrange(Boundary, mean_rank)






########
## SDM ranks
########
# Table 1: Mean ranks (med TSS) across all scenarios and species
TSS.ci %>% ungroup %>% 
  filter(!issue %in% c("No Seedbank", "Under dispersal", "Over dispersal")) %>% 
  group_by(SDM, Boundary) %>% 
  summarise(min_rank=min(rank_SDM), mean_rank=mean(rank_SDM), 
            med_rank=median(rank_SDM), max_rank=max(rank_SDM)) %>%
  arrange(Boundary, mean_rank)

# Supplemental figure: Mean ranks (med TSS) across all scenarios BY species
ggplot(SDM.ranks, aes(x=Boundary, y=mean_rank, colour=SDM, shape=sp)) + 
  geom_hline(yintercept=seq(1, 4, 1), colour="gray90", linetype=2, size=0.2) +
  geom_point(position=position_jitter(width=0.1, height=0)) + 
  scale_colour_manual("SDM\nMethod", values=SDM_col, 
                      labels=c(expression(CA[p], CA[i], IPM, MaxEnt))) +
  scale_x_discrete("True range", labels=expression(lambda>1, N>0)) +
  labs(y="Mean rank across scenarios") +
  theme(panel.grid=element_blank(), 
        legend.text.align=0,
        legend.key.size=unit(0.3, 'cm'))
ggsave("figs/Supp_Fig_Ranks.jpg", width=3, height=3)






########
## No issue
########
tss.none <- filter(TSS.df, issue=="Ideal")
tss.none.ci <- filter(TSS.ci, issue=="Ideal")

TukeyHSD(aov(TSS ~ SDM*sp, 
             data=filter(tss.none, Boundary=="True~range:~lambda > 1")))
TukeyHSD(aov(TSS ~ SDM*sp, 
             data=filter(tss.none, Boundary=="True~range:~N > 0")))

ggplot(tss.none.ci, aes(x=SDM, y=mn, min=mn-.196*sd, ymax=mn+.196*sd)) +
  geom_point() + geom_errorbar() +
  facet_grid(Boundary~sp) + labs(y="Mean TSS ± 1.96*SE")
ggplot(tss.none.ci, aes(x=SDM, y=med, min=q025, ymax=q975)) +
  geom_point() + geom_errorbar() +
  facet_grid(Boundary~sp)

p <- ggplot(filter(P.df, sp=="shrub" & issue!="Ideal"), 
            aes(x=prP.none, y=prP)) + 
  geom_abline(colour="gray70") + geom_point(alpha=0.01) + 
  facet_grid(SDM~issue)
ggsave("~/Desktop/shrub_prP_vIdeal.jpg", p, height=6, width=10)
p <- ggplot(filter(P.df, sp=="biennial" & issue!="Ideal"), 
            aes(x=prP.none, y=prP)) + 
  geom_abline(colour="gray70") + geom_point(alpha=0.01) + 
  facet_grid(SDM~issue)
ggsave("~/Desktop/biennial_prP_vIdeal.jpg", p, height=6, width=10)





########
## TSS/Sensitivity/Specificity: median + 95% CI
########
# MS Figure: TSS
ggplot(TSS.df, aes(x=issue, y=TSS, colour=SDM)) + 
  geom_hline(yintercept=seq(0.6, 1, 0.1), colour="gray90", linetype=2, size=0.2) +
  geom_boxplot(outlier.shape=NA) + ylim(.7, 1) +
  facet_grid(Boundary~., labeller=labeller(Boundary=label_parsed)) +
  scale_colour_manual("SDM", values=SDM_col, 
                      labels=c(expression(CA[p], CA[i], IPM, MaxEnt))) +
  theme(panel.grid=element_blank(), 
        legend.text.align=0,
        axis.text.x=element_text(angle=270, vjust=0.5, hjust=0),
        legend.key.size=unit(0.3, 'cm')) + 
  labs(x="Scenario", y="TSS")
ggsave("figs/Fig_TSS_box.jpg", width=6, height=4)

ggplot(TSS.ci, aes(x=issue, y=med, colour=SDM, group=SDM)) + 
  geom_hline(yintercept=seq(0.4, 1, 0.1), colour="gray90", linetype=2, size=0.2) +
  geom_line(size=0.2, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=q025, ymax=q975), size=0.3, width=0.75, 
                position=position_dodge(width=0.4)) + 
  geom_linerange(aes(ymin=q25, ymax=q75), size=0.75, 
                 position=position_dodge(width=0.4)) + 
  facet_grid(Boundary~sp, labeller=labeller(Boundary=label_parsed)) +
  scale_colour_manual("SDM", values=SDM_col, 
                      labels=c(expression(CA[p], CA[i], IPM, MaxEnt))) +
  theme(panel.grid=element_blank(), 
        legend.text.align=0,
        axis.text.x=element_text(angle=270, vjust=0.5, hjust=0),
        legend.key.size=unit(0.3, 'cm')) + 
  labs(title="TSS across 100 sampled datasets", x="Scenario", y="TSS")
ggsave("figs/Fig_TSS_med.jpg", width=5, height=5)

ggplot(TSS.ci, aes(x=issue, y=mn, colour=SDM, group=SDM)) + 
  geom_hline(yintercept=seq(0.7, 1, 0.1), colour="gray90", linetype=2, size=0.2) +
  geom_line(size=0.2, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=mn-.196*sd, ymax=mn+.196*sd), size=0.3, width=0.75, 
                position=position_dodge(width=0.4)) + 
  facet_grid(Boundary~sp, labeller=labeller(Boundary=label_parsed)) +
  scale_colour_manual("SDM", values=SDM_col, 
                      labels=c(expression(CA[p], CA[i], IPM, MaxEnt))) +
  theme(panel.grid=element_blank(), 
        legend.text.align=0,
        axis.text.x=element_text(angle=270, vjust=0.5, hjust=0),
        legend.key.size=unit(0.3, 'cm')) + 
  labs(#title="TSS across 100 sampled datasets", 
       x="Scenario", y="TSS")
ggsave("figs/Fig_TSS_mn+CI.jpg", width=5, height=5)

ggplot(TSS.ci, aes(x=issue, y=mn, colour=SDM, group=SDM)) + 
  geom_hline(yintercept=seq(0.7, 1, 0.1), colour="gray90", linetype=2, size=0.2) +
  geom_line(size=0.2, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=mn-.196*sd, ymax=mn+.196*sd), size=0.3, width=0.75, 
                position=position_dodge(width=0.4)) + 
  facet_grid(Boundary~sp, labeller=labeller(Boundary=label_parsed)) +
  scale_colour_manual("SDM", values=SDM_col, 
                      labels=c(expression(CA[p], CA[i], IPM, MaxEnt))) +
  theme(panel.grid=element_blank(), 
        legend.text.align=0,
        legend.key.size=unit(0.3, 'cm')) + 
  labs(#title="TSS across 100 sampled datasets", 
    x="Scenario", y="TSS") + coord_flip()
ggsave("figs/Fig_TSS_mn+CI_flip.jpg", width=6, height=4)

ggplot(TSS.ci, aes(x=issue, y=mn, colour=SDM, group=SDM)) + 
  geom_hline(yintercept=seq(0.7, 1, 0.1), colour="gray90", linetype=2, size=0.2) +
  geom_line(size=0.2, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=q25, ymax=q75), size=0.3, width=0.75, 
                position=position_dodge(width=0.4)) + 
  facet_grid(Boundary~sp, labeller=labeller(Boundary=label_parsed)) +
  scale_colour_manual("SDM", values=SDM_col, 
                      labels=c(expression(CA[p], CA[i], IPM, MaxEnt))) +
  theme(panel.grid=element_blank(), 
        legend.text.align=0,
        axis.text.x=element_text(angle=270, vjust=0.5, hjust=0),
        legend.key.size=unit(0.3, 'cm')) + 
  labs(title="TSS across 100 sampled datasets", x="Scenario", y="TSS")
ggsave("figs/Fig_TSS_mn+50.jpg", width=5, height=5)

# Supplemental Figure: Sensitivity
ggplot(senspe.ci, aes(x=issue, y=med.sen, colour=SDM, group=SDM)) + 
  geom_hline(yintercept=seq(0.4, 1, 0.1), colour="gray90", linetype=2, size=0.2) +
  geom_line(size=0.2, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=q025.sen, ymax=q975.sen), size=0.3, width=0.75, 
                position=position_dodge(width=0.4)) + 
  geom_linerange(aes(ymin=q25.sen, ymax=q75.sen), size=0.75, 
                 position=position_dodge(width=0.4)) + 
  facet_grid(Boundary~sp, labeller=labeller(Boundary=label_parsed)) +
  scale_colour_manual("SDM", values=SDM_col, 
                      labels=c(expression(CA[p], CA[i], IPM, MaxEnt))) +
  theme(panel.grid=element_blank(), 
        legend.text.align=0,
        axis.text.x=element_text(angle=270, vjust=0.5, hjust=0),
        legend.key.size=unit(0.3, 'cm')) + 
  labs(title="Sensitivity across 100 sampled datasets", 
       x="Scenario", y="Sensitivity\nPr(Correctly predicted presence)")
ggsave("figs/Supp_Sens_med.jpg", width=5, height=5)

# Supplemental Figure: Specificity
ggplot(senspe.ci, aes(x=issue, y=med.spe, colour=SDM, group=SDM)) + 
  geom_hline(yintercept=seq(0.4, 1, 0.1), colour="gray90", linetype=2, size=0.2) +
  geom_line(size=0.2, position=position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=q025.spe, ymax=q975.spe), size=0.3, width=0.75, 
                position=position_dodge(width=0.4)) + 
  geom_linerange(aes(ymin=q25.spe, ymax=q75.spe), size=0.75, 
                 position=position_dodge(width=0.4)) + 
  facet_grid(Boundary~sp, labeller=labeller(Boundary=label_parsed)) +
  scale_colour_manual("SDM\nMethod", values=SDM_col, 
                      labels=c(expression(CA[p], CA[i], IPM, MaxEnt))) +
  theme(panel.grid=element_blank(), 
        legend.text.align=0,
        axis.text.x=element_text(angle=270, vjust=0.5, hjust=0),
        legend.key.size=unit(0.3, 'cm')) + 
  labs(title="Specificity across 100 sampled datasets", 
       x="Scenario", y="Specificity\nPr(Correctly predicted absence)")
ggsave("figs/Supp_Spec_med.jpg", width=5, height=5)






########
## TSS/Sensitivity/Specificity vs. Ideal
########
TSS.df <- full_join(TSS.df, 
                    filter(TSS.df, issue=="Ideal") %>% 
                      select(TSS, sp, SDM, Boundary), 
                    by=c("sp", "SDM", "Boundary"), 
                      suffix=c("", ".none"))

# MS Figure: TSS vs. none
ggplot(filter(TSS.ci, issue != "Ideal"), 
       aes(x=issue, y=med-med.none, colour=SDM, group=SDM, shape=sp)) + 
  geom_hline(yintercept=seq(-0.1, 0.05, 0.05), colour="gray85", 
             linetype=2, size=0.2) + ylim(-0.12, 0.08) +
  geom_hline(yintercept=0, colour="gray80", linetype=1, size=0.2) +
  geom_point(aes(fill=SDM), position=position_dodge(width=1), alpha=0.5, size=2) +
  geom_point(position=position_dodge(width=1), fill=NA, size=2) +
  facet_grid(Boundary~issue, scales="free", 
             labeller=labeller(Boundary=label_parsed, issue=label_wrap_gen(11))) +
  scale_colour_manual("SDM:", values=SDM_col, 
                      labels=c(expression(CA[p], CA[i], IPM, MaxEnt))) +
  scale_fill_manual("SDM:", values=SDM_col, 
                      labels=c(expression(CA[p], CA[i], IPM, MaxEnt))) +
  scale_shape_manual("Species:", values=c(21,24)) +
  # scale_shape_manual("Species", values=1:2) +
  theme(panel.grid=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.key.size=unit(0.25, 'cm'), 
        legend.box.margin=margin(-20,0,0,0),
        legend.text.align=0,
        legend.position="bottom",
        panel.spacing.x=unit(c(rep(-.1, 6)), "cm"),
        panel.spacing.y=unit(0.1, "cm"),
        panel.border=element_rect(colour="gray40")) + 
  labs(#title="Effect of scenarios on median TSS", 
       x="", 
       y=expression(TSS[i]-TSS['Ideal']))
ggsave("figs/Fig_TSSvIdeal.jpg", width=8, height=3.75)

ggplot(filter(senspe.ci, issue != "Ideal"), 
       aes(x=issue, y=med.sen-med.sen.none, colour=SDM, group=SDM, shape=sp)) + 
  geom_hline(yintercept=seq(-0.15, 0.05, 0.05), colour="gray90",
             linetype=2, size=0.2) +
  geom_hline(yintercept=0, colour="gray80", linetype=1, size=0.2) +
  geom_point(aes(fill=SDM), position=position_dodge(width=1), alpha=0.5, size=2) +
  geom_point(position=position_dodge(width=1), fill=NA, size=2) +
  facet_grid(Boundary~issue, scales="free", 
             labeller=labeller(Boundary=label_parsed, issue=label_wrap_gen(11))) +
  scale_colour_manual("SDM:", values=SDM_col, 
                      labels=c(expression(CA[p], CA[i], IPM, MaxEnt))) +
  scale_fill_manual("SDM:", values=SDM_col, 
                    labels=c(expression(CA[p], CA[i], IPM, MaxEnt))) +
  scale_shape_manual("Species:", values=c(21,24)) +
  theme(panel.grid=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.key.size=unit(0.25, 'cm'), 
        legend.box.margin=margin(-20,0,0,0),
        legend.text.align=0,
        legend.position="bottom",
        panel.spacing.x=unit(c(rep(-.1, 6)), "cm"),
        panel.spacing.y=unit(0.1, "cm"),
        panel.border=element_rect(colour="gray40")) + 
  labs(title="Effect of scenarios on median sensitivity", x="", 
       y=expression(Sensitivity[i]-Sensitivity['ideal']))
ggsave("figs/Supp_SensvIdeal.jpg", width=8, height=3.75)

ggplot(filter(senspe.ci, issue != "Ideal"), 
       aes(x=issue, y=med.spe-med.spe.none, colour=SDM, group=SDM, shape=sp)) + 
  geom_hline(yintercept=seq(-0.15, 0.1, 0.05), colour="gray90",
             linetype=2, size=0.2) +
  geom_hline(yintercept=0, colour="gray80", linetype=1, size=0.2) +
  geom_point(aes(fill=SDM), position=position_dodge(width=1), alpha=0.5, size=2) +
  geom_point(position=position_dodge(width=1), fill=NA, size=2) +
  facet_grid(Boundary~issue, scales="free", 
             labeller=labeller(Boundary=label_parsed, issue=label_wrap_gen(11))) +
  scale_colour_manual("SDM:", values=SDM_col, 
                      labels=c(expression(CA[p], CA[i], IPM, MaxEnt))) +
  scale_fill_manual("SDM:", values=SDM_col, 
                    labels=c(expression(CA[p], CA[i], IPM, MaxEnt))) +
  scale_shape_manual("Species:", values=c(21,24)) +
  theme(panel.grid=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.key.size=unit(0.25, 'cm'), 
        legend.box.margin=margin(-20,0,0,0),
        legend.text.align=0,
        legend.position="bottom",
        panel.spacing.x=unit(c(rep(-.1, 6)), "cm"),
        panel.spacing.y=unit(0.1, "cm"),
        panel.border=element_rect(colour="gray40")) + 
  labs(title="Effect of scenarios on median specificity", x="", 
       y=expression(Specificity[i]-Specificity['ideal']))
ggsave("figs/Supp_SpecvIdeal.jpg", width=8, height=3.75)






########
## Log Likelihood
########
ggplot(LL.df, aes(x=LogLik, y=fct_rev(issue), colour=SDM)) + 
  facet_grid(Boundary~sp, labeller=labeller(Boundary=label_parsed)) +
  geom_point(alpha=0.75, size=3) + 
  scale_colour_manual("SDM\nMethod", values=SDM_col, 
                      labels=c(expression(CA[p], CA[i], IPM, MaxEnt))) +
  theme(legend.text.align=0) +
  labs(title="Log Likelihood", x="Log likelihood: S ~ Binom(prP)", y="")
ggsave("figs/Supp_LogLik.jpg", width=7, height=4)






########
## Prediction maps
########
P.df <- P.df %>% 
  mutate(lam_err=lambda.f-lambda,
         llam_err=log(lambda.f)-log(lambda),
         N_err=Surv.S.f-Surv.S,
         lN_err=log(Surv.S.f+1)-log(Surv.S+1)) %>%
  group_by(SDM, issue, sp) %>%
  mutate(lam_err_Z=scale(lam_err),
         llam_err_Z=scale(llam_err),
         N_err_Z=scale(N_err),
         lN_err_Z=scale(lN_err))
# Lambda
lam.p <- ggplot(filter(P.df, SDM=="IPM")) + 
  geom_tile(aes(lon, lat, fill=llam_err_Z)) + 
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank(),
        panel.spacing.x=unit(c(rep(-.1, 7)), "cm"),
        panel.spacing.y=unit(-0.1, "cm"),
        legend.position="bottom",
        legend.key.height=unit(0.25, "cm")) + 
  scale_fill_gradient2(expression(lambda:Predicted-True~(z-transformed~log-scale)), 
                       midpoint=0) + 
  facet_grid(sp~issue, labeller=labeller(issue=label_wrap_gen(11))) +
  labs(title=expression(Error~'in'~lambda~predictions))
ggsave("figs/map_llam_err_Z.jpg", lam.p, width=9, height=4)

# N
N.p <- ggplot(filter(P.df, grepl("CA", SDM))) + 
  geom_tile(aes(lon, lat, fill=lN_err_Z)) + 
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank(),
        panel.spacing.x=unit(c(rep(-.1, 7)), "cm"),
        panel.spacing.y=unit(c(-0.1, 0.1, -0.1), "cm"),
        legend.position="bottom",
        legend.key.height=unit(0.25, "cm")) + 
  scale_fill_gradient2("N: Predicted - True (z-transformed log scale)", midpoint=0) + 
  facet_nested(SDM + sp ~ issue, scales="free_x", 
               labeller=labeller(issue=label_wrap_gen(11), SDM=label_parsed)) +
  labs(title="Error in abundance predictions")
ggsave("figs/map_lN_err_Z.jpg", N.p, width=11, height=7)




p <- ggplot() + 
  geom_tile(data=filter(P.df, SDM %in% c("CA[p]", "CA[i]")), 
            aes(lon, lat, fill=Surv.S.f - Surv.S)) +
  scale_fill_gradient2("N: Predicted - True", midpoint=0, na.value="black",
                       breaks=c(-5000,0,5000), limits=c(-5000,5000),
                       low="#40004b", high="#00441b", mid="white",
                       guide=guide_colourbar(title.position="top", 
                                             title.hjust=0.5,
                                             barheight=unit(0.3, "cm"),
                                             barwidth=unit(5, "cm"))) +
  new_scale("fill") +
  geom_tile(data=filter(P.df, SDM=="IPM"),
            aes(lon, lat, fill=lambda.f - lambda)) +
  scale_fill_gradient2(expression(lambda:Predicted-True), midpoint=0, 
                       breaks=c(-3,0,3), limits=c(-3,3), na.value="black",
                       low="#67001f", high="#053061", mid="white",
                       guide=guide_colourbar(title.position="top", 
                                             title.hjust=0.5,
                                             barheight=unit(0.3, "cm"),
                                             barwidth=unit(5, "cm"))) +
  facet_nested(SDM + sp ~ issue, scales="free_x", 
               labeller=labeller(issue=label_wrap_gen(11), SDM=label_parsed)) +
  theme(panel.grid=element_blank(), axis.text=element_blank(), 
        axis.ticks=element_blank(), axis.title=element_blank(),
        panel.spacing.x=unit(0, "mm"),
        panel.spacing.y=unit(c(0, 0.5, 0, 0.5, 0), "mm"),
        legend.position="bottom", 
        legend.box.margin=margin(0, 0, 0, 0, "cm"),
        panel.background=element_rect(fill="gray35"))
ggsave("figs/map_N_Lambda.jpg", p, width=7.25, height=6, dpi=500)





dunn.list <- filter(TSS.df, issue == "Ideal" & sp == "shrub" &
                      Boundary == "True~range:~lambda > 1") %>%
  list(CAp=filter(., SDM=="CA[p]")$TSS,
       CAp=filter(., SDM=="CA[i]")$TSS,
       CAp=filter(., SDM=="IPM")$TSS,
       CAp=filter(., SDM=="MaxEnt")$TSS)
dunn.test(dunn.list[-1])
# For N > 0, 






################################################################################
## SCRAPPED FIGURE DRAFTS
########
# Scrapped figure: TSS caterpillar plot
ggplot(TSS.ci, aes(x=SDM, y=med, colour=issue, group=issue)) + 
  geom_hline(yintercept=seq(0.4, 1, 0.1), colour="gray90", linetype=2, size=0.2) +
  geom_errorbar(aes(ymin=q025, ymax=q975), size=0.3, width=0.7, 
                position=position_dodge(width=1)) + 
  geom_linerange(aes(ymin=q25, ymax=q75), size=1.2, 
                 position=position_dodge(width=1)) + 
  geom_point(shape="-", size=1.75, colour="black", 
             position=position_dodge(width=1)) +
  facet_nested(Boundary ~ sp + SDM, scales="free_x", 
               labeller=labeller(Boundary=label_parsed, SDM=label_parsed)) +
  scale_colour_brewer("Scenario", type="qual", palette=2) +
  theme(panel.grid=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.key.size=unit(0.3, 'cm'), 
        panel.spacing.x=unit(c(-.1,-.1,-.1,0.1,-0.1,-0.1,-0.1), "cm"), 
        panel.spacing.y=unit(0.1, "cm"),
        panel.border=element_rect(colour="gray40")) + 
  labs(title="TSS across 100 sampled datasets", x="", y="TSS")
ggsave("figs/TSS_issueCaterpillar_Sp-Boundary.jpg", width=7.5, height=3.75)


# Scrapped figure: TSS CI widths facets = SDM x Boundary
ggplot(TSS.ci, aes(x=issue, y=q975-q025, colour=issue, shape=sp)) +
  geom_hline(yintercept=seq(0, 0.6, 0.1), colour="gray90", linetype=2, size=0.2) +
  geom_point() + 
  facet_grid(Boundary ~ SDM, scales="free_x", 
             labeller=labeller(Boundary=label_parsed, SDM=label_parsed)) +
  scale_colour_brewer("Scenario", type="qual", palette=2) +
  scale_shape_manual("Species", values=c(19,1)) +
  theme(panel.grid=element_blank(), 
        axis.text.x=element_text(angle=270, vjust=0.5, hjust=0),
        legend.key.size=unit(0.3, 'cm'), 
        panel.border=element_rect(colour="gray40")) + 
  labs(title="TSS variability across 100 sampled datasets", 
       x="Scenario", y="TSS 95% CI width")
ggsave("figs/TSS_CIwidths_Boundary-SDM.jpg", width=6.5, height=4.25)


# Scrapped figure: TSS CI widths facets = Scenarios x Boundary
ggplot(TSS.ci, aes(x=SDM, y=q975-q025, colour=SDM, shape=sp)) +
  geom_hline(yintercept=seq(0, 0.6, 0.1), colour="gray90", linetype=2, size=0.2) +
  geom_point() + 
  facet_grid(Boundary ~ issue, scales="free_x", 
             labeller=labeller(Boundary=label_parsed, 
                               issue=label_wrap_gen(width=11))) +
  scale_colour_manual("SDM\nMethod", values=SDM_col, 
                      labels=c(expression(CA[p], CA[i], IPM, MaxEnt))) +
  scale_shape_manual("Species", values=c(19,1)) +
  theme(panel.grid=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.text.align=0,
        legend.key.size=unit(0.3, 'cm'), 
        panel.spacing.x=unit(-0.1, "cm"), 
        panel.border=element_rect(colour="gray40")) + 
  labs(title="TSS variability across 100 sampled datasets", 
       x="", y="TSS 95% CI width")
ggsave("figs/TSS_CIwidths_Boundary-Issue.jpg", width=9, height=3.25)


# Scrapped figure: Proportional disagreement facets = Scenarios x Species
ggplot(P.sum, aes(x=SDM, y=propDisputedCells, fill=SDM)) + 
  geom_hline(yintercept=0, colour="gray30") + 
  geom_bar(stat="identity", position="dodge", colour="gray30") + 
  scale_fill_manual("SDM\nMethod", values=SDM_col, 
                    labels=c(expression(CA[p], CA[i], IPM, MaxEnt))) +
  facet_grid(sp~issue, labeller=labeller(issue=label_wrap_gen(width=11))) + 
  ylim(0, 0.6) + 
  labs(x="SDM", y="Proportion of cells with\ndisagreement among samples") +
  theme(panel.grid=element_blank(), 
        axis.text.x=element_text(angle=270, vjust=0.5, hjust=0),
        legend.text.align=0,
        legend.key.size=unit(0.3, 'cm'))
ggsave("figs/bar_sampleDisagree_Issue-Sp.jpg", width=9, height=4)


# Scrapped figure: Proportional disagreement facets = Scenario
ggplot(P.sum, aes(x=sp, y=propDisputedCells, fill=SDM)) + 
  geom_hline(yintercept=0, colour="gray30") + 
  geom_bar(stat="identity", position="dodge", colour="gray30") + 
  scale_fill_manual("SDM\nMethod", values=SDM_col, 
                    labels=c(expression(CA[p], CA[i], IPM, MaxEnt))) +
  facet_grid(.~issue, labeller=labeller(issue=label_wrap_gen(width=11))) + 
  ylim(0, NA) + 
  labs(x="Species", y="Proportion of cells with\ndisagreement among samples") +
  theme(panel.grid=element_blank(), 
        axis.text.x=element_text(angle=270, vjust=0.5, hjust=0),
        legend.text.align=0,
        legend.key.size=unit(0.3, 'cm'))
ggsave("figs/bar_sampleDisagree_Issue.jpg", width=10, height=2.5)


# Scrapped figure: Proportional disagreement facets = SDM x Species
ggplot(P.sum, aes(x=issue, y=propDisputedCells, fill=issue)) + 
  geom_hline(yintercept=0, colour="gray30") + 
  geom_bar(stat="identity", position="dodge", colour="gray30") + 
  scale_fill_brewer("Scenario", type="qual", palette=2) +
  facet_grid(sp~SDM, labeller=labeller(SDM=label_parsed)) + 
  ylim(0, NA) + 
  labs(x="Scenario", y="Proportion of cells with\ndisagreement among samples") +
  theme(panel.grid=element_blank(), 
        axis.text.x=element_text(angle=270, vjust=0.5, hjust=0),
        legend.key.size=unit(0.3, 'cm'))
ggsave("figs/bar_sampleDisagree_SDM-Sp.jpg", width=9, height=4)


# Scrapped figure: Proportional disagreement facets = SDM
ggplot(P.sum, aes(x=sp, y=propDisputedCells, fill=issue)) + 
  geom_hline(yintercept=0, colour="gray30") + 
  geom_bar(stat="identity", position="dodge", colour="gray30") + 
  scale_fill_brewer("Scenario", type="qual", palette=2) +
  facet_grid(.~SDM, labeller=labeller(SDM=label_parsed)) + 
  ylim(0, NA) + 
  labs(x="Species", y="Proportion of cells with\ndisagreement among samples") +
  theme(panel.grid=element_blank(), 
        axis.text.x=element_text(angle=270, vjust=0.5, hjust=0),
        legend.key.size=unit(0.3, 'cm'))
ggsave("figs/bar_sampleDisagree_SDM.jpg", width=10, height=2.5)


# Scrapped figure: Proportion of cells with S in CIs facets = Scenario x Species
ggplot(P.sum, aes(x=SDM, y=propSinCI, fill=SDM)) + 
  geom_hline(yintercept=seq(0, 1, 0.25), colour="gray90", linetype=2, size=0.2) +
  geom_hline(yintercept=0, colour="gray30") + 
  geom_bar(stat="identity", position="dodge", colour="gray30") + 
  scale_fill_manual("SDM\nMethod", values=SDM_col, 
                    labels=c(expression(CA[p], CA[i], IPM, MaxEnt))) +
  facet_grid(sp~issue, labeller=labeller(issue=label_wrap_gen(width=11))) + 
  ylim(0, 1) + 
  labs(x="SDM", y="Proportion of cells with\nS in 95% CIs") +
  theme(panel.grid=element_blank(), 
        axis.text.x=element_text(angle=270, vjust=0.5, hjust=0),
        legend.text.align=0,
        legend.key.size=unit(0.3, 'cm'))


# Scrapped figure: Proportion of cells with S in CIs facets = SDM
ggplot(P.sum, aes(x=sp, y=propSinCI, fill=issue)) + 
  geom_hline(yintercept=0, colour="gray30") + 
  geom_bar(stat="identity", position="dodge", colour="gray30") + 
  scale_fill_brewer("Scenario", type="qual", palette=2) +
  facet_grid(.~SDM, labeller=labeller(SDM=label_parsed)) + 
  ylim(0, NA) + 
  labs(x="Species", y="Proportion of cells with\nS in 95% CIs") +
  theme(panel.grid=element_blank(), 
        axis.text.x=element_text(angle=270, vjust=0.5, hjust=0),
        legend.key.size=unit(0.3, 'cm'))


# Scrapped figure: Mean 95% CI width facets = Scenario x Species
ggplot(P.sum, aes(x=SDM, y=mean95CI, fill=SDM)) + 
  geom_bar(stat="identity", position="dodge", colour="gray30") + 
  scale_fill_manual("SDM\nMethod", values=SDM_col, 
                    labels=c(expression(CA[p], CA[i], IPM, MaxEnt))) +
  facet_grid(sp~issue, labeller=labeller(issue=label_wrap_gen(width=11)), 
             scales="free_y") + 
  labs(x="SDM", y="Mean 95% CI width") +
  theme(panel.grid=element_blank(), 
        axis.text.x=element_text(angle=270, vjust=0.5, hjust=0),
        legend.text.align=0,
        legend.key.size=unit(0.3, 'cm'))

# Scrapped figure: Mean 95% CI width facets = SDM
ggplot(P.sum, aes(x=sp, y=mean95CI, fill=issue)) + 
  geom_bar(stat="identity", position="dodge", colour="gray30") + 
  scale_fill_brewer("Scenario", type="qual", palette=2) +
  facet_grid(SDM~., labeller=labeller(SDM=label_parsed), scales="free") + 
  labs(x="Species", y="Mean 95% CI width") +
  theme(panel.grid=element_blank(), 
        axis.text.x=element_text(angle=270, vjust=0.5, hjust=0),
        legend.key.size=unit(0.3, 'cm'))









