# 6: Figures for manuscript
# Comparison of SDM approaches
# Tim Szewczyk

# This script plots the output from `3_fitModels.R` and `4_analysis.R`.

########
## Setup
########

# load workspace
pkgs <- c("tidyverse", "magrittr", "stringr", "here", "viridis", "grid", "gtable")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "fn", full.names=T), source)
theme_set(theme_bw())
sp.names <- c("shrub", "biennial")
SDM_col <- c(MaxEnt="#fd8d3c", IPM="#08306b", 
             'CA[i]'="#2171b5", 'CA[p]'="#6baed6")
model.issues <- c("No Seedbank", "Under dispersal", "Over dispersal")

# load TSS, etc
TSS.df <- rbind(read.csv("out/sp1/out_TSS.csv"),
                read.csv("out/sp2/out_TSS.csv")) %>% 
  rbind(., expand.grid(TSS=NA, sensitivity=NA, specificity=NA,
                       Boundary=unique(.$Boundary),
                       SDM="MxE", issue=c("noSB", "underDisp", "overDisp"),
                       sp=unique(.$sp))) 
TSS.df$issue <- factor(TSS.df$issue, 
                       levels=c("none", "noise", "sampBias", "nonEq", 
                                "wrongCov", "noSB", "underDisp", "overDisp"),
                       labels=c("None", "Measurement error", "Sampling bias", 
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
            se=sd/10,
            q025=quantile(TSS, probs=0.025, na.rm=T), 
            q25=quantile(TSS, probs=0.25, na.rm=T), 
            q75=quantile(TSS, probs=0.75, na.rm=T), 
            q975=quantile(TSS, probs=0.975, na.rm=T)) %>%
  arrange(Boundary, issue, sp, desc(med)) %>%
  ungroup %>% group_by(Boundary, issue, sp) %>%
  mutate(rank_SDM=row_number())  %>%
  full_join(., filter(., issue=="None"), by=c("sp", "SDM", "Boundary"), 
            suffix=c("", ".none"))
TSS.ci$rank_SDM[TSS.ci$SDM=="MaxEnt" & TSS.ci$issue %in% model.issues] <- NA
SDM.ranks <- TSS.ci %>% ungroup %>% 
  filter(!issue %in% model.issues) %>% 
  group_by(SDM, Boundary, sp) %>% 
  summarise(min_rank=min(rank_SDM), mean_rank=mean(rank_SDM), 
            med_rank=median(rank_SDM), max_rank=max(rank_SDM)) %>%
  arrange(Boundary, mean_rank)
rank.bar <- expand.grid(Boundary=unique(TSS.ci$Boundary), 
                        SDM=unique(TSS.ci$SDM), 
                        sp=unique(TSS.ci$sp), 
                        rank_SDM=unique(TSS.ci$rank_SDM)) %>%
  left_join(., TSS.ci %>% 
              filter(!issue %in% model.issues) %>%
              ungroup %>% group_by(Boundary, SDM, rank_SDM, sp) %>%
              summarise(ct=n()), 
            by=c("Boundary", "SDM", "rank_SDM", "sp"))
rank.bar$ct[is.na(rank.bar$ct)] <- 0



# ISEM 2019 specifications
dir.isem <- "~/Documents/conferences/201910_ISEM/figs/"
iss.isem <- c("None", "Non equilibrium", "Incorrect covariates", 
              "No Seedbank", "Under dispersal", "Over dispersal")
fonts.isem <- theme(plot.title=element_text(size=20),
                    strip.text=element_text(size=14),
                    axis.title=element_text(size=14),
                    axis.text=element_text(size=12),
                    legend.text=element_text(size=12))



####----
## Ranks
####----

## Mean rank across all scenarios
# blank
p <- ggplot(SDM.ranks, aes(x=sp, y=mean_rank, colour=SDM)) + 
  geom_hline(yintercept=seq(1, 4, 1), colour="gray90", linetype=2, size=0.2) +
  geom_point(size=4, colour=NA) + 
  scale_colour_manual("SDM\nMethod", values=SDM_col, 
                      labels=c(expression(CA[p], CA[i], IPM, MaxEnt))) +
  facet_grid(Boundary~., labeller=labeller(Boundary=label_parsed)) +
  labs(x="", y="Mean rank across scenarios") + ylim(1, 4.25) +
  theme(panel.grid=element_blank(), 
        legend.position="none", 
        strip.text.y=element_text(hjust=0)) +
  fonts.isem
ggsave(paste0(dir.isem, "Ranks_mn_pt_0.jpg"), p, width=3, height=5)

# points
p <- ggplot(SDM.ranks, aes(x=sp, y=mean_rank, colour=SDM)) + 
  geom_hline(yintercept=seq(1, 4, 1), colour="gray90", linetype=2, size=0.2) +
  geom_point(size=4, position=position_jitter(width=0.02, height=0), alpha=0.8) + 
  scale_colour_manual("SDM\nMethod", values=SDM_col, 
                      labels=c(expression(CA[p], CA[i], IPM, MaxEnt))) +
  facet_grid(Boundary~., labeller=labeller(Boundary=label_parsed)) +
  labs(x="", y="Mean rank across scenarios") + ylim(1, 4.25) +
  theme(panel.grid=element_blank(), 
        legend.position="none", 
        strip.text.y=element_text(hjust=0)) +
  fonts.isem
ggsave(paste0(dir.isem, "Ranks_mn_pt.jpg"), p, width=3, height=5)


## Rank for each scenario
# lines by scenario
ggplot(TSS.ci, aes(x=issue, y=rank_SDM, colour=SDM, group=SDM)) +
  geom_hline(yintercept=seq(1, 4, 1), colour="gray90", linetype=2, size=0.2) +
  geom_line(size=2) + 
  scale_colour_manual("SDM\nMethod", values=SDM_col, 
                      labels=c(expression(CA[p], CA[i], IPM, MaxEnt))) +
  facet_grid(Boundary~sp, labeller=labeller(Boundary=label_parsed)) +
  labs(x="", y="Rank") + ylim(1, 4.25) +
  theme(panel.grid=element_blank(), 
        axis.text.x=element_text(angle=270, vjust=0.5, hjust=0),
        legend.position="none", 
        strip.text.y=element_text(hjust=0)) +
  fonts.isem
  
# barplot
p <- ggplot(rank.bar, aes(x=rank_SDM, y=ct, fill=SDM)) + 
  geom_hline(yintercept=seq(1, 4, 1), colour="gray90", linetype=2, size=0.2) +
  geom_bar(stat="identity", position="dodge", colour="gray30") +
  scale_fill_manual("SDM\nMethod", values=SDM_col, 
                      labels=c(expression(CA[p], CA[i], IPM, MaxEnt))) +
  facet_grid(Boundary~sp, labeller=labeller(Boundary=label_parsed)) +
  labs(x="Rank", y="Count across scenarios") +
  theme(panel.grid=element_blank(), 
        legend.position="none", 
        strip.text.y=element_text(hjust=0)) +
  fonts.isem
ggsave(paste0(dir.isem, "Ranks_bar.jpg"), p, width=7, height=5)




####----
## Mean TSS, etc
####----

# Mean TSS ± 2*SE
mn_jitt <- position_dodge(width=0.4)
for(i in 1:length(iss.isem)) {
  p <- ggplot(filter(TSS.ci, issue %in% iss.isem), 
         aes(x=issue, y=mn, colour=SDM, group=SDM)) + 
    geom_point(colour=NA) +
    geom_hline(yintercept=seq(0.7, 1, 0.1), colour="gray80", 
               linetype=2, size=0.5) +
    geom_point(data=filter(TSS.ci, issue %in% iss.isem[1:i]), shape=1, 
               size=4, position=mn_jitt) +
    geom_point(data=filter(TSS.ci, issue %in% iss.isem[1:i]), alpha=0.7,#shape=95, 
               size=4, position=mn_jitt) +
    geom_errorbar(data=filter(TSS.ci, issue %in% iss.isem[1:i]),
                  aes(ymin=mn-2*se, ymax=mn+2*se), size=0.5, width=0.25,
                  position=mn_jitt) + 
    facet_grid(Boundary~sp, labeller=labeller(Boundary=label_parsed)) +
    scale_x_discrete(labels=scales::wrap_format(10)) +
    scale_colour_manual("SDM\nMethod", values=SDM_col, 
                        labels=c(expression(CA[p], CA[i], IPM, MaxEnt))) +
    theme(panel.grid=element_blank(), 
          strip.text.y=element_text(hjust=0),
          axis.text.x=element_text(angle=270, vjust=0.5, hjust=0),
          legend.position="none") + 
    fonts.isem +
    labs(x="", y="TSS (mean ± 2 SE)")
  ggsave(paste0(dir.isem, "TSS_mn+SE_", i, ".jpg"), p, width=8, height=6)
}

# Mean TSS ± SD
for(i in 1:length(iss.isem)) {
  p <- ggplot(filter(TSS.ci, issue %in% iss.isem), 
              aes(x=issue, y=mn, colour=SDM, group=SDM)) + 
    geom_point(colour=NA) +
    geom_hline(yintercept=seq(0.7, 1, 0.1), colour="gray80", 
               linetype=2, size=0.5) +
    geom_point(data=filter(TSS.ci, issue %in% iss.isem[1:i]), shape=1, 
               size=4, position=mn_jitt) +
    geom_point(data=filter(TSS.ci, issue %in% iss.isem[1:i]), alpha=0.7,#shape=95, 
               size=4, position=mn_jitt) +
    geom_errorbar(data=filter(TSS.ci, issue %in% iss.isem[1:i]),
                  aes(ymin=mn-sd, ymax=mn+sd), size=0.5, width=0.25,
                  position=mn_jitt) + 
    facet_grid(Boundary~sp, labeller=labeller(Boundary=label_parsed)) +
    scale_x_discrete(labels=scales::wrap_format(10)) +
    scale_colour_manual("SDM\nMethod", values=SDM_col, 
                        labels=c(expression(CA[p], CA[i], IPM, MaxEnt))) +
    theme(panel.grid=element_blank(), 
          strip.text.y=element_text(hjust=0),
          axis.text.x=element_text(angle=270, vjust=0.5, hjust=0),
          legend.position="none") + 
    fonts.isem +
    labs(x="", y="TSS (mean ± 2 SE)")
  ggsave(paste0(dir.isem, "TSS_mn+SD_", i, ".jpg"), p, width=8, height=6)
}




####----
## TSS vs. none
####----

## Scenario effects on median TSS
# blank
p <- ggplot(filter(TSS.ci, issue != "None" & issue %in% iss.isem), 
            aes(x=issue, y=med-med.none, colour=SDM, group=SDM, shape=sp)) + 
  geom_hline(yintercept=seq(-0.15, 0.05, 0.05), colour="gray90", 
             linetype=2, size=0.2) + 
  geom_hline(yintercept=0, colour="gray80", linetype=2, size=1) +
  geom_point(position=position_dodge(width=1), size=4, colour=NA) +
  facet_grid(Boundary~issue, scales="free_x", 
             labeller=labeller(Boundary=label_parsed, issue=label_wrap_gen(11))) +
  scale_colour_manual("SDM\nMethod", values=SDM_col, 
                      labels=c(expression(CA[p], CA[i], IPM, MaxEnt))) +
  scale_shape_manual("Species", values=c(19,17)) +
  scale_y_continuous(breaks=c(-0.1, 0, 0.1), limits=c(-0.15,0.1)) +
  theme(panel.grid=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text.y=element_text(hjust=0),
        legend.position="none",
        panel.spacing.x=unit(c(rep(-.1, 4)), "cm"),
        panel.spacing.y=unit(0.5, "cm"),
        panel.border=element_rect(colour="gray40")) + 
  fonts.isem +
  labs(x="", y=expression(TSS[i]-TSS['None']))
ggsave(paste0(dir.isem, "TSS_md_v_None_0.jpg"), p, width=7.5, height=5.25)

# points
for(i in 2:length(iss.isem)) {
  p <- ggplot(filter(TSS.ci, issue != "None" & issue %in% iss.isem), 
              aes(x=issue, y=med-med.none, colour=SDM, group=SDM, shape=sp)) + 
    geom_point(position=position_dodge(width=1), size=4, colour=NA) +
    geom_hline(yintercept=seq(-0.15, 0.05, 0.05), colour="gray90", 
               linetype=2, size=0.2) + 
    geom_hline(yintercept=0, colour="gray80", linetype=2, size=1) +
    geom_point(data=filter(TSS.ci, issue %in% iss.isem[2:i]),
               position=position_dodge(width=1), size=4, alpha=0.7) +
    facet_grid(Boundary~issue, scales="free_x", 
               labeller=labeller(Boundary=label_parsed, issue=label_wrap_gen(11))) +
    scale_colour_manual("SDM\nMethod", values=SDM_col, 
                        labels=c(expression(CA[p], CA[i], IPM, MaxEnt))) +
    scale_shape_manual("Species", values=c(19,17)) +
    scale_y_continuous(breaks=c(-0.1, 0, 0.1), limits=c(-0.15,0.1)) +
    theme(panel.grid=element_blank(), 
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.text.y=element_text(hjust=0),
          legend.position="none",
          panel.spacing.x=unit(c(rep(-.1, 4)), "cm"),
          panel.spacing.y=unit(0.5, "cm"),
          panel.border=element_rect(colour="gray40")) + 
    fonts.isem +
    labs(x="", y=expression(TSS[i]-TSS['None']))
  ggsave(paste0(dir.isem, "TSS_md_v_None_", i-1, ".jpg"), p, width=7.5, height=5.25)
}


ggplot(filter(TSS.ci, issue != "None" & issue %in% iss.isem), 
       aes(x=issue, y=(med-med.none)/med.none, colour=SDM, group=SDM)) + 
  # geom_hline(yintercept=seq(-0.15, 0.05, 0.05), colour="gray90", 
             # linetype=2, size=0.2) + 
  geom_hline(yintercept=0, colour="gray80", linetype=2, size=1) +
  geom_point(position=position_jitter(width=0.2), size=4, alpha=0.7) +
  facet_grid(Boundary~sp, 
             labeller=labeller(Boundary=label_parsed, issue=label_wrap_gen(11))) +
  scale_colour_manual("SDM\nMethod", values=SDM_col, 
                      labels=c(expression(CA[p], CA[i], IPM, MaxEnt))) +
  scale_y_continuous(breaks=c(-0.1, 0, 0.1), limits=c(-0.15,0.1)) +
  theme(panel.grid=element_blank(), 
        strip.text.y=element_text(hjust=0),
        legend.position="none",
        axis.text.x=element_text(angle=270, vjust=0.5, hjust=0),
        # panel.spacing.x=unit(c(rep(-.1, 4)), "cm"),
        # panel.spacing.y=unit(0.5, "cm"),
        panel.border=element_rect(colour="gray40")) + 
  fonts.isem +
  labs(x="", y=expression(TSS[i]-TSS['None']))




