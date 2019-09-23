# 5: Plot output
# Comparison of SDM approaches
# Tim Szewczyk

# This script plots the output from `3_fitModels.R` and `4_analysis.R`.

########
## Setup
########
# file specifications
sp <- "sp2"
plots <- TRUE
sp.names <- c("barberry", "garlic mustard")[ifelse(sp=="sp1", 1, 2)]

# load workspace
pkgs <- c("tidyverse", "magrittr", "stringr", "here", "viridis")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "fn", full.names=T), source)
theme_set(theme_bw())
lam.df <- readRDS(here("vs", sp, "lam_df.rds"))
out_P <- read.csv(here("out", sp, "out_P.csv"))
out_TSS <- read.csv(here("out", sp, "out_TSS.csv"))
out_LL <- read.csv(here("out", sp, "out_LL.csv"))
out_P$issue <- factor(out_P$issue, 
                    levels=c("none", "noise", "sampBias", "nonEq",
                             "noSB", "underDisp", "overDisp", "wrongCov"),
                    labels=c("None", "Measurement error", "Sampling bias", 
                             "Non-equilibrium", "No Seedbank", "Under dispersal",
                             "Over dispersal", "Incorrect covariates"))
out_TSS$issue <- factor(out_TSS$issue, 
                      levels=c("none", "noise", "sampBias", "nonEq",
                               "noSB", "underDisp", "overDisp", "wrongCov"),
                      labels=c("None", "Measurement error", "Sampling bias", 
                               "Non-equilibrium", "No Seedbank", "Under dispersal",
                               "Over dispersal", "Incorrect covariates"))
out_LL$issue <- factor(out_LL$issue, 
                        levels=c("none", "noise", "sampBias", "nonEq",
                                 "noSB", "underDisp", "overDisp", "wrongCov"),
                        labels=c("None", "Measurement error", "Sampling bias", 
                                 "Non-equilibrium", "No Seedbank", "Under dispersal",
                                 "Over dispersal", "Incorrect covariates"))
mn_TSS <- out_TSS %>% group_by(Boundary, SDM, issue) %>% 
  summarise(TSS=mean(TSS), 
            sensitivity=mean(sensitivity), 
            specificity=mean(specificity))
SDM_col <- c(MxE="#7b3294", IPM="#00441b", CAi="#045a8d", CAd="#74a9cf")


# TSS, Sensitivity, Specificity: Density Plots
ggplot(out_TSS, aes(x=TSS, fill=SDM, colour=SDM)) + 
  geom_vline(xintercept=c(0,1), colour="grey90") +
  geom_vline(xintercept=0.5, colour="grey90", linetype=3) +
  geom_density(alpha=0.75) + 
  geom_rug(data=mn_TSS, size=1) +
  facet_grid(issue~Boundary, scales="free_y", 
             labeller=labeller(issue=label_wrap_gen(10), Boundary=label_parsed)) +
  scale_fill_manual("SDM\nMethod", values=SDM_col) + 
  scale_colour_manual(values=SDM_col, guide=F) +
  scale_x_continuous("TSS", limits=c(0,1), breaks=c(0, 0.5, 1)) + ylab("Density") +
  theme(panel.grid=element_blank()) + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/eval/", sp, "_TSS.jpg"), width=6, height=8)

ggplot(out_TSS, aes(x=sensitivity, fill=SDM, colour=SDM)) + 
  geom_vline(xintercept=c(0,1), colour="grey90") +
  geom_vline(xintercept=0.5, colour="grey90", linetype=3) +
  geom_density(alpha=0.75) + 
  geom_rug(data=mn_TSS, size=1) +
  facet_grid(issue~Boundary, scales="free_y", 
             labeller=labeller(issue=label_wrap_gen(10), Boundary=label_parsed)) +
  scale_fill_manual("SDM\nMethod", values=SDM_col) + 
  scale_colour_manual(values=SDM_col, guide=F) +
  scale_x_continuous("Sensitivity", c(0,0.5,1)) + ylab("Density") +
  theme(panel.grid=element_blank()) + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/eval/", sp, "_Sensitivity.jpg"), width=6, height=8)

ggplot(out_TSS, aes(x=specificity, fill=SDM, colour=SDM)) + 
  geom_vline(xintercept=c(0,1), colour="grey90") +
  geom_vline(xintercept=0.5, colour="grey90", linetype=3) +
  geom_density(alpha=0.75) + 
  geom_rug(data=mn_TSS, size=1) +
  facet_grid(issue~Boundary, scales="free_y", 
             labeller=labeller(issue=label_wrap_gen(10), Boundary=label_parsed)) +
  scale_fill_manual("SDM\nMethod", values=SDM_col) + 
  scale_colour_manual(values=SDM_col, guide=F) +
  scale_x_continuous("Specificity", c(0,0.5,1)) + ylab("Density") +
  theme(panel.grid=element_blank()) + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/eval/", sp, "_Specificity.jpg"), width=6, height=8)

out_TSS %>% group_by(SDM, issue, Boundary) %>% 
  summarise(med=median(TSS, na.rm=T), 
            q025=quantile(TSS, probs=0.025, na.rm=T), 
            q975=quantile(TSS, probs=0.975, na.rm=T)) %>% 
  ggplot(aes(x=SDM, colour=SDM)) + 
  geom_pointrange(aes(ymin=q025, y=med, ymax=q975), fatten=0.75) + 
  facet_grid(Boundary~issue,
             labeller=labeller(issue=label_wrap_gen(10), Boundary=label_parsed)) +
  scale_fill_manual("SDM\nMethod", values=SDM_col) + 
  scale_colour_manual(values=SDM_col, guide=F) +
  theme(panel.grid=element_blank(), 
        axis.text.x=element_text(angle=270, vjust=0.5)) + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/eval/", sp, "_TSS_95CI_byIssue.jpg"), width=9, height=4)

out_TSS %>% group_by(SDM, issue, Boundary) %>% 
  summarise(med=median(TSS, na.rm=T), 
            q025=quantile(TSS, probs=0.025, na.rm=T), 
            q975=quantile(TSS, probs=0.975, na.rm=T)) %>% 
  ggplot(aes(x=issue, colour=issue)) + 
  geom_pointrange(aes(ymin=q025, y=med, ymax=q975), fatten=0.75) + 
  facet_grid(Boundary~SDM,
             labeller=labeller(Boundary=label_parsed)) +
  scale_fill_brewer(type="qual", palette=2) +
  scale_colour_brewer(type="qual", palette=2) +
  theme(panel.grid=element_blank(), 
        axis.text.x=element_text(angle=270, vjust=0.5, hjust=0)) + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/eval/", sp, "_TSS_95CI_bySDM.jpg"), width=9, height=4)



# TSS, Sensitivity, Specificity: Boxplots
ggplot(out_TSS, aes(y=TSS, x=SDM, fill=SDM, colour=SDM)) + 
  geom_hline(yintercept=c(0,1), colour="grey90") +
  geom_hline(yintercept=0.5, colour="grey90", linetype=3) +
  geom_boxplot(alpha=0.75, outlier.size=0.1) + 
  facet_grid(Boundary~issue, scales="free_y", 
             labeller=labeller(issue=label_wrap_gen(10), Boundary=label_parsed)) +
  scale_fill_manual("SDM\nMethod", values=SDM_col) + 
  scale_colour_manual(values=SDM_col, guide=F) +
  scale_y_continuous("TSS", limits=c(0,1), breaks=c(0, 0.5, 1)) + 
  theme(panel.grid=element_blank(), 
        axis.text.x=element_text(angle=270, vjust=0.5)) + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/eval/", sp, "_TSS_box_byIssue.jpg"), width=9, height=4)

ggplot(out_TSS, aes(y=TSS, x=issue, fill=issue, colour=issue)) + 
  geom_hline(yintercept=c(0,1), colour="grey90") +
  geom_hline(yintercept=0.5, colour="grey90", linetype=3) +
  geom_boxplot(alpha=0.75, outlier.size=0.1) + 
  facet_grid(Boundary~SDM, scales="free_y", 
             labeller=labeller(Boundary=label_parsed)) +
  scale_fill_brewer(type="qual", palette=2) +
  scale_colour_brewer(type="qual", palette=2) +
  scale_y_continuous("TSS", limits=c(0,1), breaks=c(0, 0.5, 1)) + 
  theme(panel.grid=element_blank(), 
        axis.text.x=element_text(angle=270, vjust=0.5, hjust=0)) + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/eval/", sp, "_TSS_box_bySDM.jpg"), width=9, height=4)



# Log Likelihood
ggplot(out_LL, aes(x=LogLik, y=fct_rev(issue), colour=SDM)) + 
  facet_grid(Boundary~.) +
  geom_point(alpha=0.75, size=3) + 
  scale_colour_manual(values=SDM_col) +
  labs(x="Log likelihood: S ~ Binom(prP)", y="") + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/eval/", sp, "_LogLik.jpg"), width=5, height=4)



# S
if(plots) jpeg(paste0("eda/S/", sp, "_scatter.jpg"), 
               width=8, height=8, units="in", res=300)
par(mfrow=c(3,3))
plot(lam.df$bio10_1, log(lam.df$lambda), col=rgb(0,0,0,0.25))
plot(lam.df$bio10_12, log(lam.df$lambda), col=rgb(0,0,0,0.25))
plot(lam.df$bio10_6, log(lam.df$lambda), col=rgb(0,0,0,0.25))
plot(lam.df$bio10_prMay, log(lam.df$lambda), col=rgb(0,0,0,0.25))
hist(log(lam.df$lambda), main="")
plot(lam.df$bio10_1, log(lam.df$Surv.S), col=rgb(0,0,0,0.25))
plot(lam.df$bio10_12, log(lam.df$Surv.S), col=rgb(0,0,0,0.25))
plot(lam.df$bio10_6, log(lam.df$Surv.S), col=rgb(0,0,0,0.25))
plot(lam.df$bio10_prMay, log(lam.df$Surv.S), col=rgb(0,0,0,0.25))
if(plots) dev.off()

ggplot() + 
  geom_tile(data=lam.df, aes(lon, lat, fill=lambda)) +
  geom_tile(data=filter(lam.df, lambda<1), aes(lon, lat), fill="gray", alpha=0.2) +
  scale_fill_viridis(option="B") + 
  ggtitle(sp.names)
if(plots) ggsave(paste0("eda/S/", sp, "_lambda.jpg"), width=9, height=8)

ggplot() + geom_tile(data=lam.df, aes(lon, lat), fill="gray30") +
  geom_tile(data=filter(lam.df, Surv.S>0), aes(lon, lat, fill=Surv.S)) +
  scale_fill_viridis("N (t.max)", option="B") + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/S/", sp, "_N_Eq.jpg"), width=9, height=8)

ggplot() + geom_tile(data=lam.df, aes(lon, lat), fill="gray30") +
  geom_tile(data=filter(lam.df, Surv.S_nonEq>0), aes(lon, lat, fill=Surv.S_nonEq)) +
  scale_fill_viridis("N (t.1/3)", option="B") + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/S/", sp, "_N_nonEq.jpg"), width=9, height=8)

ggplot() + theme_bw() +
  geom_tile(data=lam.df, aes(x=lon, y=lat), fill="gray90") +
  geom_tile(data=filter(lam.df, lambda>1), aes(x=lon, y=lat), 
            fill="darkblue", alpha=0.5) +
  geom_tile(data=filter(lam.df, Surv.S>0), aes(x=lon, y=lat), 
            fill="red", alpha=0.5) +
  theme(axis.text=element_blank()) + labs(x="", y="") +
  ggtitle(paste0(sp, ": 5km x 5km, favorable habitat"))



# P: prP
p <- ggplot(out_P, aes(x=lon, y=lat, fill=prP)) +
  geom_tile() + facet_grid(SDM~issue) + 
  scale_fill_viridis("prob(P)", option="B") +
  theme(axis.text=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank()) + 
  ggtitle(sp.names) + labs(x="", y="")
if(plots) ggsave(paste0("eda/P/", sp, "_prP_map.jpg"), p, width=12, height=7)



# P: lambda, N, etc
p <- ggplot() + geom_tile(data=filter(out_P, SDM!="MxE"),
                          aes(lon, lat), fill="gray30") +
  geom_tile(data=filter(out_P, SDM=="IPM"), 
            aes(lon, lat, fill=lambda.sd.f/lambda.f)) +
  geom_tile(data=filter(out_P, SDM %in% c("CAi", "CAd")), 
            aes(lon, lat, fill=Surv.S.sd.f/Surv.S.f)) +
  facet_grid(SDM~issue) + 
  scale_fill_viridis("CV") +
  theme(axis.text=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank()) + 
  ggtitle(sp.names) + labs(x="", y="")
if(plots) ggsave(paste0("eda/P/", sp, "_map_CV.jpg"), p, width=12, height=6.5)

p <- ggplot() + geom_tile(data=filter(out_P, SDM!="MxE"),
                          aes(lon, lat), fill="gray30") +
  geom_tile(data=filter(out_P, SDM=="IPM"), 
            aes(lon, lat, fill=lambda.sd.f/max(lambda.sd.f, na.rm=T))) +
  geom_tile(data=filter(out_P, SDM %in% c("CAi", "CAd")), 
            aes(lon, lat, fill=Surv.S.sd.f/max(Surv.S.sd.f, na.rm=T))) +
  facet_grid(SDM~issue) + 
  scale_fill_viridis("scaled SD") +
  theme(axis.text=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank()) + 
  ggtitle(sp.names) + labs(x="", y="")
if(plots) ggsave(paste0("eda/P/", sp, "_map_SD.jpg"), p, width=12, height=6.5)

p <- ggplot() + geom_tile(data=filter(out_P, SDM!="MxE"),
                          aes(lon, lat), fill="gray30") +
  geom_tile(data=filter(out_P, SDM=="IPM"),
            aes(lon, lat, fill=(lambda.975.f-lambda.025.f)/lambda.f)) +
  geom_tile(data=filter(out_P, SDM %in% c("CAi", "CAd")),
            aes(lon, lat, fill=(Surv.S.975.f-Surv.S.025.f)/Surv.S.f)) +
  facet_grid(SDM~issue) + 
  scale_fill_viridis("95%CI width / mean") +
  theme(axis.text=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank()) + 
  ggtitle(sp.names) + labs(x="", y="")
if(plots) ggsave(paste0("eda/P/", sp, "_map_95CI_std.jpg"), p, width=12, height=6.5)

p <- ggplot() + geom_tile(data=filter(out_P, SDM!="MxE"),
                          aes(lon, lat), fill="gray30") +
  geom_tile(data=filter(out_P, SDM=="IPM"),
            aes(lon, lat, fill=lambda.975.f-lambda.025.f)) +
  geom_tile(data=filter(out_P, SDM %in% c("CAi", "CAd")),
            aes(lon, lat, fill=Surv.S.975.f-Surv.S.025.f)) +
  facet_grid(SDM~issue) + 
  scale_fill_viridis("95%CI width") +
  theme(axis.text=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank()) + 
  ggtitle(sp.names) + labs(x="", y="")
if(plots) ggsave(paste0("eda/P/", sp, "_map_95CI.jpg"), p, width=12, height=6.5)

p <- ggplot() + geom_tile(data=filter(out_P, SDM=="IPM"), 
                     aes(lon, lat), fill="gray30") +
  geom_tile(data=filter(out_P, SDM=="IPM" & lambda.f>1), 
            aes(lon, lat, fill=lambda.f)) +
  facet_grid(SDM~issue) + 
  scale_fill_viridis("lambda", option="B") +
  theme(axis.text=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank()) + 
  ggtitle(sp.names) + labs(x="", y="")
if(plots) ggsave(paste0("eda/P/", sp, "_map_lambda.jpg"), p, width=12, height=2.5)

p <- ggplot() + geom_tile(data=filter(out_P, grepl("CA", SDM)), 
                     aes(x=lon, y=lat), fill="gray30") +
  geom_tile(data=filter(out_P, grepl("CA", SDM) & Surv.S.f>0.5), 
            aes(lon, lat, fill=Surv.S.f)) +
  facet_grid(SDM~issue) + 
  scale_fill_viridis("N", option="B") +
  theme(axis.text=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank()) + 
  ggtitle(sp.names) + labs(x="", y="")
if(plots) ggsave(paste0("eda/P/", sp, "_map_N.jpg"), p, width=12, height=4)

p <- ggplot() + geom_tile(data=filter(out_P, grepl("CA", SDM)), 
                     aes(x=lon, y=lat), fill="gray30") +
  geom_tile(data=filter(out_P, grepl("CA", SDM) & (Rcr.S.f)>0.5), 
            aes(lon, lat, fill=Rcr.S.f)) +
  facet_grid(SDM~issue) + 
  scale_fill_viridis("Recruits", option="B") +
  theme(axis.text=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank()) + 
  ggtitle(sp.names) + labs(x="", y="")
if(plots) ggsave(paste0("eda/P/", sp, "_map_Rcr.jpg"), p, width=12, height=4)

p <- ggplot() + geom_tile(data=filter(out_P, grepl("CA", SDM)), 
                     aes(x=lon, y=lat), fill="gray30") +
  geom_tile(data=filter(out_P, grepl("CA", SDM) & B.f>0.5), 
            aes(lon, lat, fill=B.f)) +
  facet_grid(SDM~issue) + 
  scale_fill_viridis("B", option="B") +
  theme(axis.text=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank()) + 
  ggtitle(sp.names) + labs(x="", y="")
if(plots) ggsave(paste0("eda/P/", sp, "_map_B.jpg"), p, width=12, height=4)

p <- ggplot() +
  geom_tile(data=filter(out_P, grepl("CA", SDM)), 
            aes(x=lon, y=lat), fill="gray30") +
  geom_tile(data=filter(out_P, grepl("CA", SDM) & D.f>0.5), 
            aes(lon, lat, fill=D.f)) +
  facet_grid(SDM~issue) + 
  scale_fill_viridis("D", option="B") +
  theme(axis.text=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank()) + 
  ggtitle(sp.names) + labs(x="", y="")
if(plots) ggsave(paste0("eda/P/", sp, "_map_D.jpg"), p, width=12, height=4)

p <- ggplot() +
  geom_tile(data=filter(out_P, grepl("CA", SDM)), 
            aes(x=lon, y=lat), fill="gray30") +
  geom_tile(data=filter(out_P, grepl("CA", SDM) & (nSdStay.f+nSdLeave.f)>0.5), 
            aes(lon, lat, fill=nSdStay.f+nSdLeave.f)) +
  facet_grid(SDM~issue) + 
  scale_fill_viridis("SdTot", option="B") +
  theme(axis.text=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank()) + 
  ggtitle(sp.names) + labs(x="", y="")
if(plots) ggsave(paste0("eda/P/", sp, "_map_nSeed.jpg"), p, width=12, height=4)

p <- ggplot() +
  geom_tile(data=filter(out_P, grepl("CA", SDM)), 
            aes(x=lon, y=lat), fill="gray30") +
  geom_tile(data=filter(out_P, grepl("CA", SDM) & (nSdStay.f)>0.5), 
            aes(lon, lat, fill=nSdStay.f)) +
  facet_grid(SDM~issue) + 
  scale_fill_viridis("SdStay", option="B") +
  theme(axis.text=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank()) + 
  ggtitle(sp.names) + labs(x="", y="")
if(plots) ggsave(paste0("eda/P/", sp, "_map_nSeedLocal.jpg"), p, width=12, height=4)

p <- ggplot() +
  geom_tile(data=filter(out_P, grepl("CA", SDM)), 
            aes(x=lon, y=lat), fill="gray30") +
  geom_tile(data=filter(out_P, grepl("CA", SDM) & (nSdStay.f+D.f)>0.5), 
            aes(lon, lat, fill=nSdStay.f+D.f)) +
  facet_grid(SDM~issue) + 
  scale_fill_viridis("Seed\nPressure", option="B") +
  theme(axis.text=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank()) + 
  ggtitle(sp.names) + labs(x="", y="")
if(plots) ggsave(paste0("eda/P/", sp, "_map_propagule.jpg"), p, width=12, height=4)



# Diff: P - S
p <-ggplot(filter(out_P, grepl("CA", SDM))) + 
  geom_tile(aes(lon, lat, fill=Surv.S.f - Surv.S)) + 
  theme(axis.text=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank()) + 
  scale_fill_gradient2("N\nFit-True", midpoint=0) + 
  facet_grid(SDM~issue) + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/diff_S/", sp, "_map_N.jpg"), p, width=12, height=4)

p <- ggplot(filter(out_P, grepl("CA", SDM))) + 
  geom_tile(aes(lon, lat, fill=B.f - B)) + 
  theme(axis.text=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank()) + 
  scale_fill_gradient2("B\nFit-True", midpoint=0) + 
  facet_grid(SDM~issue) + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/diff_S/", sp, "_map_B.jpg"), p, width=12, height=4)

p <- ggplot(filter(out_P, grepl("CA", SDM))) + 
  geom_tile(aes(lon, lat, fill=D.f - D)) + 
  theme(axis.text=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank()) + 
  scale_fill_gradient2("D\nFit-True", midpoint=0) + 
  facet_grid(SDM~issue) + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/diff_S/", sp, "_map_D.jpg"), p, width=12, height=4)

p <- ggplot(filter(out_P, grepl("CA", SDM))) + 
  geom_tile(aes(lon, lat, fill=nSeed.f - nSeed)) + 
  theme(axis.text=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank()) + 
  scale_fill_gradient2("nSeed\nFit-True", midpoint=0) + 
  facet_grid(SDM~issue) + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/diff_S/", sp, "_map_nSeed.jpg"), p, width=12, height=4)

p <- ggplot(filter(out_P, grepl("CA", SDM))) + 
  geom_tile(aes(lon, lat, fill=(nSdStay.f+D.f) - (nSdStay+D))) + 
  theme(axis.text=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank()) + 
  scale_fill_gradient2("Seed\nPressure\nFit-True", midpoint=0) + 
  facet_grid(SDM~issue) + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/diff_S/", sp, "_map_propagule.jpg"), p, width=12, height=4)

p <- ggplot(filter(out_P, SDM=="IPM")) + 
  geom_tile(aes(lon, lat, fill=lambda.f - lambda)) + 
  theme(axis.text=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank()) + 
  scale_fill_gradient2("lambda\nFit-True", midpoint=0) + 
  facet_grid(SDM~issue) + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/diff_S/", sp, "_map_lambda.jpg"), width=12, height=2.5)



# Diff: P - P.none
p <- ggplot(out_P) + 
  geom_tile(aes(lon, lat, fill=prP - prP.none)) + 
  theme(axis.text=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank()) + 
  scale_fill_gradient2("prP\nFit-None", midpoint=0) + 
  facet_grid(SDM~issue) + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/diff_none/", sp, "_map_prP.jpg"), p, width=12, height=7)

p <- ggplot(out_P) + 
  geom_tile(data=filter(out_P, SDM %in% c("CAi", "CAd")),
                        aes(lon, lat, fill=Surv.S.sd.f - Surv.S.sd.f.none)) + 
  geom_tile(data=filter(out_P, SDM=="IPM"),
            aes(lon, lat, fill=lambda.sd.f - lambda.sd.f.none)) +               
  theme(axis.text=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank()) + 
  scale_fill_gradient2("SD\nFit-None", midpoint=0) + 
  facet_grid(SDM~issue) + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/diff_none/", sp, "_map_SD.jpg"), p, width=12, height=7)

p <- ggplot(out_P) + 
  geom_tile(data=filter(out_P, SDM %in% c("CAi", "CAd")),
            aes(lon, lat, fill=(Surv.S.975.f-Surv.S.025.f) - (Surv.S.975.f.none-Surv.S.025.f.none))) + 
  geom_tile(data=filter(out_P, SDM=="IPM"),
            aes(lon, lat, fill=(lambda.975.f-lambda.025.f) - (lambda.975.f.none-lambda.025.f.none))) + 
  theme(axis.text=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank()) + 
  scale_fill_gradient2("95% CI width\nFit-None", midpoint=0) + 
  facet_grid(SDM~issue) + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/diff_none/", sp, "_map_95CI.jpg"), p, width=12, height=7)

p <- ggplot(out_P) + 
  geom_tile(data=filter(out_P, SDM %in% c("CAi", "CAd")),
            aes(lon, lat, fill=(Surv.S.975.f-Surv.S.025.f)/Surv.S.f - (Surv.S.975.f.none-Surv.S.025.f.none)/Surv.S.f.none)) + 
  geom_tile(data=filter(out_P, SDM=="IPM"),
            aes(lon, lat, fill=(lambda.975.f-lambda.025.f)/lambda.f - (lambda.975.f.none-lambda.025.f.none)/lambda.f.none)) + 
  theme(axis.text=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank()) + 
  scale_fill_gradient2("95% CI width\nFit-None", midpoint=0) + 
  facet_grid(SDM~issue) + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/diff_none/", sp, "_map_95CI_std.jpg"), p, width=12, height=7)

p <- ggplot(out_P) + 
  geom_tile(data=filter(out_P, SDM %in% c("CAi", "CAd")),
            aes(lon, lat, fill=Surv.S.sd.f/Surv.S.f - Surv.S.sd.f.none/Surv.S.f.none)) + 
  geom_tile(data=filter(out_P, SDM=="IPM"),
            aes(lon, lat, fill=lambda.sd.f/lambda.f - lambda.sd.f.none/lambda.f.none)) +
  theme(axis.text=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank()) + 
  scale_fill_gradient2("CV\nFit-None", midpoint=0) + 
  facet_grid(SDM~issue) + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/diff_none/", sp, "_map_CV.jpg"), p, width=12, height=7)

p <- ggplot(filter(out_P, grepl("CA", SDM))) + 
  geom_tile(aes(lon, lat, fill=Surv.S.f - Surv.S.f.none)) + 
  theme(axis.text=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank()) + 
  scale_fill_gradient2("N\nFit-None", midpoint=0) + 
  facet_grid(SDM~issue) + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/diff_none/", sp, "_map_N.jpg"), p, width=12, height=4)

p <- ggplot(filter(out_P, grepl("CA", SDM))) + 
  geom_tile(aes(lon, lat, fill=B.f - B.f.none)) + 
  theme(axis.text=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank()) + 
  scale_fill_gradient2("B\nFit-None", midpoint=0) + 
  facet_grid(SDM~issue) + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/diff_none/", sp, "_map_B.jpg"), p, width=12, height=4)

p <- ggplot(filter(out_P, grepl("CA", SDM))) + 
  geom_tile(aes(lon, lat, fill=D.f - D.f.none)) + 
  theme(axis.text=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank()) + 
  scale_fill_gradient2("D\nFit-None", midpoint=0) + 
  facet_grid(SDM~issue) + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/diff_none/", sp, "_map_D.jpg"), p, width=12, height=4)

p <- ggplot(filter(out_P, grepl("CA", SDM))) + 
  geom_tile(aes(lon, lat, fill=nSeed.f - nSeed.f.none)) + 
  theme(axis.text=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank()) + 
  scale_fill_gradient2("nSeed\nFit-None", midpoint=0) + 
  facet_grid(SDM~issue) + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/diff_none/", sp, "_map_nSeed.jpg"), p, width=12, height=4)

p <- ggplot(filter(out_P, grepl("CA", SDM))) + 
  geom_tile(aes(lon, lat, fill=(nSdStay.f+D.f) - (nSdStay.f.none+D.f.none))) + 
  theme(axis.text=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank()) + 
  scale_fill_gradient2("Seed\nPressure\nFit-None", midpoint=0) + 
  facet_grid(SDM~issue) + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/diff_none/", sp, "_map_propagule.jpg"), p, width=12, height=4)

p <- ggplot(filter(out_P, SDM=="IPM")) + 
  geom_tile(aes(lon, lat, fill=lambda.f - lambda.f.none)) + 
  theme(axis.text=element_blank(),
        panel.background=element_rect(fill="gray80"),
        panel.grid=element_blank()) + 
  scale_fill_gradient2("lambda\nFit-None", midpoint=0) + 
  facet_grid(SDM~issue) + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/diff_none/", sp, "_map_lambda.jpg"), width=12, height=2.5)


# 
# # Scatter plots
# ggplot(out_P, aes(prP, log(Surv.S+1))) + 
#   labs(x="Probability of presence", y="True log(N+1)") + 
#   geom_point(alpha=0.02) + 
#   facet_grid(SDM~issue) + ggtitle(sp.names)
# if(plots) ggsave(paste0("eda/scatter/", sp, "_prP_N.jpg"), width=12, height=7)
# 
# ggplot(out_P, aes(prP, lambda)) + geom_point(alpha=0.02) + 
#   labs(x="Probability of presence", y="True lambda") + 
#   geom_hline(colour="red", linetype=3, yintercept=1) + 
#   facet_grid(SDM~issue) + ggtitle(sp.names)
# if(plots) ggsave(paste0("eda/scatter/", sp, "_prP_lambda.jpg"), width=12, height=7)
# 
# ggplot(filter(out_P, grepl("CA", SDM)), aes(Surv.S.f, Surv.S)) + 
#   labs(x="Predicted N", y="True N") + 
#   geom_point(alpha=0.02) + geom_abline(colour="red", linetype=3) +
#   facet_grid(SDM~issue) + xlim(range(c(out_P$Surv.S, out_P$Surv.S.f), na.rm=T)) +
#   ylim(range(c(out_P$Surv.S, out_P$Surv.S.f), na.rm=T)) + ggtitle(sp.names)
# if(plots) ggsave(paste0("eda/scatter/", sp, "_N.jpg"), width=12, height=4)
# 
# ggplot(filter(out_P, SDM=="IPM"), aes(lambda.f, lambda)) + 
#   labs(x="Predicted lambda", y="True lambda") + 
#   geom_point(alpha=0.02) + geom_abline(colour="red", linetype=3) +
#   facet_grid(SDM~issue) + xlim(range(c(out_P$lambda, out_P$lambda.f), na.rm=T)) +
#   ylim(range(c(out_P$lambda, out_P$lambda.f), na.rm=T)) + ggtitle(sp.names)
# if(plots) ggsave(paste0("eda/scatter/", sp, "_lambda.jpg"), width=12, height=2.5)
# 


# prP: Boxplots
ggplot(out_P, aes(lambda>1, prP, fill=SDM, colour=SDM)) +
  geom_boxplot(alpha=0.75, outlier.size=0.25, outlier.alpha=0.1) + 
  facet_wrap(.~issue) +
  scale_fill_manual("SDM\nMethod", values=SDM_col) + 
  scale_colour_manual(values=SDM_col, guide=F) +
  theme(legend.position=c(0.85, 0.15)) +
  labs(x="Presence/Absence: lambda > 1", y="Probability of presence") + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/eval/", sp, "_prP_lambda_box_byIssue.jpg"), width=6, height=6)

ggplot(out_P, aes(Surv.S>0, prP, fill=SDM, colour=SDM)) +
  geom_boxplot(alpha=0.75, outlier.size=0.25, outlier.alpha=0.1) + 
  facet_wrap(.~issue) +
  scale_fill_manual("SDM\nMethod", values=SDM_col) + 
  scale_colour_manual(values=SDM_col, guide=F) +
  theme(legend.position=c(0.85, 0.15)) +
  labs(x="Presence/Absence: N > 0", y="Probability of presence") + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/eval/", sp, "_prP_N_box_byIssue.jpg"), width=6, height=6)

ggplot(out_P, aes(lambda>1, prP, fill=issue, colour=issue)) +
  geom_boxplot(alpha=0.75, outlier.size=0.25, outlier.alpha=0.1) + 
  facet_wrap(.~SDM) +
  scale_fill_brewer(type="qual", palette=2) +
  scale_colour_brewer(type="qual", palette=2) +
  labs(x="Presence/Absence: lambda > 1", y="Probability of presence") + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/eval/", sp, "_prP_lambda_box_bySDM.jpg"), width=6, height=4)

ggplot(out_P, aes(Surv.S>0, prP, fill=issue, colour=issue)) +
  geom_boxplot(alpha=0.75, outlier.size=0.25, outlier.alpha=0.1) + 
  facet_wrap(.~SDM) +
  scale_fill_brewer(type="qual", palette=2) +
  scale_colour_brewer(type="qual", palette=2) +
  labs(x="Presence/Absence: N > 0", y="Probability of presence") + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/eval/", sp, "_prP_N_box_bySDM.jpg"), width=6, height=4)



# violin plots of prP differences: P vs. P.none
p <- ggplot(out_P, aes(issue, prP-prP.none, fill=issue, colour=issue)) + 
  geom_hline(yintercept=0, linetype=3, colour="grey80") +
  geom_violin(size=0.1, alpha=0.75) + 
  facet_wrap(~SDM) + ylim(-1,1) +
  scale_fill_brewer(type="qual", palette=2) +
  scale_colour_brewer(type="qual", palette=2) +  
  theme(panel.grid=element_blank(), 
        axis.text.x=element_text(angle=270, vjust=0.5, hjust=0)) +
  labs(x="Issue", y="Probability of presence") + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/diff_none/", sp, "_prP_violin_bySDM.jpg"), p, width=6, height=6)

p <- ggplot(out_P, aes(SDM, prP-prP.none, fill=SDM, colour=SDM)) + 
  geom_hline(yintercept=0, linetype=3, colour="grey80") +
  geom_violin(size=0.1, alpha=0.75) + 
  facet_wrap(~issue) + ylim(-1,1) +
  scale_fill_manual("SDM\nMethod", values=SDM_col) + 
  scale_colour_manual(values=SDM_col, guide=F) +
  theme(panel.grid=element_blank(), 
        axis.text.x=element_text(angle=270, vjust=0.5, hjust=0)) +
  labs(x="SDM", y="Probability of presence") + ggtitle(sp.names)
if(plots) ggsave(paste0("eda/diff_none/", sp, "_prP_violin_byIssue.jpg"), p, width=6, height=6)



# MaxEnt diagnostics
ggplot(ROC_MxE, aes(FPR, TPR, group=obs)) + 
  geom_abline(linetype=3, colour="red", alpha=0.5) + 
  geom_line(alpha=0.1) + facet_grid(Boundary~issue)
ggplot(AUC_MxE, aes(AUC, colour=issue)) + geom_density() + facet_wrap(~Boundary)

