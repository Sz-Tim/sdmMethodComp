# 5: Plot output
# Comparison of SDM approaches
# Tim Szewczyk

# This script plots the output from `3_fitModels.R` and `4_analysis.R`.

########
## Setup
########
# file specifications
sp <- "sp1"
plots <- TRUE

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


ggplot(out_TSS, aes(x=TSS, fill=SDM, colour=SDM)) + 
  geom_vline(xintercept=c(-1,0,1), colour="grey90") +
  geom_vline(xintercept=c(-0.5,0.5), colour="grey90", linetype=3) +
  geom_density(alpha=0.75) + 
  geom_rug(data=mn_TSS, size=1) +
  facet_grid(issue~Boundary, scales="free_y", 
             labeller=labeller(issue=label_wrap_gen(10), Boundary=label_parsed)) +
  scale_fill_manual("SDM\nMethod", values=SDM_col) + 
  scale_colour_manual(values=SDM_col, guide=F) +
  scale_x_continuous("TSS", c(-1,0,1)) + ylab("Density") +
  theme(panel.grid=element_blank())
if(plots) ggsave(paste0("figs/eval/", sp, "_TSS.jpg"), width=6, height=8)
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
  theme(panel.grid=element_blank())
if(plots) ggsave(paste0("figs/eval/", sp, "_Sensitivity.jpg"), width=6, height=8)
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
  theme(panel.grid=element_blank())
if(plots) ggsave(paste0("figs/eval/", sp, "_Specificity.jpg"), width=6, height=8)


ggplot(out_LL, aes(x=LogLik, y=fct_rev(issue), colour=SDM)) + 
  facet_grid(Boundary~.) +
  geom_point(alpha=0.75, size=3) + 
  scale_colour_manual(values=SDM_col) +
  labs(x="Log likelihood: S ~ Binom(prP)", y="")
if(plots) ggsave(paste0("figs/eval/", sp, "_LogLik.jpg"), width=5, height=4)


if(plots) jpeg(paste0("figs/S/", sp, "_scatter.jpg"), 
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


ggplot() + geom_tile(data=lam.df, aes(lon, lat), fill="gray30") +
  geom_tile(data=filter(lam.df, lambda>1), aes(lon, lat, fill=lambda)) +
  scale_fill_viridis(option="B") + 
  ggtitle(paste0(ifelse(sp=="sp1", "barberry", "garlic mustard"), ": lambda"))
if(plots) ggsave(paste0("figs/S/", sp, "_lambda.jpg"), width=9, height=8)
ggplot() + geom_tile(data=lam.df, aes(lon, lat), fill="gray30") +
  geom_tile(data=filter(lam.df, Surv.S>0), aes(lon, lat, fill=Surv.S)) +
  scale_fill_viridis("N", option="B") + 
  ggtitle(paste0(ifelse(sp=="sp1", "barberry", "garlic mustard"), ": N.tmax"))
if(plots) ggsave(paste0("figs/S/", sp, "_N_Eq.jpg"), width=9, height=8)
ggplot() + geom_tile(data=lam.df, aes(lon, lat), fill="gray30") +
  geom_tile(data=filter(lam.df, Surv.S_nonEq>0), aes(lon, lat, fill=Surv.S_nonEq)) +
  scale_fill_viridis("N", option="B") + 
  ggtitle(paste0(ifelse(sp=="sp1", "barberry", "garlic mustard"), ": N.nonEq"))
if(plots) ggsave(paste0("figs/S/", sp, "_N_nonEq.jpg"), width=9, height=8)

ggplot() + theme_bw() +
  geom_tile(data=lam.df, aes(x=lon, y=lat), fill="gray90") +
  geom_tile(data=filter(lam.df, lambda>1), aes(x=lon, y=lat), 
            fill="darkblue", alpha=0.5) +
  geom_tile(data=filter(lam.df, Surv.S>0), aes(x=lon, y=lat), 
            fill="red", alpha=0.5) +
  theme(axis.text=element_blank()) + labs(x="", y="") +
  ggtitle(paste0(sp, ": 5km x 5km, favorable habitat"))


ggplot(out_P, aes(fill=fate_lam, x=issue)) + geom_bar(position="fill") + 
  facet_wrap(~SDM) + scale_fill_brewer(name="", type="div") + 
  ylab("Proportion of cells") + coord_flip() + ggtitle(sp, "lambda-based range")
if(plots) ggsave(paste0("figs/thresh/", sp, "_pr50_bar_lambda.jpg"), width=7, height=8)
ggplot(out_P, aes(fill=fate_S, x=issue)) + geom_bar(position="fill") + 
  facet_wrap(~SDM) + scale_fill_brewer(name="", type="div") + 
  ylab("Proportion of cells") + coord_flip() + ggtitle(sp, "abundance-based range")
if(plots) ggsave(paste0("figs/thresh/", sp, "_pr50_bar_N.jpg"), width=7, height=8)


ggplot(out_P, aes(x=lon, y=lat, fill=fate_lam)) +
  geom_tile() + facet_grid(SDM~issue) + 
  scale_fill_brewer(name="Boundary:\nlambda ≥ 1", type="div") +
  theme(axis.text=element_blank()) + ggtitle(sp)
if(plots) ggsave(paste0("figs/thresh/", sp, "_pr50_map_lambda.jpg"), width=12, height=7)
ggplot(out_P, aes(x=lon, y=lat, fill=fate_S)) +
  geom_tile() + facet_grid(SDM~issue) + 
  scale_fill_brewer(name="Boundary:\nN > 0", type="div") +
  theme(axis.text=element_blank()) + ggtitle(sp)
if(plots) ggsave(paste0("figs/thresh/", sp, "_pr50_map_N.jpg"), width=12, height=7)


ggplot(out_P, aes(x=lon, y=lat, fill=prP)) +
  geom_tile() + facet_grid(SDM~issue) + 
  scale_fill_viridis("prob(P)", option="B") +
  theme(axis.text=element_blank()) + ggtitle(sp)
if(plots) ggsave(paste0("figs/P/", sp, "_prP_map.jpg"), width=12, height=7)
ggplot() + geom_tile(data=filter(out_P, SDM=="IPM"), 
                     aes(lon, lat), fill="gray30") +
  geom_tile(data=filter(out_P, SDM=="IPM" & lambda.f>1), 
            aes(lon, lat, fill=lambda.f)) +
  facet_grid(SDM~issue) + 
  scale_fill_viridis("lambda", option="B") +
  theme(axis.text=element_blank())
if(plots) ggsave(paste0("figs/P/", sp, "_map_lambda.jpg"), width=12, height=2.5)
ggplot() + geom_tile(data=filter(out_P, grepl("CA", SDM)), 
                     aes(x=lon, y=lat), fill="gray30") +
  geom_tile(data=filter(out_P, grepl("CA", SDM) & Surv.S.f>0.5), 
            aes(lon, lat, fill=Surv.S.f)) +
  facet_grid(SDM~issue) + 
  scale_fill_viridis("N", option="B") +
  theme(axis.text=element_blank())
if(plots) ggsave(paste0("figs/", sp, "_P_map_N.jpg"), width=12, height=4)
ggplot() + geom_tile(data=filter(out_P, grepl("CA", SDM)), 
                     aes(x=lon, y=lat), fill="gray30") +
  geom_tile(data=filter(out_P, grepl("CA", SDM) & (Rcr.S.f)>0.5), 
            aes(lon, lat, fill=Rcr.S.f)) +
  facet_grid(SDM~issue) + 
  scale_fill_viridis("Recruits", option="B") +
  theme(axis.text=element_blank())
if(plots) ggsave(paste0("figs/P/", sp, "_map_Rcr.jpg"), width=12, height=4)
ggplot() + geom_tile(data=filter(out_P, grepl("CA", SDM)), 
                     aes(x=lon, y=lat), fill="gray30") +
  geom_tile(data=filter(out_P, grepl("CA", SDM) & B.f>0.5), 
            aes(lon, lat, fill=B.f)) +
  facet_grid(SDM~issue) + 
  scale_fill_viridis("B", option="B") +
  theme(axis.text=element_blank())
if(plots) ggsave(paste0("figs/P/", sp, "_map_B.jpg"), width=12, height=4)
ggplot() +
  geom_tile(data=filter(out_P, grepl("CA", SDM)), 
            aes(x=lon, y=lat), fill="gray30") +
  geom_tile(data=filter(out_P, grepl("CA", SDM) & D.f>0.5), 
            aes(lon, lat, fill=D.f)) +
  facet_grid(SDM~issue) + 
  scale_fill_viridis("D", option="B") +
  theme(axis.text=element_blank())
if(plots) ggsave(paste0("figs/P/", sp, "_map_D.jpg"), width=12, height=4)
ggplot() +
  geom_tile(data=filter(out_P, grepl("CA", SDM)), 
            aes(x=lon, y=lat), fill="gray30") +
  geom_tile(data=filter(out_P, grepl("CA", SDM) & (nSdStay.f+nSdLeave.f)>0.5), 
            aes(lon, lat, fill=nSdStay.f+nSdLeave.f)) +
  facet_grid(SDM~issue) + 
  scale_fill_viridis("SdTot", option="B") +
  theme(axis.text=element_blank())
if(plots) ggsave(paste0("figs/P/", sp, "_map_nSeed.jpg"), width=12, height=4)
ggplot() +
  geom_tile(data=filter(out_P, grepl("CA", SDM)), 
            aes(x=lon, y=lat), fill="gray30") +
  geom_tile(data=filter(out_P, grepl("CA", SDM) & (nSdStay.f)>0.5), 
            aes(lon, lat, fill=nSdStay.f)) +
  facet_grid(SDM~issue) + 
  scale_fill_viridis("SdStay", option="B") +
  theme(axis.text=element_blank())
if(plots) ggsave(paste0("figs/P/", sp, "_map_nSeedLocal.jpg"), width=12, height=4)
ggplot() +
  geom_tile(data=filter(out_P, grepl("CA", SDM)), 
            aes(x=lon, y=lat), fill="gray30") +
  geom_tile(data=filter(out_P, grepl("CA", SDM) & (nSdStay.f+D.f)>0.5), 
            aes(lon, lat, fill=nSdStay.f+D.f)) +
  facet_grid(SDM~issue) + 
  scale_fill_viridis("Seed\nPressure", option="B") +
  theme(axis.text=element_blank())
if(plots) ggsave(paste0("figs/P/", sp, "_map_propagule.jpg"), width=12, height=4)


ggplot(filter(out_P, grepl("CA", SDM))) + 
  geom_tile(aes(lon, lat, fill=Surv.S.f - Surv.S)) + 
  scale_fill_gradient2("N\nFit-True", midpoint=0) + facet_grid(SDM~issue)
if(plots) ggsave(paste0("figs/", sp, "_map_diff_N.jpg"), width=12, height=4)
ggplot(filter(out_P, grepl("CA", SDM))) + 
  geom_tile(aes(lon, lat, fill=nSeed.f - nSeed)) + 
  scale_fill_gradient2("nSeed\nFit-True", midpoint=0) + facet_grid(SDM~issue)
if(plots) ggsave(paste0("figs/", sp, "_map_diff_nSeed.jpg"), width=12, height=4)
ggplot(filter(out_P, SDM=="IPM")) + 
  geom_tile(aes(lon, lat, fill=lambda.f - lambda)) + 
  scale_fill_gradient2("lambda\nFit-True", midpoint=0) + facet_grid(SDM~issue)
if(plots) ggsave(paste0("figs/P/", sp, "_map_diff_lambda.jpg"), width=12, height=2.5)


ggplot(out_P, aes(prP, log(Surv.S+1))) + 
  labs(x="Probability of presence", y="True log(N+1)") + 
  geom_point(alpha=0.02) + facet_grid(SDM~issue)
if(plots) ggsave(paste0("figs/scatter/", sp, "_prP_N.jpg"), width=12, height=7)
ggplot(out_P, aes(prP, lambda)) + geom_point(alpha=0.02) + 
  labs(x="Probability of presence", y="True lambda") + 
  geom_hline(colour="red", linetype=3, yintercept=1) + facet_grid(SDM~issue)
if(plots) ggsave(paste0("figs/scatter/", sp, "_prP_lambda.jpg"), width=12, height=7)
ggplot(filter(out_P, grepl("CA", SDM)), aes(Surv.S.f, Surv.S)) + 
  labs(x="Predicted N", y="True N") + 
  geom_point(alpha=0.02) + geom_abline(colour="red", linetype=3) +
  facet_grid(SDM~issue) + xlim(range(c(out_P$Surv.S, out_P$Surv.S.f), na.rm=T)) +
  ylim(range(c(out_P$Surv.S, out_P$Surv.S.f), na.rm=T))
if(plots) ggsave(paste0("figs/scatter/", sp, "_N.jpg"), width=12, height=4)
ggplot(filter(out_P, SDM=="IPM"), aes(lambda.f, lambda)) + 
  labs(x="Predicted lambda", y="True lambda") + 
  geom_point(alpha=0.02) + geom_abline(colour="red", linetype=3) +
  facet_grid(SDM~issue) + xlim(range(c(out_P$lambda, out_P$lambda.f), na.rm=T)) +
  ylim(range(c(out_P$lambda, out_P$lambda.f), na.rm=T))
if(plots) ggsave(paste0("figs/scatter/", sp, "_lambda.jpg"), width=12, height=2.5)





break()
################################################################################
par(mfrow=c(5,6))
map(list.files("out", "Diag_Mx", full.names=T), readRDS) %>%
  walk(~walk(., ~plot(., 'ROC')))
ggplot(out_P, aes(x=prP, y=1*(Surv.S>0))) + geom_point(alpha=0.05) + 
  facet_grid(SDM~issue)



ggplot(out_P, aes(fill=sign(log(lambda.f))==sign(log(lambda)), x=SDM)) + 
  geom_bar(position="fill") + facet_wrap(~issue)

ggplot(out_P, aes(x=lon, y=lat, fill=Surv.S.f-Surv.S)) +
  geom_tile() + facet_grid(SDM~issue) + scale_fill_gradient2()
ggplot(out_P, aes(x=lon, y=lat, fill=(Surv.S.f-Surv.S)/Surv.S)) +
  geom_tile() + facet_grid(SDM~issue) + scale_fill_gradient()



ggplot(out_P, aes(x=lon, y=lat)) + facet_grid(SDM~issue) +
  geom_tile(aes(fill=Surv.S.f>1)) + 
  scale_fill_manual(values=c(NA, "dodgerblue")) +
  geom_point(aes(colour=Surv.S>0), alpha=0.3, size=0.6) + 
  scale_colour_manual(values=c(NA, "black"))
ggplot(out_P, aes(x=lon, y=lat, fill=log(lam.S.f))) +
  geom_tile() + facet_grid(SDM~issue) + scale_fill_gradient2()
ggplot(out_P, aes(x=lon, y=lat, fill=log(lam.S.f)>=0)) +
  geom_tile() + facet_grid(SDM~issue)

ggplot(out_P, aes(x=lon, y=lat, fill=log(Surv.S.f))) + geom_tile() + 
  facet_grid(SDM~issue) + scale_fill_gradient(low="white", high="red")

ggplot(out_P, aes(x=lon, y=lat, fill=sign(log(lam.S.f))==sign(log(lam.S)))) +
  geom_tile() + facet_grid(SDM~issue)

ggplot(out_P, aes(x=lon, y=lat, fill=sign(log(lambda.f))==sign(log(lambda)))) +
  geom_tile() + facet_grid(SDM~issue)

ggplot(filter(out_P, SDM=="IPM"), aes(x=lon, y=lat, fill=lambda.f-lambda)) +
  geom_tile() + facet_wrap(~issue) + 
  scale_fill_gradient2(low="blue", high="red")
ggplot(filter(out_P, SDM=="IPM"), aes(x=lon, y=lat, fill=lambda.f)) +
  geom_tile() + facet_wrap(~issue) + 
  scale_fill_gradient2(low="blue", high="red", midpoint=1, limits=c(0,NA))
ggplot(filter(out_P, SDM=="IPM"), aes(x=lon, y=lat, fill=lambda)) +
  geom_tile() + facet_wrap(~issue) + 
  scale_fill_gradient2(low="blue", high="red", midpoint=1, limits=c(0,4))

ggplot(out_P, aes(x=lon, y=lat, fill=D.f-D)) +
  geom_tile() + facet_grid(SDM~issue) + scale_fill_gradient2()

ggplot(out_P, aes(x=lon, y=lat, fill=Btmax.f-Btmax)) +
  geom_tile() + facet_grid(SDM~issue) + scale_fill_gradient2()

ggplot(out_P, aes(x=lon, y=lat, fill=nSeed.f-nSeed)) +
  geom_tile() + facet_grid(SDM~issue) + scale_fill_gradient2()

ggplot(out_P, aes(x=lon, y=lat, fill=Rcr.S.f-Rcr.S)) +
  geom_tile() + facet_grid(SDM~issue) + scale_fill_gradient2()

ggplot(out_P, aes(x=N.S, y=N.S.f, colour=SDM, group=SDM)) + 
  geom_point(alpha=0.5) + facet_wrap(~issue) +
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)

ggplot(out_P, aes(x=Surv.S, y=Surv.S.f, colour=SDM, group=SDM)) + 
  geom_point(alpha=0.5) + facet_wrap(~issue, scales="free") +
  scale_colour_manual(values=SDM_col) + 
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)

ggplot(out_P, aes(x=Rcr.S, y=Rcr.S.f, colour=SDM, group=SDM)) + 
  geom_point(alpha=0.5) + facet_wrap(~issue, scales="free") +
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)
ggplot(out_P, aes(x=Rcr.S.f-Rcr.S, colour=SDM)) + geom_density() + facet_wrap(~issue)

ggplot(filter(out_P, SDM=="IPM"), aes(x=lambda, y=lambda.f)) + 
  geom_point(alpha=0.5) + facet_wrap(~issue) +
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)

ggplot(out_P, aes(x=lam.S, y=lam.S.f)) + geom_point(alpha=0.5) + 
  facet_grid(SDM~issue) + scale_colour_manual(values=SDM_col) + 
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)

ggplot(filter(out_P, SDM=="IPM"), aes(x=lam.U, y=lam.U.f)) + 
  geom_point(alpha=0.5) + facet_wrap(~issue) +
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)

ggplot(filter(out_P, SDM=="IPM"), aes(x=D, y=D.f)) + geom_point(alpha=0.5) + 
  facet_wrap(~issue) +
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)

ggplot(filter(out_P, SDM=="IPM"), aes(x=nSdStay, y=nSdStay.f)) + geom_point(alpha=0.5) + 
  facet_wrap(~issue) +
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)
