# 5: Plot output
# Comparison of SDM approaches
# Tim Szewczyk

# This script plots the output from `3_fitModels.R` and `4_analysis.R`.

########
## Setup
########
# file specifications
sp <- "sp1"

# load workspace
pkgs <- c("tidyverse", "magrittr", "stringr", "here", "viridis")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "fn", full.names=T), source)
theme_set(theme_bw())
lam.df <- readRDS(here("vs", sp, "lam_df.rds"))
out <- read.csv(here("out", sp, "out.csv"))
out$issue <- factor(out$issue, 
                    levels=c("none", "noise", "sampBias", "nonEq",
                             "noSB", "underDisp", "overDisp", "wrongCov"),
                    labels=c("None", "Measurement error", "Sampling bias", 
                             "Non-equilibrium", "No Seedbank", "Under dispersal",
                             "Over dispersal", "Incorrect covariates"))
SDM_col <- c(MxE="#3f007d", IPM="#014636", CAi="#02818a", CAd="#67a9cf")

par(mfrow=c(3,3))
plot(lam.df$bio10_1, log(lam.df$lambda), col=rgb(0,0,0,0.75))
plot(lam.df$bio10_12, log(lam.df$lambda), col=rgb(0,0,0,0.75))
plot(lam.df$bio10_5, log(lam.df$lambda), col=rgb(0,0,0,0.75))
plot(lam.df$bio10_prMay, log(lam.df$lambda), col=rgb(0,0,0,0.75))
hist(log(lam.df$lambda))
plot(lam.df$bio10_1, log(lam.df$Surv.S), col=rgb(0,0,0,0.75))
plot(lam.df$bio10_12, log(lam.df$Surv.S), col=rgb(0,0,0,0.75))
plot(lam.df$bio10_5, log(lam.df$Surv.S), col=rgb(0,0,0,0.75))
plot(lam.df$bio10_prMay, log(lam.df$Surv.S), col=rgb(0,0,0,0.75))

ggplot(lam.df, aes(x=lon, y=lat, fill=lambda>1)) + geom_tile() + ggtitle(sp)
ggplot(lam.df, aes(x=lon, y=lat, fill=log(lambda))) + geom_tile() + ggtitle(sp) +
  scale_fill_viridis(option="B")
ggplot(lam.df, aes(x=lon, y=lat, fill=Surv.S>0)) + geom_tile() + ggtitle(sp)
ggplot(lam.df, aes(x=lon, y=lat, fill=log(Surv.S))) + geom_tile() + ggtitle(sp) +
  scale_fill_viridis(option="B")

ggplot(out, aes(fill=fate_lam, x=SDM)) + geom_bar(position="fill") + 
  facet_wrap(~issue) + scale_fill_brewer(name="", type="div") + 
  ylab("Proportion of cells") + ggtitle(sp, "lambda-based range")
ggplot(out, aes(fill=fate_S, x=SDM)) + geom_bar(position="fill") + 
  facet_wrap(~issue) + scale_fill_brewer(name="", type="div") + 
  ylab("Proportion of cells") + ggtitle(sp, "abundance-based range")
ggplot(out, aes(fill=fate_lam, x=issue)) + geom_bar(position="fill") + 
  facet_wrap(~SDM) + scale_fill_brewer(name="", type="div") + 
  ylab("Proportion of cells") + coord_flip() + ggtitle(sp, "lambda-based range")
ggplot(out, aes(fill=fate_S, x=issue)) + geom_bar(position="fill") + 
  facet_wrap(~SDM) + scale_fill_brewer(name="", type="div") + 
  ylab("Proportion of cells") + coord_flip() + ggtitle(sp, "abundance-based range")

ggplot(out, aes(x=lon, y=lat, fill=fate_lam)) +
  geom_tile() + facet_grid(SDM~issue) + 
  scale_fill_brewer(name="Boundary:\nlambda ≥ 1", type="div") +
  theme(axis.text=element_blank()) + ggtitle(sp)
ggplot(out, aes(x=lon, y=lat, fill=fate_S)) +
  geom_tile() + facet_grid(SDM~issue) + 
  scale_fill_brewer(name="Boundary:\nN > 0", type="div") +
  theme(axis.text=element_blank()) + ggtitle(sp)

ggplot(out, aes(x=lon, y=lat, fill=prP)) +
  geom_tile() + facet_grid(SDM~issue) + 
  scale_fill_viridis("prob(P)", option="B") +
  theme(axis.text=element_blank()) + ggtitle(sp)
ggplot(out, aes(x=lon, y=lat, fill=prP>0.5)) +
  geom_tile() + facet_grid(SDM~issue) + 
  theme(axis.text=element_blank())
ggplot(out, aes(x=lon, y=lat, fill=lambda.f>=1)) +
  geom_tile() + facet_grid(SDM~issue) + 
  theme(axis.text=element_blank())
ggplot(out, aes(x=lon, y=lat, fill=log(Surv.S.f))) +
  geom_tile() + facet_grid(SDM~issue) + 
  scale_fill_viridis("log(N)", option="B") +
  theme(axis.text=element_blank())
ggplot(out, aes(x=lon, y=lat, fill=log(lambda.f))) +
  geom_tile() + facet_grid(SDM~issue) + 
  scale_fill_viridis("log(lambda)", option="B") +
  theme(axis.text=element_blank())
ggplot(out, aes(x=lon, y=lat, fill=log(B.f))) +
  geom_tile() + facet_grid(SDM~issue) + 
  scale_fill_viridis("log(B)", option="B") +
  theme(axis.text=element_blank())
ggplot(out, aes(x=lon, y=lat, fill=log(D.f))) +
  geom_tile() + facet_grid(SDM~issue) + 
  scale_fill_viridis("log(D)", option="B") +
  theme(axis.text=element_blank())
ggplot(out, aes(x=lon, y=lat, fill=log(nSdStay.f + nSdLeave.f))) +
  geom_tile() + facet_grid(SDM~issue) + 
  scale_fill_viridis("log(Sd)", option="B") +
  theme(axis.text=element_blank())

lam.sum <- out %>% group_by(SDM, issue, fate_lam) %>%
  summarise(ct=n()) %>%
  mutate(rate=case_when(fate_lam=="S:0 P:0" ~ ct/sum(lam.df$lambda<1),
                            fate_lam=="S:0 P:1" ~ ct/sum(lam.df$lambda<1),
                            fate_lam=="S:1 P:0" ~ ct/sum(lam.df$lambda>=1),
                            fate_lam=="S:1 P:1" ~ ct/sum(lam.df$lambda>=1)))
S.sum <- out %>% group_by(SDM, issue, fate_S) %>%
  summarise(ct=n()) %>%
  mutate(rate=case_when(fate_S=="S:0 P:0" ~ ct/sum(lam.df$Surv.S==0),
                            fate_S=="S:0 P:1" ~ ct/sum(lam.df$Surv.S==0),
                            fate_S=="S:1 P:0" ~ ct/sum(lam.df$Surv.S>0),
                            fate_S=="S:1 P:1" ~ ct/sum(lam.df$Surv.S>0)))
tss.lam <- lam.sum %>% ungroup %>% group_by(SDM, issue) %>%
  summarise(TSS=sum(rate[fate_lam %in% c("S:0 P:0", "S:1 P:1")])-1) %>%
  ungroup %>% mutate(issue=fct_rev(issue), metric="lambda > 1")
tss.S <- S.sum %>% ungroup %>% group_by(SDM, issue) %>%
  summarise(TSS=sum(rate[fate_S %in% c("S:0 P:0", "S:1 P:1")])-1) %>%
  ungroup %>% mutate(issue=fct_rev(issue), metric="N > 0")
tss.df <- rbind(tss.lam, tss.S)
ggplot(tss.df, aes(x=TSS, y=issue, shape=metric)) +
  labs(title=sp, y="") + xlim(0,1) + 
  geom_point(size=5, alpha=0.9, fill=NA, colour="black") +
  geom_point(size=5, alpha=0.7, aes(fill=SDM, colour=SDM)) +
  geom_vline(xintercept=1, colour="gray") + 
  geom_vline(xintercept=0, colour="gray", linetype=2) +
  theme(panel.grid.major.y=element_line(colour="gray")) +
  scale_fill_manual(values=SDM_col) + scale_colour_manual(values=SDM_col) +
  scale_shape_manual(values=c(21, 23))


par(mfrow=c(5,6))
map(list.files("out", "Diag_Mx", full.names=T), readRDS) %>%
  walk(~walk(., ~plot(., 'ROC')))
ggplot(out, aes(x=prP, y=1*(Surv.S>0))) + geom_point(alpha=0.05) + 
  facet_grid(SDM~issue)
  


out %>% filter(fate_lam == "S:0 P:1") %>% group_by(SDM, issue) %>% 
  summarise(prop=round(n()/1084, 3)) %>% 
  ggplot(aes(x=issue, y=prop, colour=SDM)) + geom_point(size=3) +
  ggtitle("Commission Rate") + ylim(0,1) + ylab("Prop P=1 | S=0")
out %>% filter(fate_lam == "S:1 P:0") %>%  group_by(SDM, issue) %>% 
  summarise(prop=round(n()/1354, 3)) %>%
  ggplot(aes(x=issue, y=prop, colour=SDM)) + geom_point(size=3) + 
  ggtitle("Ommission Rate") + ylim(0, 1) + ylab("Prop P=0 | S=1")
out %>% filter(fate_lam %in% c("S:1 P:1", "S:0 P:0")) %>% 
  group_by(SDM, issue) %>% 
  summarise(prop=round(n()/2438, 3)) %>% 
  ggplot(aes(x=issue, y=prop, colour=SDM)) + geom_point(size=3) + 
  ggtitle("TSS") + ylim(0, 1) + ylab("Prop P=S")

ggplot(out, aes(fill=sign(log(lambda.f))==sign(log(lambda)), x=SDM)) + 
  geom_bar(position="fill") + facet_wrap(~issue)

ggplot(out, aes(x=lon, y=lat, fill=Surv.S.f-Surv.S)) +
  geom_tile() + facet_grid(SDM~issue) + scale_fill_gradient2()
ggplot(out, aes(x=lon, y=lat, fill=(Surv.S.f-Surv.S)/Surv.S)) +
  geom_tile() + facet_grid(SDM~issue) + scale_fill_gradient()



ggplot(out, aes(x=lon, y=lat)) + facet_grid(SDM~issue) +
  geom_tile(aes(fill=Surv.S.f>1)) + 
  scale_fill_manual(values=c(NA, "dodgerblue")) +
  geom_point(aes(colour=Surv.S>0), alpha=0.3, size=0.6) + 
  scale_colour_manual(values=c(NA, "black"))
ggplot(out, aes(x=lon, y=lat, fill=log(lam.S.f))) +
  geom_tile() + facet_grid(SDM~issue) + scale_fill_gradient2()
ggplot(out, aes(x=lon, y=lat, fill=log(lam.S.f)>=0)) +
  geom_tile() + facet_grid(SDM~issue)

ggplot(out, aes(x=lon, y=lat, fill=log(Surv.S.f))) + geom_tile() + 
  facet_grid(SDM~issue) + scale_fill_gradient(low="white", high="red")

ggplot(out, aes(x=lon, y=lat, fill=sign(log(lam.S.f))==sign(log(lam.S)))) +
  geom_tile() + facet_grid(SDM~issue)

ggplot(out, aes(x=lon, y=lat, fill=sign(log(lambda.f))==sign(log(lambda)))) +
  geom_tile() + facet_grid(SDM~issue)

ggplot(filter(out, SDM=="IPM"), aes(x=lon, y=lat, fill=lambda.f-lambda)) +
  geom_tile() + facet_wrap(~issue) + 
  scale_fill_gradient2(low="blue", high="red")
ggplot(filter(out, SDM=="IPM"), aes(x=lon, y=lat, fill=lambda.f)) +
  geom_tile() + facet_wrap(~issue) + 
  scale_fill_gradient2(low="blue", high="red", midpoint=1, limits=c(0,NA))
ggplot(filter(out, SDM=="IPM"), aes(x=lon, y=lat, fill=lambda)) +
  geom_tile() + facet_wrap(~issue) + 
  scale_fill_gradient2(low="blue", high="red", midpoint=1, limits=c(0,4))

ggplot(out, aes(x=lon, y=lat, fill=D.f-D)) +
  geom_tile() + facet_grid(SDM~issue) + scale_fill_gradient2()

ggplot(out, aes(x=lon, y=lat, fill=Btmax.f-Btmax)) +
  geom_tile() + facet_grid(SDM~issue) + scale_fill_gradient2()

ggplot(out, aes(x=lon, y=lat, fill=nSeed.f-nSeed)) +
  geom_tile() + facet_grid(SDM~issue) + scale_fill_gradient2()

ggplot(out, aes(x=lon, y=lat, fill=Rcr.S.f-Rcr.S)) +
  geom_tile() + facet_grid(SDM~issue) + scale_fill_gradient2()

ggplot(out, aes(x=N.S, y=N.S.f, colour=SDM, group=SDM)) + 
  geom_point(alpha=0.5) + facet_wrap(~issue) +
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)

ggplot(out, aes(x=Surv.S, y=Surv.S.f, colour=SDM, group=SDM)) + 
  geom_point(alpha=0.5) + facet_wrap(~issue, scales="free") +
  scale_colour_manual(values=SDM_col) + 
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)

ggplot(out, aes(x=Rcr.S, y=Rcr.S.f, colour=SDM, group=SDM)) + 
  geom_point(alpha=0.5) + facet_wrap(~issue, scales="free") +
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)
ggplot(out, aes(x=Rcr.S.f-Rcr.S, colour=SDM)) + geom_density() + facet_wrap(~issue)

ggplot(filter(out, SDM=="IPM"), aes(x=lambda, y=lambda.f)) + 
  geom_point(alpha=0.5) + facet_wrap(~issue) +
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)

ggplot(out, aes(x=lam.S, y=lam.S.f)) + geom_point(alpha=0.5) + 
  facet_grid(SDM~issue) + scale_colour_manual(values=SDM_col) + 
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)

ggplot(filter(out, SDM=="IPM"), aes(x=lam.U, y=lam.U.f)) + 
  geom_point(alpha=0.5) + facet_wrap(~issue) +
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)

ggplot(filter(out, SDM=="IPM"), aes(x=D, y=D.f)) + geom_point(alpha=0.5) + 
  facet_wrap(~issue) +
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)

ggplot(filter(out, SDM=="IPM"), aes(x=nSdStay, y=nSdStay.f)) + geom_point(alpha=0.5) + 
  facet_wrap(~issue) +
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)
