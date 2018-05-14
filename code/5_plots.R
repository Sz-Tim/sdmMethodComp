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
pkgs <- c("tidyverse", "magrittr", "stringr", "here")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(paste0("code/fn_", c("IPM", "aux", "sim"), ".R"), ~source(here(.)))
out <- read.csv(here(paste0("out/", sp, "_out.csv")))
out$issue <- factor(out$issue, 
                    levels=c("none", "noise", "geogBias", "sampBias",
                             "noSB", "noDisp", "overDisp", "clim", "lc"),
                    labels=c("None", "Measurement error", "Geographic bias", 
                             "Sampling bias", "No Seedbank", "No dispersal",
                             "Over dispersal", "Climate Only", "LC only"))

par(mfrow=c(3,3))
plot(lam.df$temp, log(lam.df$lambda), col=rgb(0,0,0,0.75))
plot(lam.df$prec, log(lam.df$lambda), col=rgb(0,0,0,0.75))
plot(lam.df$pOpn, log(lam.df$lambda), col=rgb(0,0,0,0.75))
plot(lam.df$pOth, log(lam.df$lambda), col=rgb(0,0,0,0.75))
plot(lam.df$pDec, log(lam.df$lambda), col=rgb(0,0,0,0.75))
plot(lam.df$pEvg, log(lam.df$lambda), col=rgb(0,0,0,0.75))
plot(lam.df$pMxd, log(lam.df$lambda), col=rgb(0,0,0,0.75))
hist(log(lam.df$lambda))

theme_set(theme_bw())
ggplot(out, aes(fill=outcome, x=SDM)) + geom_bar(position="fill") + 
  facet_wrap(~issue) + scale_fill_brewer(name="", type="div") + 
  ylab("Proportion of cells")
ggplot(out, aes(fill=outcome, x=issue)) + geom_bar(position="fill") + 
  facet_wrap(~SDM) + scale_fill_brewer(name="", type="div") + 
  ylab("Proportion of cells")
ggplot(out, aes(x=lon, y=lat, fill=prP)) +
  geom_tile() + facet_grid(SDM~issue) + 
  scale_fill_gradient(low="white", high="red") +
  theme(axis.text=element_blank())
ggplot(out, aes(x=lon, y=lat, fill=prP.sd)) +
  geom_tile() + facet_grid(SDM~issue) + 
  scale_fill_gradient(low="blue", high="white") +
  theme(axis.text=element_blank())
ggplot(out, aes(x=lon, y=lat, fill=prP>0.5)) +
  geom_tile() + facet_grid(SDM~issue) + 
  theme(axis.text=element_blank())
ggplot(out, aes(x=lon, y=lat, fill=outcome)) +
  geom_tile() + facet_grid(SDM~issue) + scale_fill_brewer(name="", type="div") +
  theme(axis.text=element_blank())
ggplot(out, aes(x=lon, y=lat, fill=log(Btmax.f))) +
  geom_tile() + facet_grid(SDM~issue) + 
  scale_fill_gradient(low="white", high="red") +
  theme(axis.text=element_blank())
ggplot(out, aes(x=lon, y=lat, fill=log(D.f))) +
  geom_tile() + facet_grid(SDM~issue) + 
  scale_fill_gradient(low="white", high="red") +
  theme(axis.text=element_blank())

out.sum <- out %>% group_by(SDM, issue, outcome) %>%
  summarise(ct=n()) %>%
  mutate(rate=case_when(outcome=="S:0 P:0" ~ ct/sum(lam.df$Surv.S==0),
                        outcome=="S:0 P:1" ~ ct/sum(lam.df$Surv.S==0),
                        outcome=="S:1 P:0" ~ ct/sum(lam.df$Surv.S > 0),
                        outcome=="S:1 P:1" ~ ct/sum(lam.df$Surv.S > 0)))
tss.df <- out.sum %>% ungroup %>% group_by(SDM, issue) %>%
  summarise(TSS=sum(rate[outcome %in% c("S:0 P:0", "S:1 P:1")])-1)
ggplot(tss.df, aes(x=issue, colour=SDM, y=TSS)) + geom_point()

out %>% filter(outcome == "S:0 P:1") %>% group_by(SDM, issue) %>% 
  summarise(prop=round(n()/1084, 3)) %>% 
  ggplot(aes(x=issue, y=prop, colour=SDM)) + geom_point() +
  ggtitle("Commission Rate") + ylim(0,1) + ylab("Prop P=1 | S=0")
out %>% filter(outcome == "S:1 P:0") %>%  group_by(SDM, issue) %>% 
  summarise(prop=round(n()/1354, 3)) %>%
  ggplot(aes(x=issue, y=prop, colour=SDM)) + geom_point() + 
  ggtitle("Ommission Rate") + ylim(0, 1) + ylab("Prop P=0 | S=1")
out %>% filter(outcome %in% c("S:1 P:1", "S:0 P:0")) %>% 
  group_by(SDM, issue) %>% 
  summarise(prop=round(n()/2438, 3)) %>% 
  ggplot(aes(x=issue, y=prop, colour=SDM)) + geom_point() + 
  ggtitle("TSS") + ylim(0, 1) + ylab("Prop P=S")

ggplot(out_IPM, aes(fill=sign(log(lambda.f))==sign(log(lambda)), x=issue)) + 
  geom_bar(position="fill")
ggplot(out_IPM, aes(fill=sign(log(lam.U.f))==sign(log(lam.U)), x=issue)) + 
  geom_bar(position="fill")
ggplot(out, aes(fill=sign(log(lam.S.f))==sign(log(lam.S)), x=SDM)) + 
  geom_bar(position="fill") + facet_wrap(~issue)

ggplot(out, aes(x=lon, y=lat, fill=Surv.S.f-Surv.S)) +
  geom_tile() + facet_grid(SDM~issue) + scale_fill_gradient2()
ggplot(out, aes(x=lon, y=lat, fill=(Surv.S.f-Surv.S)/Surv.S)) +
  geom_tile() + facet_grid(SDM~issue) + scale_fill_gradient()


ggplot(out, aes(x=lon, y=lat, fill=Surv.S.f)) +
  geom_tile() + facet_grid(SDM~issue)
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

ggplot(out, aes(x=lon, y=lat, fill=sign(log(lam.S.f+.001))==sign(log(lam.S)))) +
  geom_tile() + facet_grid(SDM~issue)

ggplot(out_IPM, aes(x=lon, y=lat, fill=sign(log(lambda.f))==sign(log(lambda)))) +
  geom_tile() + facet_wrap(~issue)

ggplot(out_IPM, aes(x=lon, y=lat, fill=N.U.f-N.U)) +
  geom_tile() + facet_wrap(~issue) + scale_fill_gradient2()

ggplot(out_IPM, aes(x=lon, y=lat, fill=lam.U.f-lam.U)) +
  geom_tile() + facet_wrap(~issue) + scale_fill_gradient2()

ggplot(out, aes(x=lon, y=lat, fill=D.f-D)) +
  geom_tile() + facet_grid(SDM~issue) + scale_fill_gradient2()

ggplot(out, aes(x=lon, y=lat, fill=Btmax.f-Btmax)) +
  geom_tile() + facet_grid(SDM~issue) + scale_fill_gradient2()

ggplot(out, aes(x=lon, y=lat, fill=nSeed.f-nSeed)) +
  geom_tile() + facet_grid(SDM~issue) + scale_fill_gradient2()

ggplot(out, aes(x=lon, y=lat, fill=Rcr.S.f-Rcr.S)) +
  geom_tile() + facet_grid(SDM~issue) + scale_fill_gradient2()

ggplot(out_IPM, aes(x=lon, y=lat, fill=N.S.f-N.S)) +
  geom_tile() + facet_wrap(~issue) + scale_fill_gradient2()

ggplot(out, aes(x=N.S, y=N.S.f, colour=SDM, group=SDM)) + 
  geom_point(alpha=0.5) + facet_wrap(~issue) +
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)

ggplot(out, aes(x=Surv.S, y=Surv.S.f, colour=SDM, group=SDM)) + 
  geom_point(alpha=0.5) + facet_wrap(~issue) +
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)

ggplot(out, aes(x=Rcr.S, y=Rcr.S.f, colour=SDM, group=SDM)) + 
  geom_point(alpha=0.5) + facet_wrap(~issue) +
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)
ggplot(out, aes(x=Rcr.S.f-Rcr.S, colour=SDM)) + geom_density() + facet_wrap(~issue)

ggplot(out_IPM, aes(x=lambda, y=lambda.f)) + geom_point(alpha=0.5) + 
  facet_wrap(~issue) +
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)

ggplot(out, aes(x=lam.S, y=lam.S.f)) + geom_point(alpha=0.5) + 
  facet_grid(SDM~issue) +
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)

ggplot(out_IPM, aes(x=lam.U, y=lam.U.f)) + geom_point(alpha=0.5) + 
  facet_wrap(~issue) +
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)

ggplot(out_IPM, aes(x=D, y=D.f)) + geom_point(alpha=0.5) + 
  facet_wrap(~issue) +
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)

ggplot(out_IPM, aes(x=nSdStay, y=nSdStay.f)) + geom_point(alpha=0.5) + 
  facet_wrap(~issue) +
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)
