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
out <- read.csv(here(paste0("out/", sp, "out.csv")))

theme_set(theme_bw())
# barplots showing % correct for presence/absence; fill=SDM; facet=issues
ggplot(out, aes(fill=outcome, x=SDM)) + geom_bar(position="fill") + 
  facet_wrap(~s_Iss) + scale_fill_brewer(type="div")
ggplot(out_IPM, aes(fill=sign(log(lambda.f))==sign(log(lambda)), x=s_Iss)) + 
  geom_bar(position="fill")
ggplot(out_IPM, aes(fill=sign(log(lam.U.f))==sign(log(lam.U)), x=s_Iss)) + 
  geom_bar(position="fill")
ggplot(out, aes(fill=sign(log(lam.S.f))==sign(log(lam.S)), x=SDM)) + 
  geom_bar(position="fill") + facet_wrap(~s_Iss)

ggplot(out, aes(x=lon, y=lat, fill=Surv.S.f-Surv.S)) +
  geom_tile() + facet_grid(SDM~s_Iss) + scale_fill_gradient2()
ggplot(out, aes(x=lon, y=lat, fill=(Surv.S.f-Surv.S)/Surv.S)) +
  geom_tile() + facet_grid(SDM~s_Iss) + scale_fill_gradient()

ggplot(out, aes(x=lon, y=lat, fill=outcome)) +
  geom_tile() + facet_grid(SDM~s_Iss) + scale_fill_brewer(type="div")
ggplot(out, aes(x=lon, y=lat, fill=Surv.S.f)) +
  geom_tile() + facet_grid(SDM~s_Iss)
ggplot(out, aes(x=lon, y=lat, fill=Surv.S.f > 1)) +
  geom_tile() + facet_grid(SDM~s_Iss)
ggplot(out, aes(x=lon, y=lat, fill=log(lam.S.f))) +
  geom_tile() + facet_grid(SDM~s_Iss) + scale_fill_gradient2()
ggplot(out, aes(x=lon, y=lat, fill=log(lam.S.f)>=0)) +
  geom_tile() + facet_grid(SDM~s_Iss)

ggplot(out, aes(x=lon, y=lat, fill=Surv.S.f)) + geom_tile() + 
  facet_grid(SDM~s_Iss) + scale_fill_gradient(low="white", high="red")

ggplot(out, aes(x=lon, y=lat, fill=sign(log(lam.S.f+.001))==sign(log(lam.S)))) +
  geom_tile() + facet_grid(SDM~s_Iss)

ggplot(out_IPM, aes(x=lon, y=lat, fill=sign(log(lambda.f))==sign(log(lambda)))) +
  geom_tile() + facet_wrap(~s_Iss)

ggplot(out_IPM, aes(x=lon, y=lat, fill=N.U.f-N.U)) +
  geom_tile() + facet_wrap(~s_Iss) + scale_fill_gradient2()

ggplot(out_IPM, aes(x=lon, y=lat, fill=lam.U.f-lam.U)) +
  geom_tile() + facet_wrap(~s_Iss) + scale_fill_gradient2()

ggplot(out_IPM, aes(x=lon, y=lat, fill=D.f-D)) +
  geom_tile() + facet_wrap(~s_Iss) + scale_fill_gradient2()

ggplot(out_IPM, aes(x=lon, y=lat, fill=Btmax.f-Btmax)) +
  geom_tile() + facet_wrap(~s_Iss) + scale_fill_gradient2()

ggplot(out, aes(x=lon, y=lat, fill=nSeed.f-nSeed)) +
  geom_tile() + facet_grid(SDM~s_Iss) + scale_fill_gradient2()

ggplot(out, aes(x=lon, y=lat, fill=Rcr.S.f-Rcr.S)) +
  geom_tile() + facet_grid(SDM~s_Iss) + scale_fill_gradient2()

ggplot(out_IPM, aes(x=lon, y=lat, fill=N.S.f-N.S)) +
  geom_tile() + facet_wrap(~s_Iss) + scale_fill_gradient2()

ggplot(out, aes(x=N.S, y=N.S.f, colour=SDM, group=SDM)) + 
  geom_point(alpha=0.5) + facet_wrap(~s_Iss) +
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)

ggplot(out, aes(x=Surv.S, y=Surv.S.f, colour=SDM, group=SDM)) + 
  geom_point(alpha=0.5) + facet_wrap(~s_Iss) +
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)

ggplot(out, aes(x=Rcr.S, y=Rcr.S.f, colour=SDM, group=SDM)) + 
  geom_point(alpha=0.5) + facet_wrap(~s_Iss) +
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)
ggplot(out, aes(x=Rcr.S.f-Rcr.S, colour=s_Iss)) + geom_density()

ggplot(out_IPM, aes(x=lambda, y=lambda.f)) + geom_point(alpha=0.5) + 
  facet_wrap(~s_Iss) +
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)

ggplot(out_IPM, aes(x=lam.S, y=lam.S.f)) + geom_point(alpha=0.5) + 
  facet_wrap(~s_Iss) +
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)

ggplot(out_IPM, aes(x=lam.U, y=lam.U.f)) + geom_point(alpha=0.5) + 
  facet_wrap(~s_Iss) +
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)

ggplot(out_IPM, aes(x=D, y=D.f)) + geom_point(alpha=0.5) + 
  facet_wrap(~s_Iss) +
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)

ggplot(out_IPM, aes(x=nSdStay, y=nSdStay.f)) + geom_point(alpha=0.5) + 
  facet_wrap(~s_Iss) +
  stat_smooth(se=F, method="loess") + geom_abline(slope=1, size=1)
