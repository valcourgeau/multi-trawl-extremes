library(rlist)
library(VineCopula)
library(magrittr)

setwd('~/GitHub/multi-trawl-extremes/r_files/')

# AIR POLLUTION
tron_ap <- rlist::list.load('~/GitHub/multi-trawl-extremes/r_files/2019-5-4-10-48-17_tron.RData')
vines_ap <- rlist::list.load('~/GitHub/multi-trawl-extremes/r_files/2019-5-4-10-48-17_vines.RData')
data_ap <- read.table('air_pollution_data')
data_ap_matrix <- rlist::list.load('air_pollution_rerun_3_matrix.RData')


# ENERGY-WEATHER
tron_ap <- rlist::list.load('~/GitHub/multi-trawl-extremes/data/merged-datasets/tron_east_ligh_1_to_72_2nd_tron.RData')
vines_ap <- rlist::list.load('~/GitHub/multi-trawl-extremes/data/merged-datasets/vine_east_light_1_to_72_2nd_vines.RData')
data_ap <- read.table('')
data_ap_matrix <- rlist::list.load('~/GitHub/multi-trawl-extremes/data/merged-datasets/matrix_east_light_1_to_72_2nd_matrix.RData')


par(mfrow=c(2,3), mar=c(1.1,0.2,0.2,0.5))
plot(vines_ap[[12]][[1]], legend.pos="bottomright", type=1,
     interactive=F, label.bg="white", label.col="black", lwd=2, label.cex = 2.2, cex.main = 2,  edge.lwd=1.1,
     edge.labels=c("family-par"), edge.len = 8,
     edge.label.cex=2.2, edge.label.col="blue", label.lwd=2)

