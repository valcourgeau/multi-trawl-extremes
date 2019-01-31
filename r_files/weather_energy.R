
merged_dataset_folder <- "C:/Users/Valentin/Documents/GitHub/multi-trawl-extremes/data/merged-datasets/"
setwd(merged_dataset_folder)
load("energy-weather-merge.Rda")

# cheching NA proportions
apply(energy_weather_merged, function(x){(is.na(x) %>% sum) / length(x)}, MARGIN = 2)

library(ggplot2)
library(tabplot)
energy_weather_merged$wind_direction.Houston

tp.houston <- tableplot(energy_weather_merged, 
          select = c(humidity.Houston, 
                     temperature.Houston,
                     pressure.Houston,
                     wind_direction.Houston,
                     wind_speed.Houston,
                     weather_description.Houston,
                     index.y),
          sortCol = index.y)

tp.new.york <- tableplot(energy_weather_merged, 
          select = c(humidity.New.York, 
                     temperature.New.York,
                     pressure.New.York,
                     wind_direction.New.York,
                     wind_speed.New.York,
                     weather_description.New.York,
                     index.y),
          sortCol = index.y)

par(mfrow=c(1, 1))

tp.energy <- tableplot(energy_weather_merged[,219:229], 
                       sortCol = index.y)

#corr_matrix <- cor(energy_weather_merged[sapply(energy_weather_merged, is.numeric)], use = "pairwise.complete.obs")
library(corrplot)
#corrplot(corr_matrix, method="color", type="upper", order="hclust")
