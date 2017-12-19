setwd("C:/Users/Valentin/Documents/MRes/Data")
source('data_extraction.R')
source('plot_data.R')
source('exploration_data.R')
source('fitting_data.R')


### Example:
## STL
plot_stl(data_leeds_1994_1998$O3, name = "O3", frequ = 365)
plot_stl(data_leeds_1994_1998$NO2x, name = "NO2x", frequ = 365)
plot_stl(data_leeds_1994_1998$SO2x, name = "SO2x", frequ = 182)
plot_stl(data_leeds_1994_1998$NO, name = "NO", frequ = 365)
plot_stl(data_leeds_1994_1998$PM10, name = "PM10", frequ = 90)

## Thresholds:
prob.quantile = 0.9
# O3:
eme_o3 <- empirical_mean_excess(data_leeds$O3,T)
threshold_o3 <- empirical_threshold_estimate(eme_o3, T)
q_o3 <- quantile(data_leeds_1994_1998$O3,prob.quantile,na.rm = T)

# NO2:
eme_no2x <- empirical_mean_excess(data_leeds_1994_1998$NO2x,T)
threshold_no2x <- empirical_threshold_estimate(eme_no2x, T)
q_no2x <- quantile(data_leeds_1994_1998$NO2x,prob.quantile,na.rm = T)

# SO2:
eme_so2x <- empirical_mean_excess(data_leeds_1994_1998$SO2x,T)
threshold_so2x <- empirical_threshold_estimate(eme_so2x, T)
q_so2x <- quantile(data_leeds_1994_1998$SO2x,prob.quantile,na.rm = T)

# NO:
eme_no <- empirical_median_excess(data_leeds_1994_1998$NO,T)
threshold_no <- empirical_threshold_estimate(eme_no, T)
q_no <- quantile(data_leeds_1994_1998$NO,prob.quantile,na.rm = T)

# PM10:
eme_pm10 <- empirical_median_excess(data_leeds_1994_1998$PM10,T)
threshold_pm10 <- empirical_threshold_estimate(eme_pm10, T)
q_pm10 <- quantile(data_leeds_1994_1998$PM10,prob.quantile,na.rm = T)

marg.threshold <- c(threshold_o3,
                    threshold_no2x,
                    threshold_so2x,
                    threshold_no,
                    threshold_pm10)
marg.threshold
quan.threshold <- c(q_o3,
                    q_no2x,
                    q_so2x,
                    q_no,
                    q_pm10)
quan.threshold
### Histograms
plot_histograms(data_to_plot = data_leeds_1994_1998,
                marg.threshold = marg.threshold, brk_hist = 40, brk_exc = 10,
                exc_without_zero = T)
plot_histograms(data_to_plot = data_leeds_1994_1998,
                marg.threshold = quan.threshold, brk_hist = 40, brk_exc = 10,
                exc_without_zero = T)

### ACF/PACF
plot_acf_pacf(data_to_plot = data_leeds)
plot_acf_pacf_exc(data_to_plot = data_leeds, marg.threshold = marg.threshold)
plot_acf_pacf_exc(data_to_plot = data_leeds, marg.threshold = quan.threshold)
