# Normality Testing

library(ggplot2)
library(readxl)

data <- c(0.042, 0.005, 0.117, 0.244, 0.16, 0.013, 0.041, 0.094, 0.161, 0.107, 0.673, 0.3, 
          1.082, 0.022, 0.149, 0.002, 0.109, 0.275, 0.062, 0.22, 0.698, 1.431, 0.351, 0.088, 
          0.034, 0.006, 0.107, 0.471, 0.029, 0.133, 0.01, 0.037, 0.026, 0.012, 0.11, 0.004, 
          0.012, 0.067, 0.08, 0.036, 0.024, 0.116, 0.005, 0.002, 0.081, 0.047, 0.412, 0.031, 
          0.002, 0.033, 0.027, 0.01, 0.024, 0.01, 0.002, 0.474, 0.002, 0.009, 0.799, 0.254, 
          0.156, 1.464, 0.002, 0.095, 0.064, 0.268, 0.178, 0.063, 0.056, 0.16, 0.162, 0.242, 
          0.885, 0.03, 0.089, 0.033, 0.003, 1.18, 0.505, 0.327, 0.133, 2.329, 0.012, 0.226, 
          0.002, 0.007, 0.133, 0.009, 1.116, 0.373, 0.234, 0.013, 0.029, 0.019, 0.014, 0.065, 
          0.006, 1.133, 0.097, 0.014, 0.011, 0.001, 0.031, 0.175, 0.069, 0.903, 2.5, 0.007, 
          0.792, 0.252, 0.01, 0.573, 0.015, 0.026, 0.026, 0.004, 0.024, 0.302, 0.239, 0.125, 
          0.641, 0.017, 0.035, 0.153, 0.171, 0.158, 0.075, 0.056, 0.007, 0.028, 0.012, 0.014, 
          0.213, 0.007, 0.635, 0.08, 0.523, 0.928, 0.019, 0.112, 0.596, 0.193, 1.003, 0.815, 
          0.095, 0.047, 0.094, 0.02, 0.01, 0.009, 0.081, 0.001, 0.1, 0.173, 0.168, 0.164, 
          1.115, 0.138, 0.014, 0.028, 0.647, 0.112, 0.207, 0.112, 1.821, 0.012, 0.051, 0.015, 
          0.004, 0.158, 0.009, 0.405, 0.003, 0.072, 0.006, 0.049, 0.059, 2.89, 0.099, 
          0.583, 3.217, 0.016, 0.609, 0.086, 0.074, 0.09, 0.008, 0.155, 0.057, 0.028, 0.092, 
          0.018, 0.005, 0.277, 0.024, 0.084, 0.008, 0.09, 0.315, 0.42, 0.141, 1.02, 0.009, 
          0.26, 0.029, 0.202, 0.004, 0.006, 0.254, 0.058, 0.003, 0.03, 0.006, 0.232, 0.005, 
          0.015, 0.01, 0.094, 0.116, 0.02, 0.056, 0.021, 0.015, 0.007, 0.037, 0.148, 0.003, 
          0.002, 0.096, 0.047, 0.042, 0.011, 0.054, 0.038, 0.224, 0.026, 0.217, 0.002, 0.131, 
          0.022, 0.003, 0.053, 0.211, 0.053, 0.002, 0.122, 0.223, 0.539, 0.102, 0.103, 0.128, 
          0.554, 1.153, 0.224, 0.012, 0.002, 0.032, 0.007, 0.063, 0.006, 0.069, 0.035, 0.119, 
          0.258, 0.004, 0.106, 0.259, 0.005, 0.06, 0.13, 0.057, 0.041, 0.052, 0.003, 0.005, 
          0.061, 0.317, 0.022, 0.078, 0.021, 0.033, 0.072, 1.865, 0.023, 0.004, 0.028, 0.097, 
          0.212, 0.272, 0.001, 0.009, 0.556, 0.014, 0.998, 1.53, 0.611, 0.044, 0.471, 0.013, 
          0.015, 0.044, 0.275, 0.007, 0.008, 2.126, 0.106, 0.206, 0.05, 0.049, 0.376, 0.01, 
          0.04, 0.044, 0.002, 0.005, 0.013, 0.037, 0.895, 5.263)

res <- hist(log(data), breaks=40, plot=FALSE)
with(res, plot(counts ~ exp(mids), type="h", lwd=10, col="skyblue", log="x", 
               lend=2, xlab="Normalisation G to A (%)", ylab="Frequency"))
  box(lwd=2)
 log_data <- log(data)
  qqnorm(log_data) 
  qqline(log_data, col = "red")
  shapiro.test(log_data) 
  
  
  # Pearson Correlation
  library(readxl)
  library(ggplot2)
  library(dplyr)
  
  
  data <- read_excel("/Users/priscilla/Documents/Dry Lab/Merged_Gene_Data.xlsx")
  head(data)
  g_count <- data$g_count 
  normalized_g_to_a <- data$`Normalised_G_to_A_(%)`
  data_clean <- na.omit(data)
  max_value_row <- which.max(data_clean$`Normalised_G_to_A_(%)`)
  data_clean <- data_clean[-max_value_row, ]
  data_clean$g_count_log <- log(data_clean$g_count)
  correlation_result <- cor(data_clean$g_count_log, data_clean$`Normalised_G_to_A_(%)`)
  cat("Pearson correlation coefficient: ", correlation_result, "\n")
  model <- lm(`Normalised_G_to_A_(%)` ~ g_count_log, data = data_clean)
  summary(model)
  r_value <- round(correlation_result, 3) 
  p_value <- signif(summary(model)$coefficients[2, 4], 3)
  ggplot(data_clean, aes(x = g_count_log, y = `Normalised_G_to_A_(%)`)) +
    geom_point() +  
    geom_smooth(method = "lm", col = "blue") + 
    labs(title = "Correlation Between Log-Transformed Gene Count and Normalized G>A(%)",
         x = "Log G Residue Count",
         y = "Normalised G>A Mutations (%)") +
    annotate("text", x = max(data_clean$g_count_log) * 0.8, 
             y = max(data_clean$`Normalised_G_to_A_(%)`) * 0.9, 
             label = paste("R =", r_value, "\np =", p_value), 
             color = "black", size = 4, hjust = 0) +
    theme_minimal() +
    theme(
      axis.line = element_line(color = "black", size = 0.5),  
      axis.ticks = element_line(color = "black"),  
      axis.text = element_text(color = "black", size = 12), 
      axis.title = element_text(color = "black", size = 14),
      plot.title = element_text(size = 11, hjust = 0.5)
    )
  
  