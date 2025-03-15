packages <- c("depmixS4", "dplyr", "lubridate", "ggplot2")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(depmixS4)  
library(dplyr)     
library(lubridate) 
library(ggplot2)   

data_raw <- read.csv("data/TermProjectData.txt", header = TRUE, stringsAsFactors = FALSE)

data_raw$DateTime <- dmy_hms(paste(data_raw$Date, data_raw$Time))
datetimes <- as.integer(data_raw$DateTime)
data_raw$Epoch <- datetimes

replace_na_with_lerped <- function(ind, dep){
  new_values <- approx(ind, dep, ind[which(is.na(dep))])
  #print(new_values)
  return(new_values)
}

cdf <- data_raw

nv = replace_na_with_lerped(data_raw$Epoch, data_raw$Global_active_power)
cdf$Global_active_power[which(is.na(data_raw$Global_active_power))] <- nv$y
nv = replace_na_with_lerped(data_raw$Epoch, data_raw$Global_reactive_power)
cdf$Global_reactive_power[which(is.na(data_raw$Global_reactive_power))] <- nv$y
nv = replace_na_with_lerped(data_raw$Epoch, data_raw$Voltage)
cdf$Voltage[which(is.na(data_raw$Voltage))] <- nv$y
nv = replace_na_with_lerped(data_raw$Epoch, data_raw$Global_intensity)
cdf$Global_intensity[which(is.na(data_raw$Global_intensity))] <- nv$y
nv = replace_na_with_lerped(data_raw$Epoch, data_raw$Sub_metering_1)
cdf$Sub_metering_1[which(is.na(data_raw$Sub_metering_1))] <- nv$y
nv = replace_na_with_lerped(data_raw$Epoch, data_raw$Sub_metering_2)
cdf$Sub_metering_2[which(is.na(data_raw$Sub_metering_2))] <- nv$y
nv = replace_na_with_lerped(data_raw$Epoch, data_raw$Sub_metering_3)
cdf$Sub_metering_3[which(is.na(data_raw$Sub_metering_3))] <- nv$y


#### STANDARDIZATION
sdf <- cdf
sdf$Global_active_power <- scale(cdf$Global_active_power, center = TRUE, scale = TRUE)
sdf$Global_reactive_power <- scale(cdf$Global_reactive_power, center = TRUE, scale = TRUE)
sdf$Voltage <- scale(cdf$Voltage, center = TRUE, scale = TRUE)
sdf$Global_intensity <- scale(cdf$Global_intensity, center = TRUE, scale = TRUE)
sdf$Sub_metering_1 <- scale(cdf$Sub_metering_1, center = TRUE, scale = TRUE)
sdf$Sub_metering_2 <- scale(cdf$Sub_metering_2, center = TRUE, scale = TRUE)
sdf$Sub_metering_3 <- scale(cdf$Sub_metering_3, center = TRUE, scale = TRUE)

if (!dir.exists("data/processed")) {
  dir.create("data/processed", recursive = TRUE)
}
write.csv(sdf, "data/processed/TermProjectData_Standardized.csv", row.names = FALSE)
cat("Processed data saved to data/processed/TermProjectData_Standardized.csv\n")

write.table(sdf, "data/processed/TermProjectData_Standardized.txt", row.names = FALSE, append = FALSE, sep = ",")
cat("Processed data saved to data/processed/TermProjectData_Standardized.txt\n")