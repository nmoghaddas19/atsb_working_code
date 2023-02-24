library(tidyverse)
setwd("Documents/GitHub/")

mali_jul_sep20 <- read.csv("211110 atsb_entomo_mali_hlc_dataset_jul2020_to_oct2021/From Jul-20 to Sep-2020-Table 1.csv")
mali_oct_jun21 <- read.csv("211110 atsb_entomo_mali_hlc_dataset_jul2020_to_oct2021/From Oct-20 to Jun-21-Table 1.csv")
mali_jul_oct21 <- read.csv("211110 atsb_entomo_mali_hlc_dataset_jul2020_to_oct2021/From Jul-21 to Oct-21-Table 1.csv")
mali_jul_oct21 <- mali_jul_oct21[,-19]
mali_data <- rbind(mali_jul_sep20, mali_oct_jun21, mali_jul_oct21)
rm("mali_jul_oct21","mali_jul_sep20","mali_oct_jun21")

mali_data$TIME_FRAME_CODE[mali_data$TIME_FRAME_CODE==7 | mali_data$TIME_FRAME_CODE==8] <- 1
mali_data$TIME_FRAME_CODE[mali_data$TIME_FRAME_CODE==9 | mali_data$TIME_FRAME_CODE==10] <- 2
mali_data$TIME_FRAME_CODE[mali_data$TIME_FRAME_CODE==11 | mali_data$TIME_FRAME_CODE==12] <- 3
mali_data$TIME_FRAME_CODE[mali_data$TIME_FRAME_CODE==13 | mali_data$TIME_FRAME_CODE==14] <- 4
mali_data$TIME_FRAME_CODE[mali_data$TIME_FRAME_CODE==15 | mali_data$TIME_FRAME_CODE==16] <- 5
mali_data$TIME_FRAME_CODE[mali_data$TIME_FRAME_CODE==17 | mali_data$TIME_FRAME_CODE==18] <- 6

mali_data <- mali_data[-which(is.na(mali_data$TIME_FRAME_CODE)),]

mali_data %>%
  group_by(date_Ident, TIME_FRAME_CODE) %>%
  summarise(Daily_Total = sum(Count)) -> df_1
unique(mali_data$TIME_FRAME_CODE)

mali_data |>
  mutate(time_frame = ifelse((TIME_FRAME_CODE==7 | TIME_FRAME_CODE==8), 1, 
                            ifelse((TIME_FRAME_CODE==9 | TIME_FRAME_CODE==10), 2,
                                   ifelse((TIME_FRAME_CODE==11 | TIME_FRAME_CODE==12), 3,
                                          ifelse((TIME_FRAME_CODE==13 | TIME_FRAME_CODE==14), 4,
                                                 ifelse((TIME_FRAME_CODE==15 | TIME_FRAME_CODE==16), 5,
                                                        ifelse((TIME_FRAME_CODE==17 | TIME_FRAME_CODE==18), 6, 
                                                               TIME_FRAME_CODE))))))) |>
  group_by(date_Ident, time_frame) |>
  summarise(Count=n()) -> df_1
|> ungroup() |>
  complete(date_Ident, time_frame, fill=list(count = 0)) -> test 

df_1 %>%
  group_by(date_Ident, time_frame) %>%
  summarise(Hourly_Count = sum(Count)) %>%
  group_by(date_Ident) %>%
  summarise(Days_with_counts = sum(Hourly_Count >= 0))

df_1 %>%
  group_by(date_Ident, time_frame) %>%
  summarise(Hourly_Count = sum(Count))

df_1 %>%
  group_by(date_Ident, time_frame) %>%
  summarise(Hourly_Count = sum(Count)) %>%
  ungroup() %>%
  complete(date_Ident, time_frame, fill = list(Hourly_Count = 0)) |>
  mutate(Date=parse_date_time(completed.df_1$date_Ident, c("mdy", "dmy"))) |>
  arrange(Date) |>
  mutate(Month_Yr = format_ISO8601(Date, precision = "ym")) -> completed.df_1

graph <- function(df) {
  plot <- plot(df$time_frame, df$Hourly_Count, frame.plot = F, ylim = c(0,700),
       pch=20, xlab = "Time Frame Code", ylab = "Count")
  
  df |> 
    group_by(time_frame) |>
    summarise(mean=mean(Hourly_Count), sd=sd(Hourly_Count)) -> out
  
  lines(out$time_frame, out$mean, col="red", lwd=2)
  arrows(out$time_frame, 
         out$mean-out$sd,
         out$time_frame,
         out$mean+out$sd,
         angle = 90,
         code = 3,
         length = 0.1, 
         col="red")
}
graph(completed.df_1)

months_ <- c("2020-07","2020-08", "2020-09", "2020-10", "2020-11", "2020-12",
             "2021-01", "2021-02", "2021-03", "2021-04", "2021-05", "2021-06",
             "2021-07", "2021-08", "2021-09", "2021-10","2021-11")
plot.new()
par(mfrow=c(3,3))
for (i in 1:length(months_)) {
  completed.df_1 |>
    filter(grepl(months_[i], Date)) |>
    graph()
  title(months_[i])
}

hist(completed.df_1$Hourly_Count, breaks=40, xlab = "Count")

completed.df_1 |>
  group_by(Date) |>
  summarise(Daily_Count = sum(Hourly_Count)) |> 
  mutate(Month_Yr = as.factor(format_ISO8601(Date, precision = "ym"))) -> Daily_Counts

plot.new()
par(mfrow=c(1,1))
hist(Daily_Counts$Daily_Count, breaks=50, xlab = "Daily Count")

plot.default(Daily_Counts$Month_Yr, Daily_Counts$Daily_Count, pch=20, cex=0.9, frame.plot = F,
             xaxt = "n", xlab = "Month", ylab = "Daily Count")
axis(1, at=1:17, labels=months_)
Daily_Counts |> 
  group_by(Month_Yr) |>
  summarise(mean=mean(Daily_Count), sd=sd(Daily_Count)) -> outt
lines(1:17, outt$mean, col="red", lwd=2)
arrows(1:17, 
       outt$mean-outt$sd,
       1:17,
       outt$mean+outt$sd,
       angle = 90,
       code = 3,
       length = 0.1, 
       col="red")
title("Daily Counts by Month")
