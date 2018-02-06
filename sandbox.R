library(pvm)

data <- readr::read_rds("~/projects/pv-monitor/data/srsim_100000_500_500_1.000_5.000_1.000_5.000_50_0.500_50_0.500_100_0.250_1.RDS")

a <- data$tables$a 
b <- data$tables$b 
c <- data$tables$c 
d <- data$tables$d 

prior <- fitPriorParametersGPS(a, b, c, d)

x <- GPS(a, b, c, d, prior = prior)

hist(x)
x <- GPS(a, b, c, d, prior = prior, alpha = 0.05)
