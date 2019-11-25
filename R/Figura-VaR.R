library(ggplot2)
x <- seq(-3, 3, by = 0.05)
y <- dnorm(x)
df <- data.frame(x, y)
p <- ggplot(df, aes(x, y)) + 
  geom_area(data = subset(df, x < -1), fill = "tomato") +
  #geom_area(data = subset(df, x < -1.96), fill = "red") +
  #geom_area(data = subset(df, x > 1), fill = "lightgreen") +
  #geom_area(data = subset(df, x > 1.96), fill = "green") +
  geom_line()+
  labs(x = "", y = "")

p+annotate("text", x = c(-1, -1.5), y = c(-0.02, 0.05), label = c("VaR[alpha]", "~alpha"), parse = TRUE, size = 6)

p+annotate("text", x = c(-2, -2), y = c(-0.02, 0), label = c("ES[alpha]", "."), parse = TRUE, size = 6)

