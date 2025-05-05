set.seed(2)
n <- 10000
counts <- sample(c(0, 1), n, replace = TRUE, prob = c(0.999, 0.001))
normalized <- ifelse(counts == 0, 0, rnorm(sum(counts), mean = 0.03, sd = 0.01))

df <- data.frame(
  value = c(counts, normalized),
  type = rep(c("counts", "normalized"), each = n)
)

# ggplot(df, aes(x = type, y = value)) +
#   geom_jitter(width = 0.2, height = 0.03, alpha = 0.5) +
#   theme_minimal()
#
#
# library(ggpubr)
# library(ggplot2)
#
# # ggboxplot (原始錯誤會發生的版本)
# p1 <- ggboxplot(df,
#                 x = "type", y = "value",
#                 add = "jitter",
#                 color = "type", palette = "jco")
#
# p1
# # ggplot 自己畫
# p2 <- ggplot(df, aes(x = type, y = value, color = type)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(width = 0.2, height = 0.03, alpha = 0.5) +
#   scale_color_brewer(palette = "Set1") +
#   theme_minimal()
#
# # 比較看看 jitter 方式是否一致
# p1; p2
#
#
# ################################################################################




library(ggpubr)
library(ggplot2)

# 將 counts 與 normalized 拆開
df_counts <- df %>% filter(type == "counts")
df_norm   <- df %>% filter(type == "normalized")

# 畫 counts 的 jitter（你原本錯誤現象的來源）
p_counts <- ggboxplot(df_counts,
                      x = "type", y = "value",
                      add = "jitter",
                      color = "type", palette = "jco") +
  ggtitle("ggboxplot - Counts") +
  theme_minimal()

# 畫 normalized 的 jitter（通常沒偏移）
p_norm <- ggboxplot(df_norm,
                    x = "type", y = "value",
                    add = "jitter",
                    color = "type", palette = "jco") +
  ggtitle("ggboxplot - Normalized") +
  theme_minimal()

# 顯示
print(p_counts)
print(p_norm)
print(p_counts+p_norm)
