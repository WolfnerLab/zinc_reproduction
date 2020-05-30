library(tidyverse)
library(gridExtra)
library(grid)
library(ggsignif)
library(ggplotify)
library(agricolae)

TPEN_data <- read_csv("TPEN_data.csv")
TPEN_F2_data <- read_csv("TPEN_F2_data.csv")
TPEN_data <- TPEN_data %>%
    mutate(Food_TPEN = as.factor(Food_TPEN), Male_TPEN = as.factor(Male_TPEN), Female_TPEN = as.factor(Female_TPEN),
           Male_Zn = as.factor(Male_Zn), Female_Zn = as.factor(Female_Zn))
TPEN_F2_data <- TPEN_F2_data %>%
    mutate(Food_TPEN = as.factor(Food_TPEN), Male_TPEN = as.factor(Male_TPEN), Female_TPEN = as.factor(Female_TPEN),
           Male_Zn = as.factor(Male_Zn), Female_Zn = as.factor(Female_Zn))

TPEN_data_50 <- filter(TPEN_data, Food_TPEN == 50)
TPEN_data_100 <- filter(TPEN_data, Food_TPEN == 100)

tx50 <- with(TPEN_data_50, interaction(Male_TPEN, Male_Zn, Female_TPEN, Female_Zn))
tx100 <- with(TPEN_data_100, interaction(Male_TPEN, Male_Zn, Female_TPEN, Female_Zn))

# Summary plots
# Figure 1A
fig1a <- ggplot(TPEN_data_50, aes(x = interaction(Male_TPEN, Male_Zn, Female_TPEN, Female_Zn), y = Mean_hatch)) +
    geom_boxplot(aes(fill = Food_TPEN), outlier.colour = NA) +
    geom_jitter(width = 0.15) +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 30, margin = margin(t = 0, r = 20, b = 0, l = 0)),
          axis.text.x = element_blank(), axis.text.y = element_text(size = 20), axis.ticks.x = element_blank(),
          panel.grid = element_blank(), legend.position = "none") +
    scale_x_discrete(limits = c("0.0.0.0", "0.0.50.0", "50.0.0.0", "50.0.50.0",
                                "0.0.50.50", "50.50.0.0", "50.50.50.50")) +
    scale_fill_manual(values = "#DDDDDD") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(y = "Average hatchability") +
    stat_summary(geom = "text", label = c("a", "a", "a", "b", "b", "a", "a"),
                 size = 8, fun.y = max, vjust = -0.3)

hatch_50_aov <- aov(Mean_hatch ~ tx50, data = TPEN_data_50)
HSD.test(hatch_50_aov, "tx50", console = TRUE)

# Figure 1B
fig1b <- ggplot(TPEN_data_50, aes(x = interaction(Male_TPEN, Male_Zn, Female_TPEN, Female_Zn), y = Mean_egg)) +
    geom_boxplot(aes(fill = Food_TPEN), outlier.colour = NA) +
    geom_jitter(width = 0.15) +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 30, margin = margin(t = 0, r = 20, b = 0, l = 0)),
          axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
          panel.grid = element_blank(), legend.position = "none") +
    scale_x_discrete(limits = c("0.0.0.0", "0.0.50.0", "50.0.0.0", "50.0.50.0",
                                "0.0.50.50", "50.50.0.0", "50.50.50.50"),
                     labels = c("0\n0\n0\n0\n14", "0\n50\n0\n0\n15", "50\n0\n0\n0\n14", "50\n50\n0\n0\n14",
                                "0\n50\n0\n50\n9", "50\n0\n50\n0\n10", "50\n50\n50\n50\n10")) +
    scale_fill_manual(values = "#DDDDDD") +
    labs(y = "Average egg number") +
    stat_summary(geom = "text", label = c("a", "ab", "ab", "a", "a", "ab", "b"),
                 size = 8, fun.y = max, vjust = -0.3)

egg_number_50_aov <- aov(Mean_egg ~ tx50, data = TPEN_data_50)
HSD.test(egg_number_50_aov, "tx50", console = TRUE)

grid.draw(rbind(ggplotGrob(fig1a), ggplotGrob(fig1b), size = "last"))

# Figure 1C
fig1c <- ggplot(TPEN_data_100, aes(x = interaction(Male_TPEN, Male_Zn, Female_TPEN, Female_Zn), y = Mean_hatch)) +
    geom_boxplot(aes(fill = Food_TPEN), outlier.colour = NA) +
    geom_jitter(width = 0.15) +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 30, margin = margin(t = 0, r = 20, b = 0, l = 0)),
          axis.text.x = element_blank(), axis.text.y = element_text(size = 20), axis.ticks.x = element_blank(),
          panel.grid = element_blank(), legend.position = "none") +
    scale_x_discrete(limits = c("0.0.0.0", "0.0.100.0", "100.0.0.0", "100.0.100.0",
                                "0.0.100.100", "100.100.0.0", "100.100.100.100")) +
    scale_fill_manual(values = "#555555") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(y = "Average hatchability") +
    stat_summary(geom = "text", label = c("a", "a", "a", "b", "b", "a", "a"),
                 size = 8, fun.y = max, vjust = -0.3)

hatch_100_aov <- aov(Mean_hatch ~ tx100, data = TPEN_data_100)
HSD.test(hatch_100_aov, "tx100", console = TRUE)

# Figure 1D
fig1d <- ggplot(TPEN_data_100, aes(x = interaction(Male_TPEN, Male_Zn, Female_TPEN, Female_Zn), y = Mean_egg)) +
    geom_boxplot(aes(fill = Food_TPEN), outlier.colour = NA) +
    geom_jitter(width = 0.15) +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 30, margin = margin(t = 0, r = 20, b = 0, l = 0)),
          axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
          panel.grid = element_blank(), legend.position = "none") +
    scale_x_discrete(limits = c("0.0.0.0", "0.0.100.0", "100.0.0.0", "100.0.100.0",
                                "0.0.100.100", "100.100.0.0", "100.100.100.100"),
                     labels = c("0\n0\n0\n0\n14", "0\n100\n0\n0\n14", "100\n0\n0\n0\n15", "100\n100\n0\n0\n15",
                                "0\n100\n0\n100\n9", "100\n0\n100\n0\n10", "100\n100\n100\n100\n10")) +
    scale_fill_manual(values = "#555555") +
    labs(y = "Average egg number") +
    stat_summary(geom = "text", label = c("ab", "b", "a", "c", "c", "a", "a"),
                 size = 8, fun.y = max, vjust = -0.3)

egg_number_100_aov <- aov(Mean_egg ~ tx100, data = TPEN_data_100)
HSD.test(egg_number_100_aov, "tx100", console = TRUE)

grid.draw(rbind(ggplotGrob(fig1c), ggplotGrob(fig1d), size = "last"))

# Figure S1
figs1a1 <- ggplot(TPEN_data_50, aes(x = interaction(Male_TPEN, Male_Zn, Female_TPEN, Female_Zn), y = Day1_hatch)) +
    geom_boxplot(aes(fill = Food_TPEN), outlier.colour = NA) +
    geom_jitter(width = 0.15) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 20),
          axis.text.x = element_blank(), axis.text.y = element_text(size = 20), axis.ticks.x = element_blank(),
          panel.grid = element_blank(), legend.position = "none") +
    scale_x_discrete(limits = c("0.0.0.0", "0.0.50.0", "50.0.0.0", "50.0.50.0",
                                "0.0.50.50", "50.50.0.0", "50.50.50.50")) +
    scale_y_continuous(limits = c(0, 1.1), labels = scales::percent_format(accuracy = 1)) +
    scale_fill_manual(values = "#DDDDDD") +
    labs(y = "Day 1 hatchability") +
    stat_summary(geom = "text", label = c("a", "a", "a", "a", "a", "a", "a"),
                 size = 8, fun.y = max, vjust = -0.3)
figs1a2 <- ggplot(TPEN_data_50, aes(x = interaction(Male_TPEN, Male_Zn, Female_TPEN, Female_Zn), y = Day2_hatch)) +
    geom_boxplot(aes(fill = Food_TPEN), outlier.colour = NA) +
    geom_jitter(width = 0.15) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 20),
          axis.text.x = element_blank(), axis.text.y = element_text(size = 20), axis.ticks.x = element_blank(),
          panel.grid = element_blank(), legend.position = "none") +
    scale_x_discrete(limits = c("0.0.0.0", "0.0.50.0", "50.0.0.0", "50.0.50.0",
                                "0.0.50.50", "50.50.0.0", "50.50.50.50")) +
    scale_y_continuous(limits = c(0, 1.1), labels = scales::percent_format(accuracy = 1)) +
    scale_fill_manual(values = "#DDDDDD") +
    labs(y = "Day 2 hatchability") +
    stat_summary(geom = "text", label = c("a", "a", "a", "b", "b", "a", "a"),
                 size = 8, fun.y = max, vjust = -0.3)
figs1a3 <- ggplot(TPEN_data_50, aes(x = interaction(Male_TPEN, Male_Zn, Female_TPEN, Female_Zn), y = Day3_hatch)) +
    geom_boxplot(aes(fill = Food_TPEN), outlier.colour = NA) +
    geom_jitter(width = 0.15) +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 20),
          axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
          panel.grid = element_blank(), legend.position = "none") +
    scale_x_discrete(limits = c("0.0.0.0", "0.0.50.0", "50.0.0.0", "50.0.50.0",
                                "0.0.50.50", "50.50.0.0", "50.50.50.50"),
                     labels = c("0\n0\n0\n0\n14", "0\n50\n0\n0\n15", "50\n0\n0\n0\n14", "50\n50\n0\n0\n14",
                                "0\n50\n0\n50\n9", "50\n0\n50\n0\n10", "50\n50\n50\n50\n10")) +
    scale_y_continuous(limits = c(0, 1.1), labels = scales::percent_format(accuracy = 1)) +
    scale_fill_manual(values = "#DDDDDD") +
    labs(y = "Day 3 hatchability") +
    stat_summary(geom = "text", label = c("a", "a", "a", "b", "b", "a", "a"),
                 size = 8, fun.y = max, vjust = -0.3)
grid.draw(rbind(ggplotGrob(figs1a1), ggplotGrob(figs1a2), ggplotGrob(figs1a3), size = "last"))

figs1b1 <- ggplot(TPEN_F2_data, aes(x = interaction(Male_TPEN, Male_Zn, Female_TPEN, Female_Zn), y = Day1_hatch)) +
    geom_boxplot(aes(fill = Food_TPEN), outlier.colour = NA) +
    geom_jitter(width = 0.15) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 20),
          axis.text.x = element_blank(), axis.text.y = element_text(size = 20), axis.ticks.x = element_blank(),
          panel.grid = element_blank(), legend.position = "none") +
    scale_x_discrete(limits = c("0.0.0.0", "0.0.50.0", "0.0.50.50")) +
    scale_y_continuous(limits = c(0, 1.1), labels = scales::percent_format(accuracy = 1)) +
    scale_fill_manual(values = "#DDDDDD") +
    labs(y = "Day 1 hatchability") +
    stat_summary(geom = "text", label = c("a", "a", "a"),
                 size = 8, fun.y = max, vjust = -0.3)
figs1b2 <- ggplot(TPEN_F2_data, aes(x = interaction(Male_TPEN, Male_Zn, Female_TPEN, Female_Zn), y = Day2_hatch)) +
    geom_boxplot(aes(fill = Food_TPEN), outlier.colour = NA) +
    geom_jitter(width = 0.15) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 20),
          axis.text.x = element_blank(), axis.text.y = element_text(size = 20), axis.ticks.x = element_blank(),
          panel.grid = element_blank(), legend.position = "none") +
    scale_x_discrete(limits = c("0.0.0.0", "0.0.50.0", "0.0.50.50")) +
    scale_y_continuous(limits = c(0, 1.1), labels = scales::percent_format(accuracy = 1)) +
    scale_fill_manual(values = "#DDDDDD") +
    labs(y = "Day 2 hatchability") +
    stat_summary(geom = "text", label = c("a", "b", "a"),
                 size = 8, fun.y = max, vjust = -0.3)
figs1b3 <- ggplot(TPEN_F2_data, aes(x = interaction(Male_TPEN, Male_Zn, Female_TPEN, Female_Zn), y = Day3_hatch)) +
    geom_boxplot(aes(fill = Food_TPEN), outlier.colour = NA) +
    geom_jitter(width = 0.15) +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 20),
          axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
          panel.grid = element_blank(), legend.position = "none") +
    scale_x_discrete(limits = c("0.0.0.0", "0.0.50.0", "0.0.50.50"),
                     labels = c("0\n0\n0\n0\n8", "0\n50\n0\n0\n10", "0\n50\n0\n50\n10")) +
    scale_y_continuous(limits = c(0, 1.1), labels = scales::percent_format(accuracy = 1)) +
    scale_fill_manual(values = "#DDDDDD") +
    labs(y = "Day 3 hatchability") +
    stat_summary(geom = "text", label = c("a", "b", "a"),
                 size = 8, fun.y = max, vjust = -0.3)
grid.draw(rbind(ggplotGrob(figs1b1), ggplotGrob(figs1b2), ggplotGrob(figs1b3), size = "last"))
