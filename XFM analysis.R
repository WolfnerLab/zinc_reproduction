library(tidyverse)
library(ggsignif)
library(RColorBrewer)
library(car)
library(agricolae)
library(PMCMRplus)
library(grid)

norm_df <- read_csv("final_data.csv")
norm_df$stage <- factor(norm_df$stage, levels = c("stage1-8_oocyte", "stage9-13_oocyte", "stage14_oocyte", "in_vivo_egg", "in_vitro_egg"))
# Zn_conc unit: fmol/um^2
norm_df <- norm_df %>%
    mutate(Zn_conc = Zn_measure / 65.38 / roi_areas,
           Cu_conc = Cu_measure / 63.546 / roi_areas,
           Fe_conc = Fe_measure / 55.845 / roi_areas)

#All three transition metals in WT EGG CHAMBERS, unit: _measure(fmol), _conc(fmol/um^2)
transition_measure <- norm_df %>%
    filter(cal == "downstream", genotype == "WT", stage != "in_vitro_egg",
           focus == "cell", !(ID %in% c("2_0148", "2_0149", "2_0150"))) %>%
    select(ID, stage, Fe_measure, Cu_measure, Zn_measure) %>%
    mutate(Fe_measure = Fe_measure / 55.845, Cu_measure = Cu_measure / 63.546, Zn_measure = Zn_measure / 65.38) %>%
    gather(key = "metal", value = "measure", -ID, -stage)
transition_conc <- norm_df %>%
    filter(cal == "downstream", genotype == "WT", stage != "in_vitro_egg",
           focus == "cell", !(ID %in% c("2_0148", "2_0149", "2_0150"))) %>%
    select(ID, stage, Fe_conc, Cu_conc, Zn_conc) %>%
    gather(key = "metal", value = "conc", -ID, -stage)

eca <- ggplot(transition_measure, aes(x = metal, y = measure, fill = stage)) +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 30),
          axis.text.x = element_text(size = 25), panel.grid = element_blank(),
          legend.position = "none") +
    geom_boxplot(outlier.colour = NA) +
    geom_point(size = 0.7, position = position_dodge(width = 0.75)) +
    scale_x_discrete(labels = c("Cu", "Fe", "Zn")) +
    scale_fill_manual(name = "Stage",
                      labels = c("stage 1-8", "stage 9-13", "stage 14", "activated egg"),
                      values = brewer.pal(4, "PuBu")) +
    ylab("Egg chamber amount (fmol)") +
    stat_summary(geom = 'text', label = c("a", "a", "a", "a", "b", "b", "b", "a", "c", "d", "cd", "a"),
                 fun.y = max, vjust = -0.3, size = 8, position = position_dodge(width = 0.75))

ecc <- ggplot(transition_conc, aes(x = metal, y = conc, fill = stage)) +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 25),
          axis.text.x = element_text(size = 25), panel.grid = element_blank(),
          legend.title = element_text(size = 20), legend.text = element_text(size = 20)) +
    geom_boxplot(outlier.colour = NA) +
    geom_point(size = 0.7, position = position_dodge(width = 0.75)) +
    scale_x_discrete(labels = c("Cu", "Fe", "Zn")) +
    scale_fill_manual(name = "Stage",
                      labels = c("stage 1-8", "stage 9-13", "stage 14", "activated egg"),
                      values = brewer.pal(4, "PuBu")) +
    ylab(expression(paste("Egg chamber concentration (fmol/", "μm"^2, ")"))) +
    stat_summary(geom = 'text', label = c("a", "a", "a", "a", "b", "b", "b", "b", "e", "d", "c", "b"),
                 fun.y = max, vjust = -0.3, size = 8, position = position_dodge(width = 0.75))

grid.draw(cbind(ggplotGrob(eca), ggplotGrob(ecc), size = "last"))

#All three transition metals in WT OOCYTE only
transition_oocyte_measure <- norm_df %>%
    filter(genotype == "WT") %>%
    filter((stage == "stage9-13_oocyte" & focus == "focus") |
               (stage == "stage14_oocyte" & focus == "cell") |
               (stage == "in_vivo_egg" & focus == "cell")) %>%
    filter(!(ID %in% c("2_0148", "2_0149", "2_0150"))) %>%
    select(ID, stage, Fe_measure, Cu_measure, Zn_measure) %>%
    mutate(Fe_measure = Fe_measure / 55.845, Cu_measure = Cu_measure / 63.546, Zn_measure = Zn_measure / 65.38) %>%
    gather(key = "metal", value = "measure", -ID, -stage)    
transition_oocyte_conc <- norm_df %>%
    filter(genotype == "WT") %>%
    filter((stage == "stage9-13_oocyte" & focus == "focus") |
               (stage == "stage14_oocyte" & focus == "cell") |
               (stage == "in_vivo_egg" & focus == "cell")) %>%
    filter(!(ID %in% c("2_0148", "2_0149", "2_0150"))) %>%
    select(ID, stage, Fe_conc, Cu_conc, Zn_conc) %>%
    gather(key = "metal", value = "conc", -ID, -stage)

ooa <- ggplot(transition_oocyte_measure, aes(x = metal, y = measure, fill = stage)) +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 30),
          axis.text.x = element_text(size = 25), panel.grid = element_blank(),
          legend.position = "none") +
    geom_boxplot(outlier.colour = NA) +
    geom_point(size = 0.7, position = position_dodge(width = 0.75)) +
    scale_x_discrete(labels = c("Cu", "Fe", "Zn")) +
    scale_fill_manual(name = "Stage",
                      labels = c("stage 9-13", "stage 14", "activated egg"),
                      values = brewer.pal(4, "PuBu")[2:4]) +
    ylab("Oocyte amount (fmol)") +
    stat_summary(geom = 'text', label = c("a", "a", "a", "c", "c", "b", "f", "e", "d"),
                 fun.y = max, vjust = -0.3, size = 8, position = position_dodge(width = 0.75))

ooc <- ggplot(transition_oocyte_conc, aes(x = metal, y = conc, fill = stage)) +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 25),
          axis.text.x = element_text(size = 25), panel.grid = element_blank(),
          legend.title = element_text(size = 20), legend.text = element_text(size = 20)) +
    geom_boxplot(outlier.colour = NA) +
    geom_point(size = 0.7, position = position_dodge(width = 0.75)) +
    scale_x_discrete(labels = c("Cu", "Fe", "Zn")) +
    scale_fill_manual(name = "Stage",
                      labels = c("stage 9-13", "stage 14", "activated egg"),
                      values = brewer.pal(4, "PuBu")[2:4]) +
    ylab(expression(paste("Oocyte concentration (fmol/", "μm"^2, ")"))) + 
    stat_summary(geom = 'text', label = c("a", "a", "a", "b", "bc", "bc", "c", "e", "d"),
                 fun.y = max, vjust = -0.3, size = 8, position = position_dodge(width = 0.75))
grid.draw(cbind(ggplotGrob(ooa), ggplotGrob(ooc), size = "last"))

#All three transition metals in KO EGG CHAMBERS, unit: _measure(fmol), _conc(fmol/um^2)
KO_transition_measure <- norm_df %>%
    filter(cal == "downstream", genotype == "KO", stage != "in_vitro_egg",
           focus == "cell", !(ID %in% c("2_0148", "2_0149", "2_0150"))) %>%
    select(ID, stage, Fe_measure, Cu_measure, Zn_measure) %>%
    mutate(Fe_measure = Fe_measure / 55.845, Cu_measure = Cu_measure / 63.546, Zn_measure = Zn_measure / 65.38) %>%
    gather(key = "metal", value = "measure", -ID, -stage)
KO_transition_conc <- norm_df %>%
    filter(cal == "downstream", genotype == "KO", stage != "in_vitro_egg",
           focus == "cell", !(ID %in% c("2_0148", "2_0149", "2_0150"))) %>%
    select(ID, stage, Fe_conc, Cu_conc, Zn_conc) %>%
    gather(key = "metal", value = "conc", -ID, -stage)

txkotrans <- with(KO_transition_measure, interaction(metal, stage))
HSD.test(aov(measure ~ txkotrans, data = KO_transition_measure), "txkotrans", console = TRUE)
txkotransc <- with(KO_transition_conc, interaction(metal, stage))
HSD.test(aov(conc ~ txkotransc, data = KO_transition_conc), "txkotransc", console = TRUE)


koeca <- ggplot(KO_transition_measure, aes(x = metal, y = measure, fill = stage)) +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 30),
          axis.text.x = element_text(size = 25), panel.grid = element_blank(),
          legend.position = "none") +
    geom_boxplot(outlier.colour = NA) +
    geom_point(size = 0.7, position = position_dodge(width = 0.75)) +
    scale_x_discrete(labels = c("Cu", "Fe", "Zn")) +
    scale_fill_manual(name = "Stage",
                      labels = c("stage 1-8", "stage 9-13", "stage 14", "activated egg"),
                      values = brewer.pal(4, "PuBu")) +
    ylab("Egg chamber amount (fmol)") +
    stat_summary(geom = 'text', label = c("a", "a", "a", "a", "b", "b", "b", "a", "d", "d", "c", "a"),
                 fun.y = max, vjust = -0.3, size = 8, position = position_dodge(width = 0.75))

koecc <- ggplot(KO_transition_conc, aes(x = metal, y = conc, fill = stage)) +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 25),
          axis.text.x = element_text(size = 25), panel.grid = element_blank(),
          legend.title = element_text(size = 20), legend.text = element_text(size = 20)) +
    geom_boxplot(outlier.colour = NA) +
    geom_point(size = 0.7, position = position_dodge(width = 0.75)) +
    scale_x_discrete(labels = c("Cu", "Fe", "Zn")) +
    scale_fill_manual(name = "Stage",
                      labels = c("stage 1-8", "stage 9-13", "stage 14", "activated egg"),
                      values = brewer.pal(4, "PuBu")) +
    ylab(expression(paste("Egg chamber concentration (fmol/", "μm"^2, ")"))) +
    stat_summary(geom = 'text', label = c("a", "a", "a", "a", "cd", "d", "e", "b", "ef", "f", "g", "bc"),
                 fun.y = max, vjust = -0.3, size = 8, position = position_dodge(width = 0.75))

grid.draw(cbind(ggplotGrob(koeca), ggplotGrob(koecc), size = "last"))

#All three transition metals in WT OOCYTE only
KO_transition_oocyte_measure <- norm_df %>%
    filter(genotype == "KO") %>%
    filter((stage == "stage9-13_oocyte" & focus == "focus") |
               (stage == "stage14_oocyte" & focus == "cell") |
               (stage == "in_vivo_egg" & focus == "cell")) %>%
    filter(!(ID %in% c("2_0148", "2_0149", "2_0150"))) %>%
    select(ID, stage, Fe_measure, Cu_measure, Zn_measure) %>%
    mutate(Fe_measure = Fe_measure / 55.845, Cu_measure = Cu_measure / 63.546, Zn_measure = Zn_measure / 65.38) %>%
    gather(key = "metal", value = "measure", -ID, -stage)    
KO_transition_oocyte_conc <- norm_df %>%
    filter(genotype == "WT") %>%
    filter((stage == "stage9-13_oocyte" & focus == "focus") |
               (stage == "stage14_oocyte" & focus == "cell") |
               (stage == "in_vivo_egg" & focus == "cell")) %>%
    filter(!(ID %in% c("2_0148", "2_0149", "2_0150"))) %>%
    select(ID, stage, Fe_conc, Cu_conc, Zn_conc) %>%
    gather(key = "metal", value = "conc", -ID, -stage)

txkotransoo <- with(KO_transition_oocyte_measure, interaction(metal, stage))
HSD.test(aov(measure ~ txkotransoo, data = KO_transition_oocyte_measure), "txkotransoo", console = TRUE)
txkotransooc <- with(KO_transition_oocyte_conc, interaction(metal, stage))
HSD.test(aov(conc ~ txkotransooc, data = KO_transition_oocyte_conc), "txkotransooc", console = TRUE)

koooa <- ggplot(KO_transition_oocyte_measure, aes(x = metal, y = measure, fill = stage)) +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 30),
          axis.text.x = element_text(size = 25), panel.grid = element_blank(),
          legend.position = "none") +
    geom_boxplot(outlier.colour = NA) +
    geom_point(size = 0.7, position = position_dodge(width = 0.75)) +
    scale_x_discrete(labels = c("Cu", "Fe", "Zn")) +
    scale_fill_manual(name = "Stage",
                      labels = c("stage 9-13", "stage 14", "activated egg"),
                      values = brewer.pal(4, "PuBu")[2:4]) +
    ylab("Oocyte amount (fmol)") +
    stat_summary(geom = 'text', label = c("a", "a", "a", "c", "c", "b", "d", "d", "d"),
                 fun.y = max, vjust = -0.3, size = 8, position = position_dodge(width = 0.75))

koooc <- ggplot(KO_transition_oocyte_conc, aes(x = metal, y = conc, fill = stage)) +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 25),
          axis.text.x = element_text(size = 25), panel.grid = element_blank(),
          legend.title = element_text(size = 20), legend.text = element_text(size = 20)) +
    geom_boxplot(outlier.colour = NA) +
    geom_point(size = 0.7, position = position_dodge(width = 0.75)) +
    scale_x_discrete(labels = c("Cu", "Fe", "Zn")) +
    scale_fill_manual(name = "Stage",
                      labels = c("stage 9-13", "stage 14", "activated egg"),
                      values = brewer.pal(4, "PuBu")[2:4]) +
    ylab(expression(paste("Oocyte concentration (fmol/", "μm"^2, ")"))) + 
    stat_summary(geom = 'text', label = c("a", "a", "a", "b", "bc", "bc", "c", "e", "d"),
                 fun.y = max, vjust = -0.3, size = 8, position = position_dodge(width = 0.75))
grid.draw(cbind(ggplotGrob(koooa), ggplotGrob(koooc), size = "last"))


#In vitro activation data
IVA_data <- norm_df %>%
    filter(stage %in% c("stage14_oocyte", "in_vivo_egg", "in_vitro_egg"),
           genotype == "WT", !(ID %in% c("2_0148", "2_0149", "2_0150"))) %>%
    select(ID, stage, Zn_measure) %>%
    mutate(Zn_measure = Zn_measure / 65.38)
 ggplot(IVA_data, aes(x = stage, y = Zn_measure)) +
     theme_bw() +
     theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 30),
           axis.text.x = element_text(size = 25, angle = 90), panel.grid = element_blank(),
           legend.title = element_text(size = 20), legend.text = element_text(size = 20)) +
     geom_boxplot(outlier.colour = NA) +
     geom_point() +
     scale_x_discrete(labels = c("stage 14 oocyte", "in vivo activated egg", "in vitro activated egg")) +
     ylab("Zn amount (fmol)") +
     stat_summary(geom = 'text', label = c("a", "b", "a"),
                  fun.y = max, vjust = -0.3, size = 8, position = position_dodge(width = 0.75))

#Figure 2
WT_Zn_EC <- norm_df %>%
     filter(genotype == "WT", focus == "cell",
            !(ID %in% c("2_0148", "2_0149", "2_0150")),
            stage != "in_vitro_egg") %>%
     mutate(Zn_mol = Zn_measure / 65.38) %>%
     select(ID, stage, Zn_mol, Zn_conc)
wteca <- ggplot(WT_Zn_EC, aes(x = stage, y = Zn_mol, fill = stage)) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 20, hjust = 1, angle = 90), axis.title.x = element_blank(),
          axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 30), panel.grid = element_blank()) +
    geom_boxplot(outlier.colour = NA) +
    geom_jitter(width = 0.2) +
    scale_x_discrete(limits = c("stage1-8_oocyte", "stage9-13_oocyte", "stage14_oocyte", "in_vivo_egg"),
                     labels = c("stage 1-8", "stage 9-13", "stage 14", "activated egg")) +
    ylab("Egg chamber Zn amount (fmol)") +
    scale_fill_manual(values = brewer.pal(4, "PuBu")) +
    stat_summary(geom = "text", label = c("a", "bc", "b", "c"), size = 8, fun.y = max, vjust = -0.3)
HSD.test(aov(Zn_mol ~ stage, data = WT_Zn_EC), "stage", console = TRUE)

wtecc <- ggplot(WT_Zn_EC, aes(x = stage, y = Zn_conc, fill = stage)) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 20, hjust = 1, angle = 90), axis.title.x = element_blank(),
          axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 25), panel.grid = element_blank()) +
    geom_boxplot(outlier.colour = NA) +
    geom_jitter(width = 0.2) +
    scale_x_discrete(limits = c("stage1-8_oocyte", "stage9-13_oocyte", "stage14_oocyte", "in_vivo_egg"),
                     labels = c("stage 1-8", "stage 9-13", "stage 14", "activated egg")) +
    ylab(expression(paste("Egg chamber Zn concentration (fmol/", "μm"^2, ")"))) +
    scale_fill_manual(values = brewer.pal(4, "PuBu")) +
    stat_summary(geom = "text", label = c("a", "b", "c", "d"), size = 8, fun.y = max, vjust = -0.3)
HSD.test(aov(Zn_conc ~ stage, data = WT_Zn_EC), "stage", console = TRUE)

WT_Zn_OO <- norm_df %>%
    filter(genotype == "WT", !(ID %in% c("2_0148", "2_0149", "2_0150"))) %>%
    filter((stage == "stage9-13_oocyte" & focus == "focus") |
               (stage == "stage14_oocyte" & focus == "cell") |
               (stage == "in_vivo_egg" & focus == "cell")) %>%
    mutate(Zn_mol = Zn_measure / 65.38) %>%
    select(ID, stage, Zn_mol, Zn_conc)
wtooa <- ggplot(WT_Zn_OO, aes(x = stage, y = Zn_mol, fill = stage)) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 20, hjust = 1, angle = 90), axis.title.x = element_blank(),
          axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 30), panel.grid = element_blank()) +
    geom_boxplot(outlier.colour = NA) +
    geom_jitter(width = 0.2) +
    scale_x_discrete(limits = c("stage9-13_oocyte", "stage14_oocyte", "in_vivo_egg"),
                     labels = c("stage 9-13", "stage 14", "activated egg")) +
    ylab("Oocyte Zn amount (fmol)") +
    scale_fill_manual(values = brewer.pal(4, "PuBu")[2:4]) +
    stat_summary(geom = "text", label = c("a", "b", "c"), size = 8, fun.y = max, vjust = -0.3)
HSD.test(aov(Zn_mol ~ stage, data = WT_Zn_OO), "stage", console = TRUE)

wtooc <- ggplot(WT_Zn_OO, aes(x = stage, y = Zn_conc, fill = stage)) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 20, hjust = 1, angle = 90), axis.title.x = element_blank(),
          axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 25), panel.grid = element_blank()) +
    geom_boxplot(outlier.colour = NA) +
    geom_jitter(width = 0.2) +
    scale_x_discrete(limits = c("stage9-13_oocyte", "stage14_oocyte", "in_vivo_egg"),
                     labels = c("stage 9-13", "stage 14", "activated egg")) +
    ylab(expression(paste("Oocyte Zn concentration (fmol/", "μm"^2, ")"))) +
    scale_fill_manual(values = brewer.pal(4, "PuBu")[2:4]) +
    stat_summary(geom = "text", label = c("a", "b", "c"), size = 8, fun.y = max, vjust = -0.3)
HSD.test(aov(Zn_conc ~ stage, data = WT_Zn_OO), "stage", console = TRUE)

grid.draw(cbind(ggplotGrob(wteca), ggplotGrob(wtecc), ggplotGrob(wtooa), ggplotGrob(wtooc), size = "last"))

#Figure 4
Zn_EC <- norm_df %>%
    filter(focus == "cell",
           !(ID %in% c("2_0148", "2_0149", "2_0150")),
           stage != "in_vitro_egg") %>%
    mutate(Zn_mol = Zn_measure / 65.38, genotype = factor(genotype, levels = c("WT", "KO")),
           Zn_normmol = (Zn_measure / 65.38) / (Fe_measure / 55.845),
           Zn_normconc = Zn_conc / Fe_conc) %>%
    select(ID, stage, genotype, Zn_mol, Zn_conc, Zn_normmol, Zn_normconc)
wtkoeca <- ggplot(Zn_EC, aes(x = stage, y = Zn_mol, fill = genotype)) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 20, hjust = 1, angle = 90), axis.title.x = element_blank(),
          axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 30), panel.grid = element_blank()) +
    geom_boxplot(position = "dodge", outlier.colour = NA) +
    geom_point(position = position_dodge(width = 0.75)) +
    scale_x_discrete(limits = c("stage1-8_oocyte", "stage9-13_oocyte", "stage14_oocyte", "in_vivo_egg"),
                     labels = c("stage 1-8", "stage 9-13", "stage 14", "activated egg")) +
    ylab("Egg chamber Zn amount (fmol)") +
    scale_fill_manual(values = brewer.pal(4, "PuBu")[2:1]) +
    stat_summary(geom = "text", label = c("a", "a", "b", "bc", "d", "b", "cd", "cd"), size = 8, fun.y = max, vjust = -0.3,
                 position = position_dodge(width = 0.75))

wtkoecc <- ggplot(Zn_EC, aes(x = stage, y = Zn_conc, fill = genotype)) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 20, hjust = 1, angle = 90), axis.title.x = element_blank(),
          axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 25), panel.grid = element_blank()) +
    geom_boxplot(position = "dodge", outlier.colour = NA) +
    geom_point(position = position_dodge(width = 0.75)) +
    scale_x_discrete(limits = c("stage1-8_oocyte", "stage9-13_oocyte", "stage14_oocyte", "in_vivo_egg"),
                     labels = c("stage 1-8", "stage 9-13", "stage 14", "activated egg")) +
    ylab(expression(paste("Egg chamber Zn concentration (fmol/", "μm"^2, ")"))) +
    scale_fill_manual(values = brewer.pal(4, "PuBu")[2:1]) +
    stat_summary(geom = "text", label = c("a", "a", "b", "cd", "de", "bc", "e", "e"), size = 8, fun.y = max, vjust = -0.3,
                 position = position_dodge(width = 0.75))

txwtkozn <- with(Zn_EC, interaction(genotype, stage))
HSD.test(aov(Zn_normmol ~ txwtkozn, data = Zn_EC), "txwtkozn", console = TRUE)

wtkoecanorm <- ggplot(Zn_EC, aes(x = stage, y = Zn_normmol, fill = genotype)) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 20, hjust = 1, angle = 90), axis.title.x = element_blank(),
          axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 20), panel.grid = element_blank()) +
    geom_boxplot(position = "dodge", outlier.colour = NA) +
    geom_point(position = position_dodge(width = 0.75)) +
    scale_x_discrete(limits = c("stage1-8_oocyte", "stage9-13_oocyte", "stage14_oocyte", "in_vivo_egg"),
                     labels = c("stage 1-8", "stage 9-13", "stage 14", "activated egg")) +
    ylab("Normalized egg chamber Zn amount (Zn/Fe)") +
    scale_fill_manual(values = brewer.pal(4, "PuBu")[2:1]) +
    stat_summary(geom = "text", label = c("a", "a", "a", "b", "a", "b", "a", "a"), size = 8, fun.y = max, vjust = -0.3,
                 position = position_dodge(width = 0.75))


Zn_OO <- norm_df %>%
    filter(!(ID %in% c("2_0148", "2_0149", "2_0150"))) %>%
    filter((stage == "stage9-13_oocyte" & focus == "focus") |
               (stage == "stage14_oocyte" & focus == "cell") |
               (stage == "in_vivo_egg" & focus == "cell")) %>%
    mutate(Zn_mol = Zn_measure / 65.38, genotype = factor(genotype, levels = c("WT", "KO")),
           Zn_normmol = (Zn_measure / 65.38) / (Fe_measure / 55.845),
           Zn_normconc = Zn_conc / Fe_conc) %>%
    select(ID, stage, genotype, Zn_mol, Zn_conc, Zn_normmol, Zn_normconc)

wtkoooa <- ggplot(Zn_OO, aes(x = stage, y = Zn_mol, fill = genotype)) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 20, hjust = 1, angle = 90), axis.title.x = element_blank(),
          axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 30), panel.grid = element_blank()) +
    geom_boxplot(position = "dodge", outlier.colour = NA) +
    geom_point(position = position_dodge(width = 0.75)) +
    scale_x_discrete(limits = c("stage9-13_oocyte", "stage14_oocyte", "in_vivo_egg"),
                     labels = c("stage 9-13", "stage 14", "activated egg")) +
    ylab("Oocyte Zn amount (fmol)") +
    scale_fill_manual(values = brewer.pal(4, "PuBu")[2:1]) +
    stat_summary(geom = "text", label = c("ab", "a", "ab", "c", "b", "b"), size = 8, fun.y = max, vjust = -0.3,
                 position = position_dodge(width = 0.75))

wtkoooc <- ggplot(Zn_OO, aes(x = stage, y = Zn_conc, fill = genotype)) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 20, hjust = 1, angle = 90), axis.title.x = element_blank(),
          axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 30), panel.grid = element_blank(),
          legend.title = element_text(size = 30), legend.text = element_text(size = 30), legend.text.align = 0) +
    geom_boxplot(position = "dodge", outlier.colour = NA) +
    geom_point(position = position_dodge(width = 0.75)) +
    scale_x_discrete(limits = c("stage9-13_oocyte", "stage14_oocyte", "in_vivo_egg"),
                     labels = c("stage 9-13", "stage 14", "activated egg")) +
    ylab(expression(paste("Oocyte Zn concentration (fmol/", "μm"^2, ")"))) +
    labs(fill = "Genotype") +
    scale_fill_manual(values = brewer.pal(4, "PuBu")[2:1], labels = c(expression("WT"^" "), expression(italic("w;znt35C"^"1")))) +
    stat_summary(geom = "text", label = c("b", "a", "d", "c", "d", "d"), size = 8, fun.y = max, vjust = -0.3,
                 position = position_dodge(width = 0.75))

txwtkoozn <- with(Zn_OO, interaction(genotype, stage))
HSD.test(aov(Zn_normmol ~ txwtkoozn, data = Zn_OO), "txwtkoozn", console = TRUE)

wtkoooanorm <- ggplot(Zn_OO, aes(x = stage, y = Zn_normmol, fill = genotype)) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 20, hjust = 1, angle = 90), axis.title.x = element_blank(),
          axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 20), panel.grid = element_blank(),
          legend.title = element_text(size = 30), legend.text = element_text(size = 30), legend.text.align = 0) +
    geom_boxplot(position = "dodge", outlier.colour = NA) +
    geom_point(position = position_dodge(width = 0.75)) +
    scale_x_discrete(limits = c("stage9-13_oocyte", "stage14_oocyte", "in_vivo_egg"),
                     labels = c("stage 9-13", "stage 14", "activated egg")) +
    ylab("Normalized oocyte Zn amount (Zn/Fe)") +
    labs(fill = "Genotype") +
    scale_fill_manual(values = brewer.pal(4, "PuBu")[2:1], labels = c(expression("WT"^" "), expression(italic("w;znt35C"^"1")))) +
    stat_summary(geom = "text", label = c("b", "a", "c", "b", "c", "c"), size = 8, fun.y = max, vjust = -0.3,
                 position = position_dodge(width = 0.75))

grid.draw(cbind(ggplotGrob(wtkoeca), ggplotGrob(wtkoecc), ggplotGrob(wtkoooa), ggplotGrob(wtkoooc), size = "last"))
grid.draw(cbind(ggplotGrob(wtkoecanorm), ggplotGrob(wtkoooanorm), size = "last"))
