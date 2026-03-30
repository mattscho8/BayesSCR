#### Sequentially run the models
## 1. DA (x2: Stan and NIMBLE)
## 2. SCDL (x4: 2xStan, 2xNIMBLE)
## 3. ODL (x1: Stan)

#### Obtain timings and posterior summaries for each model
## Timing ignores compilation time

library(here)
library(ggplot2)
library(forcats)
library(dplyr)
library(posterior)
library(cowplot)

source(here("DA","DA_Stan.R"))
source(here("DA","DA_NIMBLE.R"))
source(here("SCDL","SCDL_Stan_BPP.R"))
source(here("SCDL","SCDL_Stan_PPP.R"))
source(here("SCDL","SCDL_NIMBLE_BPP.R"))
source(here("SCDL","SCDL_NIMBLE_PPP.R"))
source(here("ODL","ODL_Stan_PPP.R"))

### comparison of timings
mod_labels = c("DA", "SCDL_BPP", "SCDL_PPP", "ODL_PPP", "DA", "SCDL_BPP", "SCDL_PPP")
mod_timings = c(partime[3], partime2[3], partime3[3], partime4[3], outtime[3], outtime2[3], outtime3[3])
mod_software = c(rep("Stan",4), rep("NIMBLE",3))

plot_df = data.frame(model = factor(mod_labels, levels = c("DA", "SCDL_BPP", "SCDL_PPP", "ODL_PPP")), time = mod_timings, software = factor(mod_software,levels = c("Stan","NIMBLE")))

plt_t1 = ggplot(plot_df, aes(x = model, y = time, fill = software)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) + 
  theme_bw() + 
#  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Model", y = "Time (seconds)", fill = "Software") +
  scale_fill_manual(values = c("Stan" = "black", "NIMBLE" = "grey50"))
ggsave(here("./output/time.pdf"), plot = plt_t1, width = 7, height = 4, units = "in")

### Effective sample size / second. Facet by parameter
## Other = psi (DA), N (BPP) and phi/lambda (PPP)
sp = summary(parout)
sp2 = summary(parout2)
sp3 = summary(parout3)
sp4 = summary(parout4)
so = summary(out_df)
so2 = summary(out_df2)
so3 = summary(out_df3)

mod_ess = c(sp$ess_bulk/partime[3],sp2$ess_bulk/partime2[3], sp3$ess_bulk[c(1,2,4)]/partime3[3], sp4$ess_bulk[c(1,2,4)]/partime4[3], so$ess_bulk[2:4]/outtime[3], so2$ess_bulk/outtime2[3], so3$ess_bulk[2:4]/outtime3[3])
mod_par = c(sp$variable, sp2$variable, sp3$variable[c(1,2,4)], sp4$variable[c(1,2,4)], so$variable[2:4], so2$variable, so3$variable[2:4])
mod_mod = rep(mod_labels, each = 3)
plot_df2 = data.frame(model = factor(mod_mod, levels = c("DA", "SCDL_BPP", "SCDL_PPP", "ODL_PPP")), parameter = factor(mod_par,levels = c("p0","sigma","psi","N","phi")), ess = mod_ess, software = factor(c(rep("Stan", 12), rep("NIMBLE",9)), levels = c("Stan", "NIMBLE")
))
plot_df2$parameter = fct_recode(plot_df2$parameter, other = "psi", other = "N", other = "phi")


plt_ess = ggplot(plot_df2, aes(x = model, y = ess, fill = software)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
  facet_wrap(~ parameter, ncol = 1, labeller = label_parsed) +
  theme_bw() +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = "Model", y = "ESS/second", fill = "Software") +
  scale_fill_manual(values = c("Stan" = "black", "NIMBLE" = "grey50"))
ggsave(here("./output/ess.pdf"), plot = plt_ess, width = 7, height = 10, units = "in")


### Confirm that the parameter estimates are similar between software
summary_df = bind_rows(
  summarise_draws(parout) |> mutate(model = "DA", software = "Stan"),
  summarise_draws(parout2) |> mutate(model = "SCDL_BPP", software = "Stan"),
  summarise_draws(parout3)[c(1,2,4),] |> mutate(model = "SCDL_PPP", software = "Stan"),
  summarise_draws(parout4)[c(1,2,4),] |> mutate(model = "ODL_PPP", software = "Stan"),
  summarise_draws(out_df)[2:4,] |> mutate(model = "DA", software = "NIMBLE"),
  summarise_draws(out_df2) |> mutate(model = "SCDL_BPP", software = "NIMBLE"),
  summarise_draws(out_df3)[2:4,] |> mutate(model = "SCDL_PPP", software = "NIMBLE"),
)
summary_df$model = factor(summary_df$model, levels = c("DA", "SCDL_BPP", "SCDL_PPP", "ODL_PPP"))
summary_df$software = factor(summary_df$software, levels = c("Stan", "NIMBLE"))

summary_p0 = summary_df |> filter(variable == "p0")
summary_sigma = summary_df |> filter(variable == "sigma")
summary_psi = summary_df |> filter(variable == "psi")
summary_N = summary_df |> filter(variable == "N")
summary_phi = summary_df |> filter(variable == "phi")

# First p0
p1 = ggplot(summary_p0, aes(x = model, y = median, ymin = q5, ymax = q95, colour = software)) +
         geom_pointrange(position = position_dodge(preserve = "single", width = 0.5)) +
         theme_bw() +
         scale_colour_manual(values = c("Stan" = "black", "NIMBLE" = "grey50")) +
         labs(x = "", y = "") +
         facet_wrap(~ "p0")

# Next sigma
p2 = ggplot(summary_sigma, aes(x = model, y = median, ymin = q5, ymax = q95, colour = software)) +
  geom_pointrange(position = position_dodge(preserve = "single", width = 0.5)) +
  theme_bw() +
  scale_colour_manual(values = c("Stan" = "black", "NIMBLE" = "grey50")) +
  labs(x = "", y = "") + 
  facet_wrap(~ "sigma", labeller = label_parsed)

# Psi
p3 = ggplot(summary_psi, aes(x = model, y = median, ymin = q5, ymax = q95, colour = software)) +
  geom_pointrange(position = position_dodge(preserve = "single", width = 0.5)) +
  theme_bw() +
  scale_colour_manual(values = c("Stan" = "black", "NIMBLE" = "grey50")) +
  labs(x = "", y = "") + 
  facet_wrap(~ "psi", labeller = label_parsed)

# N
p4 = ggplot(summary_N, aes(x = model, y = median, ymin = q5, ymax = q95, colour = software)) +
  geom_pointrange(position = position_dodge(preserve = "single", width = 0.5)) +
  theme_bw() +
  scale_colour_manual(values = c("Stan" = "black", "NIMBLE" = "grey50")) +
  labs(x = "", y = "") + 
  facet_wrap(~ "N")

# phi
p5 = ggplot(summary_phi, aes(x = model, y = median, ymin = q5, ymax = q95, colour = software)) +
  geom_pointrange(position = position_dodge(preserve = "single", width = 0.5)) +
  theme_bw() +
  scale_colour_manual(values = c("Stan" = "black", "NIMBLE" = "grey50")) +
  labs(x = "", y = "") + 
  facet_wrap(~ "phi", labeller = label_parsed)

# Remove legends from all plots 
p1_nl <- p1 + theme(legend.position = "none") 
p2_nl <- p2 + theme(legend.position = "none") 
p3_nl <- p3 + theme(legend.position = "none") 
p4_nl <- p4 + theme(legend.position = "none")
p5_nl <- p5 + theme(legend.position = "none") 

# Extract legend from one plot 
legend <- get_legend(p1) 

# Build layout 
bottom_row <- plot_grid(p3_nl, p4_nl, p5_nl, nrow = 1, rel_widths = c(1,1,2)) 
main_plot <- plot_grid(p1_nl, p2_nl, bottom_row, nrow = 3)

# Add shared axis labels 
y_label <- ggdraw() + draw_label("Value", angle = 90, size = 10, x = 0.8) 
x_label <- ggdraw() + draw_label("Model", size = 10, y = 0.98) 

# Combine y label with main plot 
with_y <- plot_grid(y_label, main_plot, nrow = 1, rel_widths = c(0.05, 1)) 

# Add x label below 
with_xy <- plot_grid(with_y, x_label, ncol = 1, rel_heights = c(1, 0.05)) 

# Add legend 
final_plot <- plot_grid(with_xy, legend, nrow = 1, rel_widths = c(1, 0.15)) 
final_plot

ggsave(here("./output/est.pdf"), plot = final_plot, width = 7, height = 10, units = "in")
