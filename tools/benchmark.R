# Run from interimlook/ after installing interimlook and ggplot2.
library(interimlook)
library(ggplot2)

comparison <- compare_interim(c(0.50, 0.75, 1))
labels <- c(obf = "O’Brien–Fleming type", pocock = "Pocock type", nominal = "Repeated nominal")
methods <- setNames(c("#4e50d7", "#008882", "#a34831"), labels)
long <- rbind(
  transform(comparison, metric = "Efficacy boundary (Z)", value = boundary_z),
  transform(comparison, metric = "Cumulative Type I error (%)", value = 100 * cumulative_alpha)
)
long$metric <- factor(long$metric, levels = c("Efficacy boundary (Z)", "Cumulative Type I error (%)"))
long$Method <- factor(labels[long$spending], levels = unname(labels))
target <- data.frame(metric = tail(levels(long$metric), 1), value = 2.5)

figure <- ggplot(long, aes(information, value, colour = Method)) +
  geom_hline(data = target, aes(yintercept = value), colour = "#8290a5", linetype = 3) +
  geom_line(aes(linetype = Method), linewidth = 1.05) + geom_point(size = 2.8) +
  facet_wrap(vars(metric), scales = "free_y") +
  scale_colour_manual(values = methods) +
  scale_linetype_manual(values = c("solid", "longdash", "solid")) +
  scale_x_continuous(breaks = c(0.5, 0.75, 1), labels = c("50%", "75%", "100%")) +
  labs(title = "Three looks at one efficacy hypothesis", subtitle = "One-sided Type I error target of 2.5%",
    x = "Information fraction", y = NULL, colour = NULL, linetype = NULL) +
  theme_minimal(base_size = 13) +
  theme(plot.background = element_rect(fill = "#f8f9fd", colour = NA),
    plot.margin = margin(17, 25, 10, 17), plot.title = element_text(face = "bold", colour = "#15243d"),
    plot.subtitle = element_text(colour = "#516178", margin = margin(b = 12)),
    panel.grid.minor = element_blank(), panel.spacing = grid::unit(1.7, "lines"),
    strip.text = element_text(face = "bold", colour = "#15243d"), legend.position = "bottom")

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
write.csv(comparison, "inst/extdata/interim-comparison.csv", row.names = FALSE)
ggsave("inst/extdata/interim-comparison.png", figure, width = 8.5, height = 4, dpi = 170)
