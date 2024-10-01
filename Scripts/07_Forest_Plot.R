library(grid)
library(magrittr)
library(checkmate)
library(forestplot)
library(forestploter)
library(dplyr)

# MR results no stratification
dt=read.csv(file="../results/data_full_model_summary.csv")
summary(dt)

dt$Variable <- ifelse(is.na(dt$estimate),
                       dt$Variable,
                       paste0("      ", dt$Variable))

dt$` ` <- paste(rep(" ", 20), collapse = " ")

dt$`OR (95% CI)` <- ifelse(is.na(dt$estimate), " ",
                           sprintf("%.2f (%.2f to %.2f)",
                                   dt$estimate, dt$X95_lo, dt$X95_up))
dt$`P value` <- ifelse(is.na(dt$estimate)|dt$estimate == 1, " ",
                           sprintf("%.4f",
                                   dt$p_value))

tm <- forest_theme(base_size = 10,
                   refline_col = "red",
                   footnote_col = "#636363",
                   footnote_fontface = "italic")

pdf(file = "../plots/data_full_model_summary.pdf") # The height of the plot in inches
forest(dt[, c(1, 6, 7, 8)], # this indicates which columns you want to plot
        est = dt$estimate,
        lower = dt$X95_lo,
        upper = dt$X95_up,
        xlim = c(0,5),
        ci_column = 2,
        ref_line = 1,
        arrow_lab = c("Negative Effect", "Positive Effect"),
        theme = tm
       )
dev.off()

# MR results stratification by gender - female
dt=read.csv(file="../results/data_female_model_summary.csv")
summary(dt)

dt$Variable <- ifelse(is.na(dt$estimate),
                      dt$Variable,
                      paste0("      ", dt$Variable))

dt$` ` <- paste(rep(" ", 20), collapse = " ")

dt$`OR (95% CI)` <- ifelse(is.na(dt$estimate), " ",
                           sprintf("%.2f (%.2f to %.2f)",
                                   dt$estimate, dt$X95_lo, dt$X95_up))
dt$`P value` <- ifelse(is.na(dt$estimate)|dt$estimate == 1, " ",
                       sprintf("%.4f",
                               dt$p_value))

tm <- forest_theme(base_size = 10,
                   refline_col = "red",
                   footnote_col = "#636363",
                   footnote_fontface = "italic")

pdf(file = "../plots/data_female_model_summary.pdf") # The height of the plot in inches
forest(dt[, c(1, 6, 7, 8)], # this indicates which columns you want to plot
       est = dt$estimate,
       lower = dt$X95_lo,
       upper = dt$X95_up,
       xlim = c(0,5),
       ci_column = 2,
       ref_line = 1,
       arrow_lab = c("Negative Effect", "Positive Effect"),
       theme = tm
)
dev.off()

# MR results stratification by gender - male
dt=read.csv(file="../results/data_male_model_summary.csv")
summary(dt)

dt$Variable <- ifelse(is.na(dt$estimate),
                      dt$Variable,
                      paste0("      ", dt$Variable))

dt$` ` <- paste(rep(" ", 20), collapse = " ")

dt$`OR (95% CI)` <- ifelse(is.na(dt$estimate), " ",
                           sprintf("%.2f (%.2f to %.2f)",
                                   dt$estimate, dt$X95_lo, dt$X95_up))
dt$`P value` <- ifelse(is.na(dt$estimate)|dt$estimate == 1, " ",
                       sprintf("%.4f",
                               dt$p_value))

tm <- forest_theme(base_size = 10,
                   refline_col = "red",
                   footnote_col = "#636363",
                   footnote_fontface = "italic")

pdf(file = "../plots/data_male_model_summary.pdf") # The height of the plot in inches
forest(dt[, c(1, 6, 7, 8)], # this indicates which columns you want to plot
       est = dt$estimate,
       lower = dt$X95_lo,
       upper = dt$X95_up,
       xlim = c(0,5),
       ci_column = 2,
       ref_line = 1,
       arrow_lab = c("Negative Effect", "Positive Effect"),
       theme = tm
)
dev.off()

# MR results stratified by gender
dt=read.csv(file="../results/data_gender_model_summary.csv")
summary(dt)

dt$Stratification <- ifelse(is.na(dt$estimate),
                   dt$Stratification,
                   paste0("      ", dt$Stratification))

dt$Variable <- ifelse(is.na(dt$estimate),
                      dt$Variable,
                      paste0("      ", dt$Variable))

dt$` ` <- paste(rep(" ", 20), collapse = " ")

dt$`OR (95% CI)` <- ifelse(is.na(dt$estimate), " ",
                           sprintf("%.2f (%.2f to %.2f)",
                                   dt$estimate, dt$X95_lo, dt$X95_up))
dt$`P value` <- ifelse(is.na(dt$estimate)|dt$estimate == 1, " ",
                       sprintf("%.4f",
                               dt$p_value))

tm <- forest_theme(base_size = 10,
                   refline_col = "red",
                   footnote_col = "#636363",
                   footnote_fontface = "italic")

pdf(file = "../plots/data_gender_model_summary.pdf", width = 10, height = 12) # The height of the plot in inches
forest(dt[, c(1, 2, 7, 8, 9)], # this indicates which columns you want to plot
       est = dt$estimate,
       lower = dt$X95_lo,
       upper = dt$X95_up,
       xlim = c(0,5),
       ci_column = 3,
       ref_line = 1,
       arrow_lab = c("Negative Effect", "Positive Effect"),
       theme = tm
)
dev.off()


# MR results conditioned
dt=read.csv(file="../results/conditioned_model_summary.csv")
summary(dt)

dt$Model <- ifelse(is.na(dt$estimate),
                      dt$Model,
                      paste0("      ", dt$Model))

dt$Variable <- ifelse(is.na(dt$estimate),
                      dt$Variable,
                      paste0("      ", dt$Variable))

dt$` ` <- paste(rep(" ", 20), collapse = " ")

dt$`OR (95% CI)` <- ifelse(is.na(dt$estimate), " ",
                           sprintf("%.2f (%.2f to %.2f)",
                                   dt$estimate, dt$X95_lo, dt$X95_up))
dt$`P value` <- ifelse(is.na(dt$estimate)|dt$estimate == 1, " ",
                       sprintf("%.4f",
                               dt$p_value))

tm <- forest_theme(base_size = 10,
                   refline_col = "red",
                   footnote_col = "#636363",
                   footnote_fontface = "italic")

pdf(file = "../plots/conditioned_model_summary.pdf", width = 10, height = 6) # The height of the plot in inches
forest(dt[, c(1, 2, 7, 8, 9)], # this indicates which columns you want to plot
       est = dt$estimate,
       lower = dt$X95_lo,
       upper = dt$X95_up,
       xlim = c(0,5),
       ci_column = 3,
       ref_line = 1,
       arrow_lab = c("Negative Effect", "Positive Effect"),
       theme = tm
)
dev.off()
