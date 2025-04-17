# --- Install required libraries if not already installed ---
required_packages <- c("ggplot2", "dplyr", "tidyr", "readr", "stringr")

installed_packages <- rownames(installed.packages())

for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# --- Load libraries ---
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)


# --- Load and combine all .windowed.pi files ---
pi_files <- list.files(pattern = "*.windowed.pi")

pi_data <- lapply(pi_files, function(file) {
  df <- read_tsv(file, col_types = cols())
  df$Population <- str_remove(file, ".windowed.pi$")
  return(df)
}) %>% bind_rows()

# --- Add π category for color coding ---
pi_data <- pi_data %>%
  mutate(pi_level = case_when(
    PI < 0.005 ~ "Low (<0.005)",
    PI > 0.015 ~ "High (>0.015)",
    TRUE ~ "Moderate"
  ))

# Now filter and plot
df <- pi_data %>% filter(Population == "Eastern_Australia")

ggplot(df, aes(x = BIN_START, y = PI)) +
  geom_line(alpha = 0.4, color = "gray50") +
  geom_point(aes(color = pi_level), size = 1.6) +
  geom_hline(yintercept = c(0.005, 0.015), linetype = "dashed", color = "gray40") +
  scale_color_manual(
    values = c("Low (<0.005)" = "blue", "Moderate" = "black", "High (>0.015)" = "red")
  ) +
  labs(title = "Nucleotide Diversity (π) — Eastern_Australia",
       x = "Genomic Window",
       y = "π (Nucleotide Diversity)",
       color = "Diversity Level") +
  theme_minimal(base_size = 14)

df <- pi_data %>% filter(Population == "Tropical_Atlantic")
ggplot(df, aes(x = BIN_START, y = PI)) +
  geom_line(alpha = 0.4, color = "gray50") +
  geom_point(aes(color = pi_level), size = 1.6) +
  geom_hline(yintercept = c(0.005, 0.015), linetype = "dashed", color = "gray40") +
  scale_color_manual(
    values = c("Low (<0.005)" = "blue", "Moderate" = "black", "High (>0.015)" = "red")
  ) +
  labs(title = "Nucleotide Diversity (π) — Tropical_Atlantic",
       x = "Genomic Window",
       y = "π (Nucleotide Diversity)",
       color = "Diversity Level") +
  theme_minimal(base_size = 14)

df <- pi_data %>% filter(Population == "Southwestern_Australia")
ggplot(df, aes(x = BIN_START, y = PI)) +
  geom_line(alpha = 0.4, color = "gray50") +
  geom_point(aes(color = pi_level), size = 1.6) +
  geom_hline(yintercept = c(0.005, 0.015), linetype = "dashed", color = "gray40") +
  scale_color_manual(
    values = c("Low (<0.005)" = "blue", "Moderate" = "black", "High (>0.015)" = "red")
  ) +
  labs(title = "Nucleotide Diversity (π) — Southwestern_Australia",
       x = "Genomic Window",
       y = "π (Nucleotide Diversity)",
       color = "Diversity Level") +
  theme_minimal(base_size = 14)

df <- pi_data %>% filter(Population == "Hawaii")
ggplot(df, aes(x = BIN_START, y = PI)) +
  geom_line(alpha = 0.4, color = "gray50") +
  geom_point(aes(color = pi_level), size = 1.6) +
  geom_hline(yintercept = c(0.005, 0.015), linetype = "dashed", color = "gray40") +
  scale_color_manual(
    values = c("Low (<0.005)" = "blue", "Moderate" = "black", "High (>0.015)" = "red")
  ) +
  labs(title = "Nucleotide Diversity (π) — Hawaii",
       x = "Genomic Window",
       y = "π (Nucleotide Diversity)",
       color = "Diversity Level") +
  theme_minimal(base_size = 14)

df <- pi_data %>% filter(Population == "Offshore")
ggplot(df, aes(x = BIN_START, y = PI)) +
  geom_line(alpha = 0.4, color = "gray50") +
  geom_point(aes(color = pi_level), size = 1.6) +
  geom_hline(yintercept = c(0.005, 0.015), linetype = "dashed", color = "gray40") +
  scale_color_manual(
    values = c("Low (<0.005)" = "blue", "Moderate" = "black", "High (>0.015)" = "red")
  ) +
  labs(title = "Nucleotide Diversity (π) — Offshore",
       x = "Genomic Window",
       y = "π (Nucleotide Diversity)",
       color = "Diversity Level") +
  theme_minimal(base_size = 14)

df <- pi_data %>% filter(Population == "Mexico_Pacific")
ggplot(df, aes(x = BIN_START, y = PI)) +
  geom_line(alpha = 0.4, color = "gray50") +
  geom_point(aes(color = pi_level), size = 1.6) +
  geom_hline(yintercept = c(0.005, 0.015), linetype = "dashed", color = "gray40") +
  scale_color_manual(
    values = c("Low (<0.005)" = "blue", "Moderate" = "black", "High (>0.015)" = "red")
  ) +
  labs(title = "Nucleotide Diversity (π) — Mexico_Pacific",
       x = "Genomic Window",
       y = "π (Nucleotide Diversity)",
       color = "Diversity Level") +
  theme_minimal(base_size = 14)

df <- pi_data %>% filter(Population == "Indian_Ocean")
ggplot(df, aes(x = BIN_START, y = PI)) +
  geom_line(alpha = 0.4, color = "gray50") +
  geom_point(aes(color = pi_level), size = 1.6) +
  geom_hline(yintercept = c(0.005, 0.015), linetype = "dashed", color = "gray40") +
  scale_color_manual(
    values = c("Low (<0.005)" = "blue", "Moderate" = "black", "High (>0.015)" = "red")
  ) +
  labs(title = "Nucleotide Diversity (π) — Indian_Ocean",
       x = "Genomic Window",
       y = "π (Nucleotide Diversity)",
       color = "Diversity Level") +
  theme_minimal(base_size = 14)


# Count high π per population and get percentages
high_pi_summary <- pi_data %>%
  group_by(Population) %>%
  summarise(
    Total = n(),
    High_Count = sum(pi_level == "High (>0.015)", na.rm = TRUE),
    Prop_High = High_Count / Total,
    SE_High = sqrt((Prop_High * (1 - Prop_High)) / Total)
  ) %>%
  mutate(
    Percent_High = Prop_High * 100,
    SE_Percent_High = SE_High * 100
  )

# Barplot with error bars (percent)
ggplot(high_pi_summary, aes(x = Population, y = Percent_High)) +
  geom_col(fill = "red") +
  geom_errorbar(aes(ymin = Percent_High - SE_Percent_High,
                    ymax = Percent_High + SE_Percent_High),
                width = 0.2, color = "black") +
  labs(title = "High π Windows (>0.015) per Ecotype",
       x = "Ecotype",
       y = "Percentage (%)") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

moderate_pi_summary <- pi_data %>%
  group_by(Population) %>%
  summarise(
    Total = n(),
    Moderate_Count = sum(pi_level == "Moderate", na.rm = TRUE),
    Prop_Moderate = Moderate_Count / Total,
    SE_Moderate = sqrt((Prop_Moderate * (1 - Prop_Moderate)) / Total)
  ) %>%
  mutate(
    Percent_Moderate = Prop_Moderate * 100,
    SE_Percent_Moderate = SE_Moderate * 100
  )

# Barplot with error bars (percent)
ggplot(moderate_pi_summary, aes(x = Population, y = Percent_Moderate)) +
  geom_col(fill = "black") +
  geom_errorbar(aes(ymin = Percent_Moderate - SE_Percent_Moderate,
                    ymax = Percent_Moderate + SE_Percent_Moderate),
                width = 0.2, color = "gray30") +
  labs(title = "Moderate π Windows (0.005–0.015) per Ecotype",
       x = "Ecotype",
       y = "Percentage (%)") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Count low π per population and get percentages
low_pi_summary <- pi_data %>%
  group_by(Population) %>%
  summarise(
    Total = n(),
    Low_Count = sum(pi_level == "Low (<0.005)", na.rm = TRUE),
    Prop_Low = Low_Count / Total,
    SE_Low = sqrt((Prop_Low * (1 - Prop_Low)) / Total)
  ) %>%
  mutate(
    Percent_Low = Prop_Low * 100,
    SE_Percent_Low = SE_Low * 100
  )

# Barplot with error bars (percent)
ggplot(low_pi_summary, aes(x = Population, y = Percent_Low)) +
  geom_col(fill = "blue") +
  geom_errorbar(aes(ymin = Percent_Low - SE_Percent_Low,
                    ymax = Percent_Low + SE_Percent_Low),
                width = 0.2, color = "black") +
  labs(title = "Low π Windows (<0.005) per Ecotype",
       x = "Ecotype",
       y = "Percentage (%)") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Count high π per population
high_pi_summary <- pi_data %>%
  group_by(Population) %>%
  summarise(
    Total = n(),
    High_Count = sum(pi_level == "High (>0.015)", na.rm = TRUE),
    SE_High = sqrt(High_Count)  # Poisson-style SE for counts
  )

# Barplot with error bars (counts)
ggplot(high_pi_summary, aes(x = Population, y = High_Count)) +
  geom_col(fill = "red") +
  geom_errorbar(aes(ymin = High_Count - SE_High,
                    ymax = High_Count + SE_High),
                width = 0.2, color = "black") +
  labs(title = "Number of High π Windows (>0.015) per Ecotype",
       x = "Ecotype",
       y = "Counts") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Count moderate π per population
moderate_pi_summary <- pi_data %>%
  group_by(Population) %>%
  summarise(
    Total = n(),
    Moderate_Count = sum(pi_level == "Moderate", na.rm = TRUE),
    SE_Moderate = sqrt(Moderate_Count)  # Poisson-style SE for count
  )

# Barplot with error bars (counts)
ggplot(moderate_pi_summary, aes(x = Population, y = Moderate_Count)) +
  geom_col(fill = "black") +
  geom_errorbar(aes(ymin = Moderate_Count - SE_Moderate,
                    ymax = Moderate_Count + SE_Moderate),
                width = 0.2, color = "gray30") +
  labs(title = "Number of Moderate π Windows (0.005–0.015) per Ecotype",
       x = "Ecotype",
       y = "Counts") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Count low π per population
low_pi_summary <- pi_data %>%
  group_by(Population) %>%
  summarise(
    Total = n(),
    Low_Count = sum(pi_level == "Low (<0.005)", na.rm = TRUE),
    SE_Low = sqrt(Low_Count)  # Poisson-style SE for count
  )

# Barplot with error bars (counts)
ggplot(low_pi_summary, aes(x = Population, y = Low_Count)) +
  geom_col(fill = "blue") +
  geom_errorbar(aes(ymin = Low_Count - SE_Low,
                    ymax = Low_Count + SE_Low),
                width = 0.2, color = "black") +
  labs(title = "Number of Low π Windows (<0.005) per Ecotype",
       x = "Ecotype",
       y = "Counts") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
