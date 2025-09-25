#Plot for distribution of protein lengths comparing positive and negative sequences
## Load tidyverse for easy data handling + plotting
library(tidyverse)

## Read the TSV file
df <- read_tsv("all_info.tsv")

df <- df %>%
  mutate(length = nchar(ProteinLength))

##creating the density plot
ggplot(df, aes(x = ProteinLength, fill = label)) +
  geom_density(alpha = 0.5) +
  labs(x = "Protein length (log scale)", y = "Density", fill = "Label") +
  theme_minimal() +
  facet_wrap(~Set) +
  scale_x_log10() +
  scale_fill_manual(values = c("positive" = "purple", "negative" = "green"))


#The distribution of SP lengths
library(tidyverse)

#read the TSV file
df <- read_tsv("all_info.tsv")

df<- df %>%
  filter(label == "positive", !is.na(SPPosition))

# calculate median SPPosition for each Set
medians <- df %>%
  group_by(Set) %>%
  summarize(median_SP = median(SPPosition))

##creating the density plot
ggplot(df, aes(x=SPPosition, fill=label)) +
  geom_density(alpha=0.5) +
  geom_vline(data = medians, aes(xintercept = median_SP), 
             color = "black", linetype = "dashed", size = 0.5) +
  labs(title="Signal Peptide Position Distribution", x="SP_
       Position", y="Density", fill="Label") +
  theme_minimal() +
  facet_wrap(~Set) +
  scale_x_log10() +
  scale_fill_manual(values = c("positive" = "red"))

