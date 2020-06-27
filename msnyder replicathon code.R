#load packages
install.packages("dplyr")
install.packages("tidyr")
install.packages("swirl")
install.packages("tidyverse")
library(tidyverse)

#plotting theme
theme_set(theme_bw())

#load data 
pharmacodata<-readRDS("./data/rawPharmacoData.rds")
head(pharmacodata)
str(pharmacodata)



#number of drugs use
length(unique(rawPharmacoData$drug))


#summarize number of drugs and number of cell lines
pharmacodata %>%
  summarize(nCellLines= n_distinct(cell_line), n_drugs=n_distinct(drug))

#count unique number of drug concentrations used
tapply(pharmacodata$concentration, pharmacodata$study,
       function(x) {length(unique(x))})

#another way to count unique drug concentrations in each study
pharmacodata %>%
  group_by(study) %>%
  summarize(n=n_distinct(concentration))

#histogram of distribution of concentrations by study
pharmacodata %>%
  ggplot(aes(x = log2(concentration))) +
  geom_histogram(fill = "gray", color = "black") +
  facet_wrap(~ study) +
  ggtitle("Distributions of concentrations by study")

pharmacodata %>%
  ggplot(aes(x = viability)) +
  geom_histogram(fill = "gray", color = "black", binwidth = 20) +
  facet_wrap(~ study) +
  ggtitle("Distributions of viability by study")

#determine the
range(pharmacodata$viability)

#number of trials with viability less than 0
sum(pharmacodata$viability<0)

#number of trials with viability greater than 100
sum(pharmacodata$viability>100)

#summarize the viability data by expected values vs. unexpected values
pharmacodata %>%
  summarize(min_viability = min(viability),
            max_viability = max(viability),
            n_too_small   = sum(viability < 0),
            n_too_big     = sum(viability > 100))

#distribution (like a normal curve) of viability scores
pharmacodata %>%
  ggplot(aes(x = viability, group = study, fill = study, color = study)) +
  geom_density(alpha = 1/4) +
  xlim(0, 170) +
  ggtitle("Distributions of viability scores by study")

#box plot distribution of viability scores with drug
gp <- pharmacodata %>%
  ggplot(aes(y = viability, x = drug, fill = study)) +
  scale_x_discrete() + 
  annotate(geom = "rect", ymin = 0, ymax = 100, xmin = -Inf, xmax = Inf,
           fill = 'black', alpha = 1/6) +
  geom_boxplot(outlier.alpha = 1/5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1/2)) +
  ggtitle("Distributions of viability scores by drug and study")
gp

#box plot of distribution of viability scores with cell line
gp <- pharmacodata %>%
  ggplot(aes(y = viability, x = concentration, fill = study)) +
  scale_x_discrete() + 
  annotate(geom = "rect", ymin = 0, ymax = 100, xmin = -Inf, xmax = Inf,
           fill = 'black', alpha = 1/6) +
  geom_boxplot(outlier.alpha = 1/5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1/2)) +
  ggtitle("Distributions of viability scores by dose and study")
gp

#drug concentration vs viability
ggplot(pharmacodata, aes(x = doseID, y = viability))+
  geom_histogram(alpha=1/2)

#load summarized data
summarizedData <- readRDS("./data/summarizedPharmacoData.rds")

#specific drug subset
azdSummary <- subset(summarizedData, drug == "AZD0530")

#ggplot of the two studies AUC data, should be a line if correlated
ggplot(azdSummary, aes(x = auc_GDSC, y = auc_CCLE)) +
  geom_point(alpha = 1/2) +
  xlab("GDSC AUC") +
  ylab("CCLE AUC") +
  ggtitle("AUC summaries of cell line response to AZD0530 across studies")

#IC50 graphed for correlation between the two studies
ggplot(azdSummary, aes(x = ic50_GDSC, y = ic50_CCLE)) +
  geom_point(alpha = 1/2) +
  xlab("GDSC IC50") +
  ylab("CCLE IC50") +
  ggtitle("IC50 summaries of cell line response to AZD0530 across studies")


#correctin gthe IC50 plot for a skewed distribution
ggplot(azdSummary, aes(x = -log10(ic50_GDSC), y = -log(ic50_CCLE))) +
  geom_point(alpha = 1/2) +
  xlab("-log10(GDSC IC50)") +
  ylab("-log10(CCLE IC50)") +
  ggtitle("IC50 summaries of cell line response to AZD0530 across studies")


#summarizing the IC50 data for all 15 drugs in each study
summarizedData %>%
  ggplot(aes(x = -log10(ic50_GDSC / 10^6),
             y = -log10(ic50_CCLE / 10^6))) +
  geom_point(cex = 1/2) + 
  facet_wrap(~ drug) +
  xlab("-log10(GDSC IC50/10^6)") +
  ylab("-log10(CCLE IC50/10^6)") +
  ggtitle("IC50 summaries of cell line response across studies")

#scattereplot for cell line response to drugs across studies

summarizedData %>%
  ggplot(aes(x = -log10(auc_GDSC / 10^6),
             y = -log10(auc_CCLE / 10^6))) +
  geom_point(cex = 1/2) + 
  facet_wrap(~ drug) +
  xlab("-log10(GDSC auc/10^6)") +
  ylab("-log10(CCLE auc/10^6)") +
  ggtitle("AUC summaries of cell line response to drugs across studies")

#summarized data correlation coefficients table by two methods
drugCorrs <- summarizedData %>% 
  group_by(drug) %>%
  summarize(Pearson_ic50  = cor(-log10(ic50_GDSC / 10^6), -log10(ic50_CCLE / 10^6), method = "pearson"),
            Spearman_ic50 = cor(-log10(ic50_GDSC / 10^6), -log10(ic50_CCLE / 10^6), method = "spearman"))
summarizedData <- readRDS("./data/summarizedPharmacoData.rds")
summarizedData %>%
  ggplot(aes(x = -log10(auc_GDSC / 10^6),
             y = -log10(auc_CCLE / 10^6))) +
  geom_point(cex = 1/2) + 
  facet_wrap(~ drug) +
  xlab("-log10(GDSC auc/10^6)") +
  ylab("-log10(CCLE auc/10^6)") +
  ggtitle("AUC summaries of cell line response to drugs across studies"


#prepping the correlation coefficient data for graph
drugCorrs

drugCorrs <- gather(drugCorrs, measure, correlation, -drug)

drugCorrs
#drug correlation cell line summaries for IC50 in both studies
drugCorrs %>%
  ggplot(aes(x = drug, y = correlation, fill = measure, group = measure)) +
  geom_bar(stat = "identity", position = position_dodge(), colour = "black") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_grey() +
  ylim(0, 1) + 
  ggtitle("Correlation of cell line IC50 summaries between studies for each drug")
