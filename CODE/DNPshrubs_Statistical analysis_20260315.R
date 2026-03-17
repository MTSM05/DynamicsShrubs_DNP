#******************************************************************************#

# AIM: Statistical analysis

# T. Sánchez-Mejía
# Last update 2026/03/15

#******************************************************************************#

library(tidyverse)
library(bbmle)
library(glmmTMB)
library(ggeffects)
library(emmeans)
library(performance)
library(DHARMa)
library(vegan)
library(car)
library(cowplot)
library(ggrepel)
library(MuMIn)


setwd("")

# 1. Load and prepare data ####

load(file = "Data/coverniche_data.RData")

View(sp_data7)                                                                  # Species data for climatic (CD) data from 1950 to 2022
View(sp_data50)                                                                 # Species data for climatic (CD) and shrub cover data from 2007 to 2022 

## 1.1. Get species data  ##### 
sp_data7 = sp_data7 %>%
  mutate(
    year= as.numeric (sp_data7$year),
    zyear = scale(sp_data7$year,center = TRUE, scale = TRUE),                   # scaled data
    sd_year = attr(sp_data7$zyear, "scaled:scale")
  )

sp_data50 = sp_data50 %>% 
  mutate(
    zyear = scale(sp_data50$year,center = TRUE, scale = TRUE),
    sd_year = attr(sp_data50$zyear, "scaled:scale"))

## 1.2. Get community data and community CD (CWCD) ####
com_data = sp_data7 %>%
  mutate(CDcover = CD * cover) %>%
  group_by(year, site) %>%
  summarise(
    covert = sum(cover, na.rm = TRUE),
    CWCD = sum(CDcover, na.rm = TRUE) / covert,
    .groups = "drop")%>%
  mutate(zyear=scale(year,center = TRUE, scale = TRUE),
         sd_year = attr(zyear, "scaled:scale"))
         
## 1.3. Common code for plots ####

species_labels = c(
  "Cistus libanotis"           = expression(italic("C. libanotis")),
  "Erica scoparia"             = expression(italic("E. scoparia")),
  "Halimium calycinum"         = expression(italic("H. calycinum")),
  "Halimium halimifolium"      = expression(italic("H. halimifolium")),
  "Helichrysum serotinum"      = expression(italic("H. serotinum")),
  "Lavandula stoechas"         = expression(italic("L. stoechas")),
  "Salvia rosmarinus"          = expression(italic("S. rosmarinus")),
  "Stauracanthus genistoides"  = expression(italic("S. genistoides")),
  "Thymus mastichina"          = expression(italic("T. mastichina")),
  "Ulex australis"             = expression(italic("U. australis"))
)

species_colors = c(
  "Cistus libanotis"           = "#FF7F0EFF",
  "Erica scoparia"             = "#C7519CFF",
  "Halimium calycinum"         = "#D63A3AFF",
  "Halimium halimifolium"      = "#FFBF50FF",  
  "Helichrysum serotinum"      = "#BCBD22FF",
  "Lavandula stoechas"         = "#8A60B0FF",
  "Salvia rosmarinus"          = "#608434", #"#78A641FF"
  "Stauracanthus genistoides"  = "#00A896FF",#"#12A2A8FF",  
  "Thymus mastichina"          = "#B9D974FF",
  "Ulex australis"             = "#2A6F97FF"#"#1F83B4FF"
)
# 2. Statistical Analysis #######
## 2.1. Temporal trend in CD & CWCD#####
### 2.1.1. spp CD ~ year (1950-2022) ####

m_5n0 = glmmTMB(CD ~ zyear , data = sp_data50)
m_5n1 = glmmTMB(CD ~ zyear + (1|species) + (1|site) , data = sp_data50) 
m_5n2 = glmmTMB(CD ~ zyear + (1|species) , data = sp_data50)
m_5n3 = glmmTMB(CD ~ zyear +(1|site) , data = sp_data50)
m_5n4 = glmmTMB(CD ~ zyear + (0+zyear | species) , data = sp_data50)                 #Random slope of zyear by species (random slope but NO random intercept)
#m_5n5 = glmmTMB(CD ~ zyear + (zyear | species) , data = sp_data50)                 # SINGULARITY! Random slope of zyear by species (random slope and random intercept)
#m_5n6 = glmmTMB(CD ~ zyear + (zyear || species) , data = sp_data50)                 # SINGULARITY!Random slope of zyear by species (random slope and random intercept) but removes correlation between intercept and slope
AICtab(m_5n0,m_5n1,m_5n2,m_5n3,m_5n4, base=TRUE, weights=TRUE, logLik=TRUE)

# FINAL MODEL
m_5nf = m_5n1                                                                 
summary_model = summary(m_5nf)
summary_model
r.squaredGLMM(m_5nf)

# CHECK MODEL
res_m_5nf = simulateResiduals(m_5nf, n = 1000)                                  # Get residuals
plot(res_m_5nf)                                                                 # Plot for visual diagnostics
testUniformity(res_m_5nf)                                                       # Uniformity test (should be non-significant)
testDispersion(res_m_5nf)                                                       # Dispersion test (should be non-significant)
testOutliers(res_m_5nf)

site = unique(sp_data50$site)
species_list = unique(sp_data50$species)
results_m_5nf = expand_grid(                                                    #Test temporal autocorrelation
   site = unique(sp_data50$site),
   species = unique(sp_data50$species)) %>%
   pmap_dfr(function(site_i, species_i)  {
     
     data_ss = sp_data50 %>% filter(site == site_i, species == species_i)
     
     res_sub = recalculateResiduals(
       res_m_5nf,sel = (sp_data50$site == site_i & sp_data50$species == species_i))
     
     test = testTemporalAutocorrelation(
       res_sub,time = sort(unique(data_ss$zyear)),plot = FALSE)
     
     data.frame(site = site_i,species = species_i,p_value = test$p.value)
   })

# PLOT MODEL 
# Predictions from the model
sp_data50_avg = sp_data50 %>%                                                   # Summarize averages by site and year
  group_by(zyear,year) %>%
  summarize(avgCD = mean(CD), .groups = "drop")

zyears_obs = sort(unique(sp_data50_avg$zyear))
df = ggpredict(m_5nf, terms = c("zyear [zyears_obs]"))
center = mean(sp_data50$year)  # 2007–2022 mean
scale_ = sd(sp_data50$year)
df$year = df$x * scale_ + center                                                # Get year values non-scaled

# Plot
sp_data50$group = sp_data50$species
pm_5nf = ggplot(df, aes(year, predicted)) + 
  geom_line(color = "gray30", size = 0.8) +
  geom_point(data = sp_data50, aes(x = year, y = CD),
              inherit.aes = FALSE, color = "gray30", size = 0.7, alpha = 0.05) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),fill = "gray30", alpha = 0.3) +
  theme_classic() + theme(text = element_text(family = "sans"),
    plot.tag = element_text(size = 11, hjust = -0.1, vjust = 1.0),plot.tag.position = c(-0.03, 1.03),    
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.text.x = element_text(size = 9,angle = 45, hjust = 1),axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 12, margin = margin(r = 4)),
    axis.title.x = element_text(size = 12, margin = margin(t = 4)),
    axis.line = element_blank(), plot.title = element_blank(),
    plot.margin = unit(c(2, 2, 2, 2), "mm")) +  
  labs(x = "Year", y = "CD", tag = "a)") +
  scale_x_continuous(breaks = seq(1950, 2022, by = 10),limits = c(1951, 2022)) +
  scale_y_continuous(breaks = seq(0, 1.25, by = 0.25),limits = c(0, 1.32)) 
pm_5nf

### 2.1.2. community CWCD ~ year (2007-2022) #####
m_7w0 = glmmTMB(CWCD ~ as.numeric(zyear) , data = com_data)
m_7w1 = glmmTMB(CWCD ~ as.numeric(zyear) + (1|site) , data = com_data)
AICtab(m_7w0,m_7w1, base=TRUE, weights=TRUE, logLik=TRUE)

# FINAL MODEL
m_7wf = m_7w1
summary_model = summary(m_7wf)
summary_model
r.squaredGLMM(m_7wf)

# CHECK MODEL
res_m_7wf = simulateResiduals(m_7wf, n = 1000)
plot(res_m_7wf)                                                                 #Plot for visual diagnostics
testUniformity(res_m_7wf)                                                       #Uniformity test (should be non-significant)
testDispersion(res_m_7wf)                                                       #Dispersion test (should be non-significant)
testOutliers(res_m_7wf)

results_m_7wf = expand_grid(                                                 
   site = unique(com_data$site)) %>%
   pmap_dfr(function(site_i)  {
     
     data_ss = com_data %>% filter(site == site_i)
     
     res_sub = recalculateResiduals(
       res_m_7wf,sel = (com_data$site == site_i))
     
     test = testTemporalAutocorrelation(
       res_sub,time = sort(unique(data_ss$zyear)),plot = FALSE)
     
     data.frame(site = site_i,p_value = test$p.value)
   })

# PLOT MODEL 
n1=ggpredict(m_7wf, terms = "zyear")
center = mean(unique(com_data$year))  # 2007–2022 mean
scale_ = sd(unique(com_data$year))
n1$year = n1$x * scale_ + center

p_w72 = ggplot(n1, aes(year, predicted)) + 
  geom_point(data = com_data, aes(x = year, y = CWCD),
             inherit.aes = FALSE, color = "gray30", size = 1.5, alpha = 0.7) +
  geom_line(color="gray30", size=0.8) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "gray30", alpha = 0.3) +  # Add ribbon with confidence intervals
   theme_classic() + theme(text = element_text(family = "sans"),
         plot.tag = element_text(size = 11, hjust = -0.1, vjust = 1.0),plot.tag.position = c(-0.03, 1.03),    
         panel.border = element_rect(color = "black", fill = NA, size = 0.5),
         axis.text.x = element_text(size = 9,angle = 45, hjust = 1),axis.text.y = element_text(size = 9),
         axis.title.y = element_text(size = 12, margin = margin(r = 4)),
         axis.title.x = element_text(size = 12, margin = margin(t = 4)),
         axis.line = element_blank(),plot.title = element_blank(),
         plot.margin = unit(c(2, 2, 2, 2), "mm")) +                             # top, right, bottom, left
  labs(x = "Year", y = expression(CWCD), tag = "b)")+
  scale_x_continuous(breaks = seq(2007, 2022, by = 2)) +
  scale_y_continuous(breaks = seq(0, 1.25, by = 0.25), limits = c(0, 1.25)
                     ) 
p_w72

### 2.1.3. spp CD ~ year (2007-2022) ####
m_7s = glmmTMB(CD ~ as.numeric(zyear) + species +as.numeric(zyear)*species +(1|site), family= gaussian(), data = sp_data7)              
m_7s1 = glmmTMB(CD ~ as.numeric(zyear) + species +as.numeric(zyear)*species, family= gaussian(), data = sp_data7)           
AICtab(m_7s,m_7s1,base=TRUE, weights=TRUE, logLik=TRUE)

# FINAL MODEL
m_7sf = m_7s                                                                   
summary_model = summary(m_7sf)
summary_model
Anova(m_7sf, type = 3)
r.squaredGLMM(m_7sf)

# CHECK MODEL
res_m_7sf = simulateResiduals(m_7sf, n = 1000)
plot(res_m_7sf)                                                                       # Plot for visual diagnostics
testUniformity(res_m_7sf)                                                             # Uniformity test (should be non-significant)
testDispersion(res_m_7sf)                                                             # Dispersion test (should be non-significant)
testOutliers(res_m_7sf)

results_m_7sf = expand_grid(                                                    # Test temporal autocorrelation
  site = unique(sp_data7$site),
  species = unique(sp_data7$species)) %>%
  pmap_dfr(function(site_i, species_i)  {
    
    data_ss = sp_data7 %>% filter(site == site_i, species == species_i)
    
    res_sub = recalculateResiduals(
      res_m_7sf,sel = (sp_data7$site == site_i & sp_data7$species == species_i))
    
    test = testTemporalAutocorrelation(
      res_sub,time = sort(unique(data_ss$zyear)),plot = FALSE)
    
    data.frame(site = site_i,species = species_i,p_value = test$p.value)
  })

# PLOT MODEL 
pred_species = ggpredict(m_7sf, terms = c("zyear", "species"))                  # Predictions
center = mean(unique(sp_data7$year))  # 2007–2022 mean
scale_ = sd(unique(sp_data7$year))
pred_species$year = pred_species$x * scale_ + center
sp_data7$group = sp_data7$species

p_sp72 = ggplot(pred_species, aes(x = year, y = predicted, color = group)) +
  geom_point(data = sp_data7, aes(x = year, y = CD, color = group),
             inherit.aes = FALSE, size = 1.5, alpha = 0.7) +
  geom_line(size = 0.8) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3, color = NA) +
    theme_classic() +theme(
    plot.tag = element_text(size = 11, hjust = -0.1, vjust = 1.0),
    plot.tag.position = c(-0.03, 1.03),
    legend.title = element_text(size = 11),
    legend.key.size = unit(0.9, "lines"),                                     # smaller key boxes
    legend.spacing.x = unit(0.55, 'cm'), legend.spacing.y = unit(0.45, 'cm'), # spacing
    legend.margin = margin(4, 6, 4, 4),                                       # legend margins
    legend.position = "right",
    legend.background = element_rect(color = "grey40", size = 0.3, fill = NA),
    legend.text = element_text(size = 11, margin = margin(l = 1.7)),
    text = element_text(family = "sans"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.text.x = element_text(size = 9,angle = 45, hjust = 1),axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 12, margin = margin(r = 4)),
    axis.title.x = element_text(size = 12, margin = margin(t = 4)),
    axis.line = element_blank(), plot.title = element_blank(),
    plot.margin = unit(c(2, 2, 4, 2), "mm")) +                                  # top, right, bottom, left
  labs(x = "Year", y = "CD", color = "Species", fill = "Species", tag = "c)") +
  scale_x_continuous(breaks = seq(2007, 2022, by = 2)) +
  scale_y_continuous(breaks = seq(0, 1.25, by = 0.25), limits = c(0, 1.25))+
  scale_color_manual(values = species_colors, labels = species_labels) +
  scale_fill_manual(values = species_colors, labels = species_labels) 
 p_sp72

#### 2.1.3.1. Slopes (emtrends) + Pearson correlation #####

# Get slopes (year * speciesyear.trend)
em_slopes = emtrends(m_7sf, ~ species, var = "zyear") 
em = summary(em_slopes, infer = TRUE) %>%                                       # Slopes and p values
  as.data.frame()

# Back-transform the trends 
em = em %>%
  mutate(
      CDT = zyear.trend / unique(sp_data7$sd_year),
      SE.z = SE,
      SE = SE.z / unique(sp_data7$sd_year),
      lower.CL.z = lower.CL,
      lower.CL = lower.CL.z / unique(sp_data7$sd_year),
      upper.CL.z = upper.CL,
      upper.CL = upper.CL.z / unique(sp_data7$sd_year)
      )

# Join slopes to sp_data7 for further analysis (Section 2.4)
em_clean = em %>%
  select("CDT","species")
sp_data7 = merge(sp_data7, em_clean, by="species")

# Get average values of observed CD in initial sampling year (2007)
avg_drodata = sp_data7%>%  
  group_by(species)%>%
  filter(year=="2007") %>%
  summarise(
    mnCD = mean(CD),                                                            # Get average across sites 
    seCD = sd(CD) / sqrt(n())                                                   # Standard error formula
  )

em2= merge(em,avg_drodata, by="species")                                        # Join slopes and initial values

# Do correlation
corr = cor.test(em2$CDT, em2$mnCD, method = "pearson")
corr

# PLOT CORRELATION
corplot = ggplot(em2, aes(x = as.numeric(mnCD), y = as.numeric(CDT), color = species)) +
  geom_errorbar(aes(ymin = as.numeric(CDT) - SE, 
                    ymax = as.numeric(CDT) + SE),
                 width = 0.0005) +
  geom_errorbarh(aes(xmin = as.numeric(mnCD) - seCD,
                      xmax = as.numeric(mnCD) + seCD), 
                 height = 0.0005)+
  geom_point(size = 2.5, alpha = 1)+
  geom_smooth(method = "lm", se = TRUE, color="gray30", size=0.8) +
  scale_color_manual(values = species_colors,labels = species_labels,name = "Species") +
  scale_fill_manual(values = species_colors,labels = species_labels,name = "Species") +
  labs(x = "Mean CD 2007",y = "CDT", tag = "d)") +
  theme_classic() +theme(
    plot.tag = element_text(size = 11, hjust = -0.1, vjust = 1.0),plot.tag.position = c(-0.03, 1.03),    
    legend.title = element_text(size = 11),
    legend.key.size = unit(0.9, "lines"),                                     # smaller key boxes
    legend.spacing.x = unit(0.55, 'cm'), legend.spacing.y = unit(0.45, 'cm'), # spacing
    legend.margin = margin(4, 6, 4, 4),                                       # legend margins
    legend.position = "bottom",
    legend.background = element_rect(color = "grey40", size = 0.3, fill = NA),
    legend.text = element_text(size = 10, margin = margin(l = 1.7)),
    axis.text.x = element_text(size = 9, hjust = 1,angle = 45),axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 12, margin = margin(r = 4)),
    axis.title.x = element_text(size = 12, margin = margin(t = 4)),
    text = element_text(family = "sans"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.line = element_blank(), plot.title = element_blank(),
    legend.text.align = 0,
    plot.margin = unit(c(2, 2, 4, 2), "mm")) +                             # top, right, bottom, left
  scale_y_continuous(breaks = seq(-0.000, 0.020, by = 0.01),labels = scales::number_format(accuracy = 0.01),
    limits = c(-0.001, 0.021))+
  scale_x_continuous(labels  = scales::number_format(accuracy = 0.01)
                     )
corplot

### 2.1.4. GRAPH -> Figure 2 #####
# Manage legends
pm_5nf_noleg = pm_5nf + theme(legend.position = "none")
p_w72_noleg  = p_w72 + theme(legend.position = "none")
p_sp72_noleg  = p_sp72 + theme(legend.position = "none")
corplot_noleg  = corplot + theme(legend.position = "none")
legend_cov = cowplot::get_legend(corplot + theme(legend.position = "bottom"))

empty = cowplot::ggdraw() 
combined_plot1 = cowplot::plot_grid(
  p_sp72_noleg,
  empty,
  corplot_noleg,
  empty,
  ncol = 4,
  rel_widths = c(1,0.03,1,0.02)
)
combined_plot2 = cowplot::plot_grid(
  pm_5nf_noleg,
  empty,
  p_w72_noleg,
  empty,
  ncol = 4,
  rel_widths = c(1,0.03,1,0.02)
)

combined_plot = (combined_plot2 / combined_plot1) +
   plot_layout(
     heights = c(1, 1),             
     guides = "keep")
 combined_plot

 combined_plot3 = cowplot::plot_grid(
   combined_plot,
   legend_cov,
   ncol = 1,
   rel_heights = c(2, 0.15)
 )
 
combined_plot_with_margin = combined_plot3 +
  theme(plot.margin = margin(t = 2, r = 2, b = 5, l = 2))  # Large bottom margin
 
#save data
pdf(file = "Figures/Figure2.CDevol.pdf",
    width = 7, height = 6)

combined_plot_with_margin
dev.off()

png(file = "Figures/Figure2.CDevol.png",
    res = 1000,                                                          
    width = 7, height = 6,                                                   
    units = "in",                                                           
    bg = "white")     

combined_plot_with_margin
dev.off()

## 2.2. Temporal trend shrub cover -> Univariate analysis ####
### 2.2.1. community cover ~ year (2007-2022)####
hist(com_data$covert)
m_covyearsite= glmmTMB(covert ~ as.numeric(zyear) +(1|site), family= nbinom2(), data = com_data)
m_covyearsitep= glmmTMB(covert ~ as.numeric(zyear) +(1|site), family= poisson(), data = com_data)              
m_covyearsite1= glmmTMB(covert ~ as.numeric(zyear), family= nbinom2(), data = com_data)   
AICtab(m_covyearsite,m_covyearsite1,m_covyearsitep,base=TRUE, weights=TRUE, logLik=TRUE)

# FINAL MODEL
m_covyearsitef = m_covyearsite
summary(m_covyearsitef)
r.squaredGLMM(m_covyearsitef)                                                 

# CHECK MODEL
res_m_covyearsitef = simulateResiduals(m_covyearsitef, n = 1000)
plot(res_m_covyearsitef)                                                        # Plot for visual diagnostics
testUniformity(res_m_covyearsitef)                                              # Uniformity test (should be non-significant)
testDispersion(res_m_covyearsitef)                                              # Dispersion test (should be non-significant)
testZeroInflation(res_m_covyearsitef)                                           # Zero-inflation test (useful if counts have many zeros), (should be non-significant)
testOutliers(res_m_covyearsitef)

results_m_covyearsitef  = expand_grid(
  site = unique(com_data$site)) %>%
  pmap_dfr(function(site_i, species_i)  {
    
    data_ss = com_data %>% filter(site == site_i)
    
    res_sub = recalculateResiduals(
      res_m_covyearsitef,sel = (com_data$site == site_i))
    
    test = testTemporalAutocorrelation(
      res_sub,time = sort(unique(data_ss$zyear)),plot = FALSE)
    
    data.frame(site = site_i,p_value = test$p.value)
  })

# PLOT MODEL 
pred_site = ggpredict(m_covyearsite, terms = c("zyear"))                        # Predictions
center = mean(unique(com_data$year))                                            # 2007–2022 mean
scale_ = sd(unique(com_data$year))
pred_site$year = pred_site$x * scale_ + center

p_totcovyear = ggplot(pred_site, aes(year, predicted)) + 
  geom_point(data = com_data, aes(x = year, y = covert), color = "gray30", size = 1.5, alpha = 0.7) +
  geom_line(color = "gray30", size = 0.8) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),fill = "gray30", alpha = 0.3) +
  theme_classic() + theme(text = element_text(family = "sans"),
                          plot.tag = element_text(size = 11, hjust = -0.1, vjust = 1.0),
                          plot.tag.position = c(-0.03, 1.03),    
                          panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                          axis.text.x = element_text(size = 9,angle = 45, hjust = 1),axis.text.y = element_text(size = 9),
                          axis.title.y = element_text(size = 12, margin = margin(r = 2)),
                          axis.title.x = element_text(size = 12, margin = margin(t = 4)),
                          axis.line = element_blank(), plot.title = element_blank(),
                          plot.margin = unit(c(2, 2, 4, 2), "mm")) +  
  labs(x = "Year", y = "Community cover", tag = "a)") +
  scale_x_continuous(breaks = seq(2007, 2022, by = 2)) +
  scale_y_continuous(breaks = seq(0, 120, by = 20),limits = c(0, 120))
p_totcovyear

### 2.2.2. spp cover ~ year (2007-2022)####

sp_data7_filtered = sp_data7 %>%                                                # Remove species-sites combinations without data
  group_by(site, species) %>%
  filter(sum(cover) > 0) %>%
  ungroup()
hist(sp_data7_filtered$cover)

# Try Poisson and Negative Binomial
m_covyearp= glmmTMB(cover ~ as.numeric(zyear) + species +as.numeric(zyear)*species +(1|site),REML=TRUE, family= poisson(), data = sp_data7_filtered)              
m_covyear= glmmTMB(cover ~ as.numeric(zyear) + species +as.numeric(zyear)*species +(1|site),REML=TRUE, family= nbinom2(), data = sp_data7_filtered)              
AICtab(m_covyearp,m_covyear,base=TRUE, weights=TRUE, logLik=TRUE)

# Try random effects
m_covyear= glmmTMB(cover ~ as.numeric(zyear) + species +as.numeric(zyear)*species +(1|site), family= nbinom2(),REML=TRUE, data = sp_data7_filtered)              
m_covyear1= glmmTMB(cover ~ as.numeric(zyear) + species +as.numeric(zyear)*species, family= nbinom2(),REML=TRUE, data = sp_data7_filtered)              
AICtab(m_covyear,m_covyear1)

# FINAL MODEL
m_covyearf = m_covyear
summary(m_covyearf)
Anova(m_covyearf, type = 3)
r.squaredGLMM(m_covyearf)


# CHECK MODEL
res_m_covyearf = simulateResiduals(m_covyearf, n = 1000)
plot(res_m_covyearf)                                                                       # Plot for visual diagnostics
testUniformity(res_m_covyearf)                                                             #Uniformity test (should be non-significant)
testDispersion(res_m_covyearf)                                                             #Dispersion test (should be non-significant)
testZeroInflation(res_m_covyearf)                                                          #Zero-inflation test (useful if counts have many zeros), (should be non-significant)
testOutliers(res_m_covyearf)

results_m_covyearf = expand_grid(
  site_i = unique(sp_data7_filtered$site),
  species_i = unique(sp_data7_filtered$species)) %>%
  pmap_dfr(function(site_i, species_i)  {
    
    data_ss = sp_data7_filtered %>% filter(site == site_i, species == species_i)
    
    res_sub = recalculateResiduals(
      res_m_covyearf,sel = (sp_data7_filtered$site == site_i & sp_data7_filtered$species == species_i))
    
    test = testTemporalAutocorrelation(
      res_sub,time = sort(unique(data_ss$zyear)),plot = FALSE)
    
    data.frame(site = site_i,species = species_i,p_value = test$p.value)
  })

# PLOT MODEL 
pred_species = ggpredict(m_covyearf, terms = c("zyear", "species"))                  # Predictions
center = mean(unique(sp_data7$year))  # 2007–2022 mean
scale_ = sd(unique(sp_data7$year))  
pred_species$year = pred_species$x * scale_ + center

pred_species$group = factor(pred_species$group, levels = names(species_colors))
sp_data7$group = factor(sp_data7$species, levels = names(species_colors))

p_covyear = ggplot(pred_species, aes(x = year, y = predicted, color = group)) +
  geom_point(data = sp_data7, aes(x = year, y = cover, color = group),
             inherit.aes = FALSE, size = 1.5, alpha = 0.7) +
  geom_line(size = 0.8) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3, color = NA) +
    theme_classic() + theme(
    plot.tag = element_text(size = 11, hjust = -0.1, vjust = 1.0),
    plot.tag.position = c(-0.03, 1.03),    
    legend.title = element_text(size=10),
    element_text(size = 10, margin = margin(l = 1.7)),
    legend.key.size = unit(0.9, "lines"),                                     # smaller key boxes
    legend.spacing.x = unit(0.55, 'cm'), legend.spacing.y = unit(0.45, 'cm'), # spacing
    legend.margin = margin(4, 6, 4, 4),                                       # legend margins
    legend.position = "right",
    legend.background = element_rect(color = "grey40", size = 0.3, fill = NA),
    text = element_text(family = "sans"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.text.x = element_text(size = 9,angle = 45, hjust = 1),axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 12, margin = margin(r = 4)),
    axis.title.x = element_text(size = 12, margin = margin(t = 4)),
    axis.line = element_blank(), plot.title = element_blank(),
    plot.margin = unit(c(2, 2, 4, 2), "mm")) +                                  # top, right, bottom, left
  labs(x = "Year", y = "Species cover", color = "Species", fill = "Species", tag = "b)") +
  scale_x_continuous(breaks = seq(2007, 2022, by = 2)) +
  scale_y_continuous(breaks = seq(0, 60, by = 10), limits = c(0,63))+
  scale_color_manual(values = species_colors, labels = species_labels) +
  scale_fill_manual(values = species_colors, labels = species_labels
                    ) 
p_covyear

#### 2.2.2.1. Slopes (emtrends) ########
# Do the emtrends
sloz = emtrends(m_covyearf, ~ species, var = "zyear", type = "response", test = TRUE)
test_results = test(sloz)
test_results

# Back-transform the trends 
slo_year = summary(sloz)
slo_year$emtrend = slo_year$zyear.trend / unique(sp_data7_filtered$sd_year)
slo_year$SE      = slo_year$SE / unique(sp_data7_filtered$sd_year)
slo_year$lower.CL.zy = slo_year$asymp.LCL
slo_year$asymp.LCL = slo_year$lower.CL.zy / unique(sp_data7_filtered$sd_year)
slo_year$upper.CL.zy = slo_year$asymp.UCL
slo_year$asymp.UCL = slo_year$upper.CL.zy / unique(sp_data7_filtered$sd_year)

### 2.2.3. GRAPH -> Figure 3 #####

# Manage legends
p_totcovyear_noleg = p_totcovyear + theme(legend.position = "none")
p_covyear_noleg = p_covyear + theme(legend.position = "none")
legend_cov = cowplot::get_legend(p_covyear + theme(legend.position = "bottom"))

empty = cowplot::ggdraw() 

top_row = cowplot::plot_grid(
  p_totcovyear_noleg,
  empty,
  p_covyear_noleg,
  empty,
  ncol = 4,
  rel_widths = c(1,0.03,1,0.02),
  align = "hv"
)

combined_plot = cowplot::plot_grid(
  top_row,
  legend_cov,
  ncol = 1,
  rel_heights = c(1, 0.15)
)

combined_plot
combined_plot_with_margin = combined_plot +
  theme(plot.margin = margin(t = 2, r = 2, b = 5, l = 2))                       # Large bottom margin


#save data
pdf(file = "Figures/Figure3.Shrubevol.pdf",
    width = 7, height = 3.5)

combined_plot_with_margin
dev.off()

png(file = "Figures/Figure3.Shrubevol.png",
    res = 1000,                                                          
    width = 7, height = 3.5,                                                   
    units = "in",                                                           
    bg = "white")     

combined_plot_with_margin
dev.off()

## 2.3. Temporal trend shrub cover -> Multivariate analysis ####
#### 2.3.1. Dissimilarity matrix #####

spdf = as.data.frame(sp_data7) %>%                              
  pivot_wider( 
    id_cols = c(site, census, year),                                            # Columns to keep
    names_from = species, 
    values_from = cover,
    values_fill = 0) %>%                                                        # fill missing species with 0
  mutate(
    site = factor(site),
    census = factor(census))

spmatrix = spdf %>%
  dplyr::select(-site, -census, -year) %>%                                      # Remove metadata columns
  as.matrix()

# Squareroot transform data (to diminish the influence of dominant species)
spmatrix_sq = sqrt(spmatrix)

# Get dissimilarity matrix
dmatrix = vegdist(spmatrix_sq, method = "bray")

# Create metadata dataframe (for PERMANOVA)
metadata = spdf %>%
  dplyr::select(site, census, year) %>%
  mutate(
    year = factor(year),                                                        
    year_numeric = as.numeric(as.character(year)),                              
    time_numeric = year_numeric - min(year_numeric)) %>%
  dplyr::select(-year_numeric)

#### 2.3.2. PERMANOVA #####

permanova_year = adonis2(dmatrix ~ time_numeric, 
                          data = metadata,
                          permutations = 9999,
                          strata = metadata$site
                          # ,method = "bray"                                    #unneeded with distance matrix (bray)
                         )
print(permanova_year)

# Dispersion check 

bd = betadisper(dmatrix, metadata$year)                                         #If non-significant, PERMANOVA result is valid
anova(bd)

#### 2.3.3. NMDS ####
nmds = metaMDS(spmatrix_sq, distance = "bray", k = 2, trymax = 100, autotransform = FALSE)             
nmds$stress
# Extract coordinates and combine with metadata
nmds_points = as.data.frame(nmds$points)
nmds_points = cbind(nmds_points, spdf)

# SCORES
site_scores = as.data.frame(scores(nmds, display = "sites"))
site_scores$site = metadata$site
site_scores$year = as.numeric(as.character(metadata$year))

# Rename for plot
colnames(site_scores)[1:2] = c("MDS1", "MDS2")

# Calculate species vector
sp_fit = envfit(nmds, spmatrix_sq, permutations = 999)

sp_scores = as.data.frame(scores(sp_fit, display = "vectors"))
sp_scores$species = rownames(sp_scores)
sp_scores$pval = sp_fit$vectors$pvals

# Keep significant species only 
sp_sig = subset(sp_scores, pval < 0.05)
sp_sig$short = c("C. libanotis","E. scoparia",
                 #"H. calycinum",
                          "H. halimifolium","H. serotinum","L. stoechas",
                          "S. rosmarinus","S. genistoides","T. mastichina",
                          "U. australis")

# NMDS_plot
site_scores$year = factor(site_scores$year,
                           levels = sort(unique(site_scores$year)))

NMDS_plot = ggplot(site_scores,
  aes(x = MDS1, y = MDS2, group = site,color = site, size = year)) +
  # Site trajectories
  geom_path(alpha = 0.4, linewidth = 0.35) +
  geom_point(alpha = 0.8) +
  # Color and size scales
  scale_color_viridis_d(option = "turbo", name = "Plot") +
  scale_size_manual(name = "Year",values = c(0.9,1.3,1.6,2,2.5,3,3.5,4,4.7)) +
  # Envfit species vectors
  geom_segment(data = sp_sig,
    aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.25, "cm")),
    linewidth = 0.9,
    color = "black") +
  # Species labels
  geom_text_repel(
    data = sp_sig,
    aes(x = NMDS1, y = NMDS2, label = short),
    inherit.aes = FALSE,
    fontface = "italic",size = 4,color = "black") +
  guides(size = guide_legend(override.aes = list(alpha = 1))) +
  theme_classic() +
  theme(
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 10),
    legend.key.size = unit(0.9, "lines"),
    legend.spacing.x = unit(0.55, "cm"),legend.spacing.y = unit(0.45, "cm"),
    legend.margin = margin(4, 6, 4, 4),
    legend.position = "right",
    legend.background = element_rect(color = "grey40", size = 0.3, fill = NA),
    text = element_text(family = "sans", size = 10),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.text.x = element_text(size = 11, angle = 0, hjust = 1),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 12, margin = margin(r = 4)),
    axis.title.x = element_text(size = 12, margin = margin(t = 4)),
    axis.line = element_blank(),
    plot.title = element_blank(),plot.margin = unit(c(2, 2, 4, 2), "mm")) +
   labs(x = "NMDS1",y = "NMDS2")
NMDS_plot

#save data
pdf(file = "Figures/FigureS6.NMDS.pdf",
    width = 8, height = 6.5)

NMDS_plot
dev.off()

png(file = "Figures/FigureS6.NMDS.png",
    res = 1000,                                                          
    width = 8, height = 6.5,                                                   
    units = "in",                                                           
    bg = "white")     

NMDS_plot
dev.off()


#### 2.3.4. SIMPER #######

years_to_compare = metadata$year %in% c(2007, 2022)  
simp = simper(spmatrix_sq[years_to_compare, ],                                  # Square-root data
               group = metadata$year[years_to_compare],
               permutations = 999)

# Get summary table
simp_df = as.data.frame(summary(simp)[[1]])
simp_df$species = rownames(simp_df)
simp$`2007_2022`$overall # Get overall dissimilarity

## 2.4. Cover ~ CDT ####
sp_data = split(sp_data7, sp_data7$species)

### 2.4.1. Get data ####
# Calculate net gain/loss of cover from 2007 to 2022

# Species
cover_diff = sp_data7 %>%
  filter(year %in% c(2007, 2022)) %>%                                           # Keep only 2007 and 2022
  dplyr::select(species, site, year, cover, CDT) %>%                            # Keep needed columns
  pivot_wider(
    names_from = year,
    values_from = cover,
    names_prefix = "cover_") %>%
  group_by(species,site)%>%
  filter(!(cover_2022==0 & cover_2007==0)) %>%                                  # Remove if cover in both 2007 and 2022 is 0 (as there is no data to calculate change)         
  mutate(
    cover_diff = cover_2022 - cover_2007)%>%
  ungroup()%>%
  mutate(
    zCDT = scale(CDT,center = TRUE, scale = TRUE),
    sd_year = attr(zCDT, "scaled:scale")
  )

# Community 

# Calculate variation of total shrub cover (2022-2007)
totcov = sp_data7 %>%                                                           # Get 2022-2007 Cover
  group_by (site,year) %>%
  summarise(totcov = sum(cover)) %>% 
  filter(year %in% c(2007,2022)) %>%
  pivot_wider(
    names_from = year, 
    values_from = totcov, 
    names_prefix = "totcov") %>%
  mutate(vartotcov2207 = totcov2022-totcov2007)
hist (totcov$vartotcov2207)

# Calculate CWCDT (community weighted mean centroid distance)
View(sp_data7)
meancov = sp_data7 %>%                                                          # Get mean cover by species and site (across years) 
  group_by (site, species) %>%
  summarise(meancov = mean(cover, na.rm=TRUE))

com_diff = merge(cover_diff,meancov, by=c("species","site")) %>%
  mutate(temp1 = CDT * meancov)%>%
  group_by(site) %>%
  mutate(
    covert = sum(meancov)) %>%
  summarise(CWCDT = sum(temp1)/covert) %>%
  ungroup() %>%
  distinct(site, .keep_all = TRUE) %>%
  mutate(
    zCWCDT = scale(CWCDT,center = TRUE, scale = TRUE),
    sd_year = attr(zCWCDT, "scaled:scale")
    )

com_diff = merge(com_diff,totcov, by=c("site"))

### 2.4.2. community cover diff(07-22) ~ CWCDT ####
hist(com_diff$vartotcov2207)

m_covdiw = glmmTMB(vartotcov2207 ~ zCWCDT,family = gaussian,data = com_diff)

# FINAL MODEL
m_covdif = m_covdiw
summary_model = summary(m_covdif)
summary_model
r.squaredGLMM(m_covdif)

# CHECK MODEL
res_m_covdif = simulateResiduals(m_covdif, n = 1000)                            # Get residuals
plot(res_m_covdif)                                                              # Plot for visual diagnostics
testUniformity(res_m_covdif)                                                    # Uniformity test (should be non-significant)
testDispersion(res_m_covdif)                                                    # Dispersion test (should be non-significant)
testOutliers(res_m_covdif)

# PLOT MODEL 
pm_covdiw = ggplot() +
  geom_point(data = com_diff, aes(x = CWCDT, y = vartotcov2207), color="grey30",
             size = 1.5, alpha = 0.7) +
  theme_classic() +theme(
    plot.tag = element_text(size = 11, hjust = -0.1, vjust = 1.0),plot.tag.position = c(-0.03, 1.03),    
    axis.text.x = element_text(size = 9, hjust = 0.5),axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 12, margin = margin(r = 4)),
    axis.title.x = element_text(size = 12, margin = margin(t = 4)),
    text = element_text(family = "sans"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.text.align = 0, axis.line = element_blank(), plot.title = element_blank(),
    plot.margin = unit(c(2, 2, 4, 2), "mm")) +
  labs(x = "CWCDT", y = expression(Delta*" Community cover"), tag = "a)") +
  scale_y_continuous(breaks = seq(-5, 55, by = 15), limits = c(-5,55))+
  scale_x_continuous(breaks = seq(0.009, 0.014, by = 0.001), limits = c(0.009, 0.0142)
  )
pm_covdiw


### 2.4.3. spp cover diff(07-22) ~ CDT ####

hist(cover_diff$cover_diff)
m_covdis0 = glmmTMB(cover_diff ~ zCDT + (1|site) + (1|species),family = gaussian,data = cover_diff)
m_covdis = glmmTMB(cover_diff ~ zCDT + (1|site),family = gaussian,data = cover_diff)
m_covdis1 = glmmTMB(cover_diff ~ zCDT, family = gaussian,data = cover_diff)
AICtab(m_covdis,m_covdis1,m_covdis0, base=TRUE, weights=TRUE, logLik=TRUE)

# FINAL MODEL
m_covdisf = m_covdis0
summary_model = summary(m_covdisf)
summary_model
r.squaredGLMM(m_covdisf)

# CHECK MODEL
res_m_covdisf = simulateResiduals(m_covdisf, n = 1000)
plot(res_m_covdisf)
testUniformity(res_m_covdisf)                                                   # Should be non-significant (p > 0.05)
testDispersion(res_m_covdisf)                                                   # Should be non-significant (p > 0.05)
testOutliers(res_m_covdisf)

# PLOT MODEL 
m_covdisf_df=ggpredict(m_covdisf, terms = c("zCDT"))
center = mean(unique(cover_diff$CDT))                                           # 2007–2022 mean
scale_ = sd(unique(cover_diff$CDT))  
m_covdisf_df$CDT = m_covdisf_df$x * scale_ + center

pm_covdisfsp = ggplot(m_covdisf_df, aes(CDT, predicted)) +
  geom_point(data = cover_diff, aes(x = CDT, y = cover_diff, color=species),size = 1.5, alpha = 0.7) +
  geom_line(color = "grey30", size = 0.8) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "grey30", alpha = 0.3) +  
  scale_color_manual(values = species_colors,labels = species_labels,name = "Species") +
  theme_classic() +theme(
    plot.tag = element_text(size = 11, hjust = -0.1, vjust = 1.0),plot.tag.position = c(-0.03, 1.03),    
    legend.title = element_text(size = 10),
    legend.key.size = unit(0.9, "lines"),                                     # smaller key boxes
    legend.spacing.x = unit(0.55, 'cm'), legend.spacing.y = unit(0.45, 'cm'), # spacing
    legend.margin = margin(4, 6, 4, 4),                                       # legend margins
    legend.position = "right",
    legend.background = element_rect(color = "grey40", size = 0.3, fill = NA),
    legend.text = element_text(size = 10,margin = margin(l = 1.7)),
    axis.text.x = element_text(size = 9, hjust = 0.5),axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 12, margin = margin(r = 4)),
    axis.title.x = element_text(size = 12, margin = margin(t = 4)),
    text = element_text(family = "sans"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.text.align = 0, 
    axis.line = element_blank(), plot.title = element_blank(),
    plot.margin = unit(c(2, 2, 4, 2), "mm")) + 
  labs(x = "CDT", y = expression(Delta*" Species cover"), tag = "b)") +
  scale_x_continuous(breaks = seq(0.005, 0.017, by = 0.003), limits = c(0.0046, 0.0174)
                     )
pm_covdisfsp

### 2.4.4. GRAPH -> Figure 4 #####

# Manage legends
pm_covdiw_noleg = pm_covdiw + theme(legend.position = "none")
pm_covdisfsp_noleg = pm_covdisfsp + theme(legend.position = "none")
legend_cov = cowplot::get_legend(pm_covdisfsp + theme(legend.position = "bottom"))
empty = cowplot::ggdraw() 

top_row = cowplot::plot_grid(
  pm_covdiw_noleg,
  empty,
  pm_covdisfsp_noleg,
  empty,
  ncol = 4,
  rel_widths = c(1,0.03,1,0.02),
  align = "hv"
)

combined_plot = cowplot::plot_grid(
  top_row,
  legend_cov,
  ncol = 1,
  rel_heights = c(1, 0.15)
)

combined_plot

# Add bottom margin after combining
combined_plot_with_margin = combined_plot +
  theme(plot.margin = margin(t = 2, r = 2, b = 5, l = 2))                       # Large bottom margin

#save data
pdf(file = "Figures/Figure4.ShrubCD.pdf",
    width = 7, height = 3.5)

combined_plot_with_margin
dev.off()

png(file = "Figures/Figure4.ShrubCD.png",
    res = 1000,                                                          
    width = 7, height = 3.5,                                                   
    units = "in",                                                           
    bg = "white")     

combined_plot_with_margin
dev.off()
