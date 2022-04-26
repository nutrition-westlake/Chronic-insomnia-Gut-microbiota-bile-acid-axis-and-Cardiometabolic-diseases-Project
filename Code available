16S rRNA sequencing profiling
Analysis of chronic insomnia with the diversity of the gut microbiome 
Analysis of chronic insomnia with gut microbiome biomarkers
Analysis of chronic insomnia with bile acid biomarkers
Analysis of association between gut microbiome biomarkers and bile acid biomarkers
Analysis of chronic insomnia-related gut microbial features and bile acids with cardiometabolic diseases and risk factors
Metaanalysis of association between chronic insomnia-related gut microbial features with cardiometabolic diseases
Analysis of prospective association of chronic insomnia-related gut microbial features with CMD
Mediation of bile acids for association between gut microbiota and cardiometabolic diseases
Analysis of prospective association between habitual dietary intakes and features of gut microbiota and bile acids

16S rRNA sequencing profiling
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.txt \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
  --i-data paired-end-demux.qza \
  --o-visualization demux.qzv

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-end-demux.qza \
  --p-trim-left-f 26 \
  --p-trim-left-r 26 \
  --p-trunc-len-f 247 \
  --p-trunc-len-r 247 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza \
  --p-n-threads 8

qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file metadata.txt

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza \
  --p-n-threads 8

qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 10000 \
  --m-metadata-file metadata.txt \
  --o-visualization alpha-rarefaction.qzv

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 6000 \
  --m-metadata-file metadata.txt \
  --output-dir core-metrics-results

qiime feature-table filter-features \
  --i-table table.qza \
  --p-min-samples 2 \
  --o-filtered-table filtered-table.qza

qiime feature-classifier classify-sklearn \
  --i-classifier classifier_silva.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza


qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file metadata.txt \
  --o-visualization taxa-bar-plots.qzv

Chronic insomnia was associated with the diversity of the gut microbiome

In the Guangzhou Nutrition and Health Study (GNHS), we examined the association of chronic insomnia with gut microbial α-diversity indices (Observed species, Shannon index and Chao 1 index)
among the four groups using a multivariable linear regression with three different statistical models. 
Model 1 was adjusted for age, sex, BMI, smoking status, alcohol status, physical activity, education, income and total energy intake at baseline. 
Model 2 was additionally controlled for hypertension, hyperlipidemia, MetS, T2D, CHD, stroke and medication for T2D. 
Model 3 was further adjusted for dietary intake of vegetables, fruits, red and processed meat, fish, dairy products, coffee and tea. 

In the GNHS, there was a significant difference in the gut microbiome structure for the New-onset group or Long-term chronic insomnia group, compared with the Long-term healthy group. 
To identify robust microbial biomarkers of chronic insomnia and increase the sample size, we combined the New-onset group and Long-term chronic insomnia group into Chronic insomnia group.

We also validated the results in an independent large cross-sectional study, the Guangdong Gut Microbiome Project (GGMP).

R code, for example:

Q<-read.table("α_diversity.txt",sep="\t",header=T,row.names=1)

hist(Q$Observed_species)

Model 1
Model_1<- lm(Observed_species~Chronic_insomnia_group+age+sex+BMI+smoking_status+alcohol_status+physical_activity+education+income+total_energy_intake,data=Q)
summary(Model_1)

Model 2
Model_2<- lm(Observed_species~Chronic_insomnia_group+age+sex+BMI+smoking_status+alcohol_status+physical_activity+education+income+total_energy_intake+hypertension+hyperlipidemia+MetS+T2D+CHD+stroke+medication_for_T2D,data=Q)
summary(Model_2)

Model 3
Model_3<- lm(Observed_species~Chronic_insomnia_group+age+sex+BMI+smoking_status+alcohol_status+physical_activity+education+income+total_energy_intake+hypertension+hyperlipidemia+MetS+T2D+CHD+stroke+medication_for_T2D+vegetables_intake+fruits_intake+red_and_processed_meat_intake+fish_intake+dairy_products_intake+coffee+tea,data=Q)
summary(Model_3)

The association between chronic insomnia β-diversity dissimilarity based on genus-level Bray-Curtis distance was examined using permutational ANOVA (PERMANOVA) (999 permutations). 

R code, for example:
otu <- read.delim('β_diversity.txt', row.names = 1, sep = '\t', stringsAsFactors = F, check.names = F)
distance <- vegdist(otu, method = 'bray')
distance <- as.matrix(distance)
write.table(cbind(rownames(distance), distance), 'distance.txt', row.names = F, sep = '\t', quote = F)
DT<-read.table("distance.txt",sep="\t",header=T,row.names=1)
group<-read.table("Group.txt",sep="\t",header=T,row.names=1)

Model 1
adonis(DT~group$Chronic_insomnia_group+group$age+group$sex+group$BMI+group$smoking_status+group$alcohol_status+group$physical_activity+group$education+group$income+group$total_energy_intake, data = group, permutations = 999,method="bray")
Model 2
adonis(DT~group$Chronic_insomnia_group+group$age+group$sex+group$BMI+group$smoking_status+group$alcohol_status+group$physical_activity+group$education+group$income+group$total_energy_intake+group$hypertension+group$hyperlipidemia+group$MetS+group$T2D+group$CHD+group$stroke+group$medication_for_T2D, data = group, permutations = 999,method="bray")
Model 3
adonis(DT~group$Chronic_insomnia_group+group$age+group$sex+group$BMI+group$smoking_status+group$alcohol_status+group$physical_activity+group$education+group$income+group$total_energy_intake+group$hypertension+group$hyperlipidemia+group$MetS+group$T2D+group$CHD+group$stroke+group$medication_for_T2D+group$vegetables_intake+group$fruits_intake+group$red_and_processed_meat_intake+group$fish_intake+dairy_products_intake+group$coffee+group$tea, data = group, permutations = 999,method="bray")


Figure 1b, Figure 2a-b and Supplementary Fig. 1-3, R code, for example:

α_diversity:

Q<-read.table("Diversity.txt",sep="\t",header=T,row.names=1)
ggboxplot(Q, "Group","Observed_species", color = "Group", palette =c('orchid1', 'tan1', 'cornflowerblue', 'springgreen3'), add = "jitter", size=0.5, add.params = list(size=0.5))

β_diversity:

dis <- read.delim('distance.txt', row.names = 1, sep = '\t', stringsAsFactors = F, check.names = F)
group <- read.delim('Group.txt', sep = '\t', stringsAsFactors = F)
pcoa <- cmdscale(as.dist(dis), k = 6, eig = T)
ordiplot(scores(pcoa)[ ,c(1, 2)], type = 't')
summary(pcoa)
pcoa$eig
point <- data.frame(pcoa$point)
write.csv(point, 'pcoa.sample_insomnia_model.csv')
species <- wascores(pcoa$points[,1:2]，)
pcoa_eig <- {pcoa$eig}[1:2]
sample_site <- data.frame({pcoa$point})[1:2]
sample_site$id <- rownames(sample_site)
names(sample_site)[1:2] <- c('PCoA1', 'PCoA2')
sample_site <- merge(sample_site, group, by = 'id', all.x = T)
sample_site$Group <- factor(sample_site$Group, levels = c('Long-term healthy group', 'Recovery group','New-onset group', 'Long-term chronic insomnia group'))
group_border <- ddply(sample_site, 'Group', function(df) df[chull(df[[2]], df[[3]]), ])
pcoa_plot <- ggplot(sample_site, aes(PCoA1, PCoA2, group = Group)) +
theme(panel.grid = element_line(color = 'gray', linetype = 2, size = 0.1), panel.background = element_rect(color = 'black', fill = 'transparent'), legend.key = element_rect(fill = 'transparent')) + 
geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
geom_point(aes(color = Group), size = 1.0, alpha = 0.8) +
stat_ellipse(level = 0.8, aes(color = Group)) +
scale_shape_manual(values = c(17, 24)) + 
scale_color_manual(values = c('orchid1', 'tan1', 'cornflowerblue', 'springgreen3')) +
guides(fill = guide_legend(order = 1), shape = guide_legend(order = 2), color = guide_legend(order = 3)) + 
labs(x = paste('PCoA1: ', round(pcoa_eig[1], 2), '%'), y = paste('PCoA2: ', round(pcoa_eig[2], 2), '%'))
ggsave('PCoA.png', pcoa_plot, width = 6, height = 5)


Analysis of Chronic insomnia with gut microbiome biomarkers

We used linear regression, adjusted for the same covariates as above model 3, to confirm the association of chronic insomnia with the gut microbiome biomarkers.
The Benjamini-Hochberg method was used to control the false discovery rate (FDR).

R code, for example:
Q<-read.table("gut_microbiome_biomarkers.txt",sep="\t",header=T,row.names=1)
for(i in 1:x) {(A<- lm(gut_microbiome_biomarker~Chronic_insomnia_group+age+sex+BMI+smoking_status+alcohol_status+physical_activity+
education+income+total_energy_intake+hypertension+hyperlipidemia+MetS+T2D+CHD+stroke+medication_for_T2D+
vegetables_intake+fruits_intake+red_and_processed_meat_intake+fish_intake+dairy_products_intake+coffee+tea,data=Q))
D[i]<-summary(A)$coefficients[2,4]
F[i]<-summary(A)$coefficients[2,1]}
write.csv(D, 'p_value.csv')
write.csv(F, 'beta_value.csv')


Analysis of Chronic insomnia with bile acids

We used linear regression, adjusted for the same covariates as above model 3, to confirm the association of chronic insomnia with the bile acid biomarkers.
The Benjamini-Hochberg method was used to control the false discovery rate (FDR).

R code, for example:
Q<-read.table("bile_acid_biomarkers.txt",sep="\t",header=T,row.names=1)
for(i in 1:x) {(A<- lm(gut_microbiome_biomarker~Chronic_insomnia_group+age+sex+BMI+smoking_status+alcohol_status+physical_activity+
education+income+total_energy_intake+hypertension+hyperlipidemia+MetS+T2D+CHD+stroke+medication_for_T2D+
vegetables_intake+fruits_intake+red_and_processed_meat_intake+fish_intake+dairy_products_intake+coffee+tea,data=Q))
D[i]<-summary(A)$coefficients[2,4]
F[i]<-summary(A)$coefficients[2,1]}
write.csv(D, 'p_value.csv')
write.csv(F, 'beta_value.csv')

Figure 2e, R code, for example:
Q<-read.table("bile_acid_biomarkers.txt",sep="\t",header=T,row.names=1)
data<-transform(Q, dist_cat_n=as.numeric(as.factor(class)), scat_adj=ifelse(Group == "Long-term healthy group", -0.2,0.2))
ggplot(data, aes(x=class, y=concentration))+
geom_boxplot(outlier.size=0, aes(fill=factor(Group)), position = position_dodge(0.8),size=0.4)+
scale_fill_manual(values= c( 'orchid1', 'springgreen3'))+
geom_jitter(aes(scat_adj+dist_cat_n, concentration, fill = factor(Group)), position = position_jitter(width=0.1,height=0),shape=21,size=1.5)+
guides(fill=guide_legend(title="Group"))+
theme_light()


Analysis of association between gut microbiome biomarkers and bile acid biomarkers

We examined the association of the above identified gut microbiota biomarkers with bile acid biomarkers using partial correlation analysis, adjusted for age, sex and BMI.

R code, for example:

data<-read.csv("Analysis_of_association_between_gut_microbiome_biomarkers_and_bile_acid_biomarkers.csv",header=T)
for(i in 1:x) { r[i] <-pcor(c(4,i+x,1,2,3), var(data)) }
write.csv(r, 'Pcorr_R_value_between_gut_microbita_and_bilc_acid_biomarkers.csv')

data<-read.csv("Analysis_of_association_between_gut_microbiome_biomarkers_and_bile_acid_biomarkers.csv",header=T)
for(i in 1:x) {r<-pcor(c(4,i+x,1,2,3), var(data))
pcor_test<- pcor.test(r, length(c(1,2,3)),dim(data)[1])
XT_p[i]<- data.frame(pcor_test$pvalue, check.names = FALSE)}
write.csv(XT_p, 'Pcorr_p_value_between_gut_microbita_and_bilc_acid_biomarkers.csv')



Analysis of chronic insomnia-related gut microbial features and bile acids with cardiometabolic diseases (CMD) and risk factors

Association of chronic insomnia-related gut microbial features and bile acids with CMD

We investigated the association of the chronic insomnia-related microbial and bile acid biomarkers with different CMD
using multivariable logistic regression in the GNHS, adjusted for age, sex, smoking status, alcohol status, physical activity, education, income, and total energy intake. 
The Benjamini-Hochberg method was used to control the false discovery rate (FDR).
We also repeated the results in the GGMP.

R code, for example:
Q<-read.table("Association_between_gut_Microbiota_biomarkers_and_Cardiometabolic_diseases.txt",sep="\t",header=T,row.names=1)
for(i in 1:x) {(A<-glm(Q[,i]~gut_Microbiota_biomarker+age+sex+BMI+smoking_status+alcohol_status+physical_activity+education+income+total_energy_intake,data=Q,family=binomial()))
D[i]<-summary(A)$coefficients[2,4]}
write.csv(D, 'P_value.csv')
for(i in 1:x) {(A<-glm(Q[,i]~gut_Microbiota_biomarker+age+sex+BMI+smoking_status+alcohol_status+physical_activity+education+income+total_energy_intake,data=Q,family=binomial()))
E[i]<-summary(A)$coefficients[2,1]}
write.csv(E, 'beta.csv')

Figure 3a-c, R code, for example:

x<-read.table("Odds_ratio_of_association_between_microbiota_biomarkers_and_cardiometabolic_diseases.txt",sep="\t",header=T,row.names=1)
bk <- c(0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4)
pheatmap(x,cellwidth = 70, cellheight = 30,  color = colorRampPalette(c("springgreen3","white","tan1"))(16), fontsize=6,  breaks=bk, cluster_rows = T, cluster_cols = F, cutree_rows=3)


Association of chronic insomnia-related gut microbial features and bile acids with CMD risk factors

We investigated the association of the chronic insomnia-related microbial and bile acid biomarkers with CMD risk factors 
(BMI, DBP, SBP, waist circumference, and fasting serum levels of TG, TC, HDL, LDL, glucose, insulin and HbA1c) using 
linear regression model in the GNHS, adjusted for age, sex, smoking status, alcohol status, physical activity, education, income, and total energy intake. 
The Benjamini-Hochberg method was used to control the false discovery rate (FDR).

R code, for example:
Q<-read.table("Association_between_gut_Microbiota_biomarkers_and_Cardiometabolic_disease_risk_factors.txt",sep="\t",header=T,row.names=1)
for(i in 1:x) {(A<-lm(Q[,i]~gut_Microbiota_biomarker+age+sex+BMI+smoking_status+alcohol_status+physical_activity+education+income+total_energy_intake,data=Q))
D[i]<-summary(A)$coefficients[2,4]}
write.csv(D, 'P_value.csv')
for(i in 1:x) {(A<-glm(Q[,i]~gut_Microbiota_biomarker+age+sex+BMI+smoking_status+alcohol_status+physical_activity+education+income+total_energy_intake,data=Q))
E[i]<-summary(A)$coefficients[2,1]}
write.csv(E, 'beta.csv')

Metaanalysis of association between chronic insomnia-related gut microbial features with cardiometabolic diseases

The effect estimates of the association between chronic insomnia-related gut microbial features with cardiometabolic diseases from the GNHS and the GGMP were pooled by random effects meta-analysis.

R code, for example:
Q<-read.table("Metaanalysis.txt",sep="\t",header=T,row.names=1)
metamod<-rma(yi=B,data=Q,sei=SE,method="DL")
summary(metamod)


Prospective association of chronic insomnia-related gut microbial features with CMD

We further examined the prospective association of the above identified gut microbiota biomarkers with the incidence of CMD outcomes at the third follow-up
using multivariable logistic regression, adjusting for age, sex, smoking status, alcohol status, physical activity, education, income, and total energy intake.

R code, for example:
Q<-read.table("Prospective_association_of_Chronic_insomnia_related_gut_microbial_features_with_CMD.txt",sep="\t",header=T,row.names=1)
for(i in 1:x) {(A<-glm(Q[,i]~gut_Microbiota_biomarker+age+sex+BMI+smoking_status+alcohol_status+physical_activity+education+income+total_energy_intake,data=Q,family=binomial()))
D[i]<-summary(A)$coefficients[2,4]}
write.csv(D, 'P_value.csv')
for(i in 1:x) {(A<-glm(Q[,i]~gut_Microbiota_biomarker+age+sex+BMI+smoking_status+alcohol_status+physical_activity+education+income+total_energy_intake,data=Q,family=binomial()))
E[i]<-summary(A)$coefficients[2,1]}
write.csv(E, 'beta.csv')


Mediation of bile acids for association between gut microbiota and cardiometabolic diseases

Based on the biological plausibility of the associations among the gut microbiota, bile acids and CMD, and our above findings, 
we performed mediation analysis to evaluate whether bile acids could mediate the association of the chronic insomnia related-gut microbiota with CMD outcomes 
(gut microbiota → bile acids → CMD). Sensitivity analysis was performed to test the robustness of the mediation effect.

Mediation analysis, R code, for example:

Q<-read.table("Mediation_data.txt",sep="\t",header=T,row.names=1)
f<- lm(bile_acid~gut_Microbiota_biomarker+age+sex+smoking_status+alcohol_status+physical_activity+education+income+total_energy_intake,data=Q)
d<- glm(CMD~gut_Microbiota_biomarker+bile_acid+age+sex+smoking_status+alcohol_status+physical_activity+education+income+total_energy_intake,data=Q,family=binomial())
CM<- mediate(f, d, treat = "gut_Microbiota_biomarker", mediator = "bile_acid", boot = “TRUE”, boot.ci.type =“perc”, conf.level = 0.95, sims=1000)
summary(CM)

Mediation sensitivity analysis, R code, for example:

WD<-read.table("Mediation_data.txt",sep="\t",header=T,row.names=1)
f<- lm(bile_acid~gut_Microbiota_biomarker+age+sex+smoking_status+alcohol_status+physical_activity+education+income+total_energy_intake,data=WD)
d<- glm(CMD~gut_Microbiota_biomarker+bile_acid+age+sex+smoking_status+alcohol_status+physical_activity+education+income+total_energy_intake,data=WD,family=binomial(probit))
SA<- mediate(f, d, treat = "gut_Microbiota_biomarker", mediator = "bile_acid", boot = TRUE, boot.ci.type = "perc", conf.level = 0.95, sims = 1000)
sens.cont<- medsens(SA, rho.by=0.01, eps=.01, effect.type="indirect", sims=1000)
summary(sens.cont)


Analysis of prospective association between habitual dietary intakes and features of gut microbiota and bile acids

we used a linear regression model to determine the prospective association of dietary factors with the gut microbial and bile acid mediators of chronic insomnia and CMD, 
adjusted for age, sex, BMI, smoking status, alcohol status, physical activity, education, income, dietary intake of vegetables/fruits/red and processed meat/fish/dairy products/coffee/tea)
 (mutual adjustment for each other) and total energy intake. 
The analyses were conducted among the GNHS participants without chronic insomnia or CMD at baseline.
We also repeated the results in the GGMP.

R code, for example:
Q<-read.table("Prospective_association_of_habitual_dietary_intakes_and_features_of_gut_microbiota_and_bile_acids.txt",sep="\t",header=T,row.names=1)
A<-lm(gut_Microbiota_biomarker~vegetables_intake+fruits_intake+red_and_processed_meat_intake+fish_intake+dairy_products_intake+coffee+tea+age+sex+BMI+smoking_status+alcohol_status+physical_activity+education+income+total_energy_intake,data=Q,family=binomial()))
summary(A)

Figure 4a-b, R code, for example:
x<-read.table("Beta_coefficient_association_of_habitual_dietary_intakes_and_features_of_gut_microbiota_and_bile_acids.txt",sep="\t",header=T,row.names=1)
bk <- c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3)
pheatmap(x,cellwidth = 70, cellheight = 30,  color = colorRampPalette(c("springgreen3","white","tan1"))(6), fontsize=6,  breaks=bk, cluster_rows = F, cluster_cols = F)
