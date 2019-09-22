###### Calculate TMB and create plots for Figures 4 and S4########
#
#     Authors: Jo Lynne Rokita, Gregory P. Way
#     Updated 2019-09-18
################################################################

###ROC CURVE
roc_file <- file.path(dataDir, "full_roc_threshold_results.tsv")
roc_df <- readr::read_tsv(roc_file,
                          col_types = readr::cols(
                            .default = readr::col_double(),
                            gene = readr::col_character(),
                            shuffled = readr::col_character()
                          )
)

curve_colors <- c(
  "#313695",
  "#35978f",
  "#c51b7d",
  "#313695",
  "#35978f",
  "#c51b7d"
)

curve_labels <- c(
  "TP53 False" = "TP53 (AUROC = 0.89)",
  "Ras False" = "Ras (AUROC = 0.55)",
  "NF1 False" = "NF1 (AUROC = 0.77)",
  "TP53 True" = "TP53 Shuffled (AUROC = 0.46)",
  "Ras True" = "Ras Shuffled (AUROC = 0.55)",
  "NF1 True" = "NF1 Shuffled (AUROC = 0.43)"
)
roc_df$model_groups <- paste(roc_df$gene, roc_df$shuffled)
roc_df$model_groups <- factor(roc_df$model_groups, levels = names(curve_labels))

pdf(paste0(results.folder, Sys.Date(), "-roc.pdf", sep = ""), width = 7, height = 5)
ggplot(roc_df,
       aes(x = fpr,
           y = tpr,
           color = model_groups)) +
  coord_fixed() +
  geom_step(aes(linetype = shuffled), size = 0.7) +
  geom_segment(aes(x = 0,
                   y = 0,
                   xend = 1,
                   yend = 1),
               linetype = "dashed",
               color = "black") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  scale_color_manual(name = "Classifier",
                     values = curve_colors,
                     labels = curve_labels) +
  scale_linetype_manual(name = "Data",
                        values = c("solid",
                                   "dashed"),
                        labels = c("True" = "Shuffled",
                                   "False" = "Real")) +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  theme_Publication()
dev.off()



##read in variants
variants <- read.delim(paste0(dataDir, "ras-tp53-nf1-alterations.txt"), sep = "\t", header = T, as.is = T)
##subset tp53 variants
p53.var <- subset(variants, classifier == "TP53")

##load and merge with clinical file
clin = read.delim(paste0(dataDir, "pptc-pdx-clinical-web.txt"), sep = "\t", header = T, as.is = T)
##select only those models with RNA-Seq
rna <- subset(clin, RNA.Part.of.PPTC == "yes")
##read in scores
scores <- read.delim(paste0(dataDir, "classifier_scores.tsv"),
                     sep = "\t", header = T, as.is = T)
colnames(scores)[colnames(scores) == "sample_id"] <- "Model"
scores.clin <- merge(scores, clin[,c("Model", "Histology.Detailed")], all.x = T)

##select only tp53 scores to work with
tp53scores <- scores.clin[,c("Model", "Histology.Detailed", "tp53_score")]

###add tp53 variants to tp53 scores df
var.score <- merge(tp53scores, p53.var[,c("Model", "Hugo_Symbol", "classifier", "Variant_Classification", "cDNA_Change", "Protein_Change", "Fusion")], all.x = T)

##add classification to WT samples
var.score$Variant_Classification <- ifelse(is.na(var.score$Hugo_Symbol), "wild-type", var.score$Variant_Classification)
var.score$Hugo_Symbol <- ifelse(is.na(var.score$Hugo_Symbol), "wild-type", var.score$Hugo_Symbol)
var.score$classifier <- ifelse(is.na(var.score$classifier), "unclassified", var.score$classifier)

###reorder
var.score$Hugo_Symbol <- factor(var.score$Hugo_Symbol, levels = c("wild-type", "TP53", "CDKN2A", "MDM2", "MDM4", "GORAB", "ATM", "ATR", "RB1", "CHEK1", "CHEK2"))

###plot for figure 4D
pdf(paste(results.folder, Sys.Date(), "-tp53-scoresbyvarianttype.pdf", sep = ""), width = 20, height = 3)
ggplot(var.score, aes(x = Variant_Classification, y = tp53_score, color = Variant_Classification, alpha = 0.5)) + 
  geom_boxplot(outlier.shape = 21, fill = 'white') +
  geom_point(pch = 21, position = position_jitterdodge())+
  geom_hline(aes(yintercept = 0.5), colour = 'black', linetype = "longdash")+
  stat_boxplot(geom ='errorbar', width = 0.5) +
  facet_grid(~Hugo_Symbol, scales = "free_x")+
  guides(alpha = FALSE, fill = FALSE) + 
  theme_Publication() + xlab("") + ylab('TP53 Classifier Score') +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  ylim(c(0,1))+
  scale_color_manual(values = colores)
dev.off()

###subset only osteosarcoma scores for supplemental figures
ost <- subset(var.score, Histology.Detailed == "Osteosarcoma")
ost$Protein_Change <- ifelse(ost$Variant_Classification == "Splice_Site", "splice", ost$Protein_Change)
ost$Alteration <- ifelse(is.na(ost$Protein_Change) & is.na(ost$Fusion),
                         paste(ost$Hugo_Symbol, ost$Variant_Classification, sep = " "),
                         ifelse(is.na(ost$Protein_Change), paste0(ost$Fusion, " Fusion"), 
                                paste(ost$Hugo_Symbol, ost$Protein_Change, sep = " ")))
table(ost$Alteration)
ost.gather <- ost[,c("Model", "tp53_score", "Alteration")]

united <- ost.gather %>% 
  group_by(Model, tp53_score) %>% 
  summarise_all(funs(trimws(paste(., collapse = ', '))))


##export for table S4
write.table(united, paste0(results.folder, Sys.Date(), "-osteo-tp53-evidence-table.txt"), sep = "\t",
            quote = F, row.names = F, col.names = T)

pdf(paste(results.folder, Sys.Date(), "-osteo-tp53scores-byvarianttype.pdf", sep = ""), width = 15, height = 3)
print(ggplot(ost, aes(x = Variant_Classification, y = tp53_score, color = Variant_Classification, alpha = 0.5)) + 
        geom_boxplot(outlier.shape = 21, fill = 'white') +
        geom_point(pch = 21, position = position_jitterdodge())+
        geom_hline(aes(yintercept = 0.5), colour = 'black', linetype = "longdash")+
        stat_boxplot(geom ='errorbar', width = 0.5) +
        facet_grid(~Hugo_Symbol, scales = "free_x")+
        guides(alpha = FALSE, fill = FALSE) + 
        theme_Publication() + xlab("") + ylab('TP53 Classifier Score') +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank()) + 
        ylim(c(0,1))+
        scale_color_manual(values = colores))
dev.off()

####remove osteosarcoma to replot 
no.ost <- subset(var.score, Histology.Detailed != "Osteosarcoma")

pdf(paste(results.folder, Sys.Date(), "-no-osteo-tp53scores-byvarianttype.pdf", sep = ""), width = 20, height = 4)
print(ggplot(no.ost, aes(x = Variant_Classification, y = tp53_score, color = Variant_Classification, alpha = 0.5)) + 
        geom_boxplot(outlier.shape = 21, fill = 'white') +
        geom_point(pch = 21, position = position_jitterdodge())+
        geom_hline(aes(yintercept = 0.5), colour = 'black', linetype = "longdash")+
        facet_grid(~Hugo_Symbol, scales = "free_x")+
        stat_boxplot(geom ='errorbar', width = 0.5) +
        guides(alpha = FALSE, fill = FALSE) + 
        theme_Publication() + xlab("") + ylab('TP53 Classifier Score') +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank()) + 
        ylim(c(0,1))+
        scale_color_manual(values = colores))

dev.off()

###define CNVs
cnvs <- c("Hemizygous_Deletion", "Amplification", "Homozygous_Deletion")

####add CNV as a variant type for comparison by type
var.score$type <- ifelse(var.score$Variant_Classification %in% cnvs, "CNV", 
                        ifelse (var.score$Variant_Classification == "Fusion", "Fusion", 
                                ifelse (var.score$Variant_Classification == "wild-type", "wild-type", "SNV")))

##relevel
var.score$type <- factor(var.score$type, levels = c("wild-type", "SNV", "CNV", "Fusion"))

###plot scores by gene for figure 4C
pdf(paste(results.folder, Sys.Date(), "-tp53scores-bygene.pdf", sep = ""), width = 12, height = 5)
print(ggplot(var.score, aes(x = Hugo_Symbol, y = tp53_score, alpha = 0.5)) + 
        geom_boxplot(outlier.shape = 21, fill = 'white') + 
        geom_jitter(position=position_jitter(width=.1, height=0), shape = 21) + 
        stat_boxplot(geom ='errorbar', width = 0.5) +
        geom_hline(aes(yintercept = 0.5), colour = 'black', linetype = "longdash")+
        guides(alpha = FALSE, fill = FALSE) + 
        theme_Publication() + xlab("Gene") + ylab('TP53 Classifier Score') +
        ylim(c(0,1))+
        scale_color_manual(values = colores)+
        stat_n_text() +
        scale_fill_manual(values = colores))
dev.off()

###plot scores by type for figure S4C
pdf(paste(results.folder, Sys.Date(), "-tp53scores-bytype.pdf", sep = ""), width = 7, height = 5)
mycomps <- list(c("wild-type", "SNV"), c("wild-type", "CNV"), c("wild-type", "Fusion"))

ggviolin(var.score, x = "type", y = "tp53_score", fill = "type",alpha = 0.5,
         palette =colores, add.params = list(fill = "white"),
         shape = 1,
         add = "boxplot", 
         bxp.errorbar = T)+
  stat_compare_means(comparisons = mycomps, label = "p.format", method = "wilcox.test", ref.group = "wild-type", step.increase = 0.2) + ##posthoc p-values
  stat_compare_means(label.y = 1.7) + #global p value
  geom_hline(aes(yintercept = 0.5), colour = 'black', linetype = "longdash")+
  theme_Publication() + 
  xlab("Type of Alteration") + 
  ylab('TP53 Classifier Score') +
  ylim(c(0,1.7))


dev.off()

###scores by gene by type - figure S4D
pdf(paste(results.folder, Sys.Date(), "-tp53scores-genebytype.pdf", sep = ""), width = 12, height = 5)
print(ggplot(var.score, aes(x = type, y = tp53_score, color = type, alpha = 0.5)) + 
        geom_boxplot(outlier.shape = 21, fill = 'white') +
        geom_jitter(position=position_jitter(width=.1, height=0), shape = 21) +   
        geom_hline(aes(yintercept = 0.5), colour = 'black', linetype = "longdash")+
        facet_wrap(~Hugo_Symbol, scales = "free_x", nrow = 2)+
        stat_boxplot(geom ='errorbar', width = 0.5) +
        guides(alpha = FALSE) + 
        theme_Publication() + xlab("Type of Alteration") + ylab('TP53 Classifier Score') +
        ylim(c(0,1))+
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+ #+ 
        scale_color_manual(values = colores)+
        scale_fill_manual(values = colores))
dev.off()


###overall classifier status
var.score$status <- ifelse(var.score$Variant_Classification == "wild-type", "wild-type", "altered")
var.score$status <- factor(var.score$status, levels = c("wild-type", "altered"))

####add only tp53 alterations
var.score$p53status <- ifelse(var.score$Hugo_Symbol == "TP53" & var.score$status != "wild-type", "altered", "wild-type")
var.score$p53status <- factor(var.score$p53status, levels = c("wild-type", "altered"))
var.score.noost <- subset(var.score, Histology.Detailed != "Osteosarcoma")


####figure 4A
table(var.score$status) ###duplicate samples included due to multiple TP53 alterations, so remove those
tp53 <- subset(var.score, Hugo_Symbol == "TP53")

#unique altered models
un.tp53 <- unique(tp53[,c("Model", "status")])
###merge with all scores
un.tp53<- merge(un.tp53, rna[,c("Model", "Histology.Detailed")], all.y = T)
un.tp53$status <-ifelse(is.na(un.tp53$status), paste0("wild-type"), paste0(un.tp53$status))
un.tp53.scores <- merge(un.tp53, tp53scores[,c("Model", "tp53_score")])
table(un.tp53.scores$status) ###new Ns for figure 4A

pdf(paste(results.folder, Sys.Date(), "-all-tp53-scores.pdf", sep = ""), width = 7, height = 5)
print(ggviolin(un.tp53.scores, x = "status", y = "tp53_score", fill = "status",
               palette = colores,alpha = 0.5,
               add = "boxplot", 
               add.params = list(fill = "white"))+
        stat_compare_means(comparisons = list(c("wild-type", "altered")), label = "p.signif")+ # Add significance levels
        stat_compare_means(label.y = 1.5) + ##global p
        theme_Publication() + 
        xlab("Mutation Status") + ylab('TP53 Classifier Score'))
wilcox.test(tp53_score ~ status, un.tp53.scores)
median(as.numeric(unlist(subset(un.tp53.scores, status=="wild-type", select=tp53_score))))
median(as.numeric(unlist(subset(un.tp53.scores, status=="altered", select=tp53_score))))
dev.off()

###subset for models not osteosarcoma
new.noost <- subset(un.tp53.scores, Histology.Detailed != "Osteosarcoma")
table(new.noost$status) ##Ns


###figure for S4A
pdf(paste(results.folder, Sys.Date(), "-no-osteo-tp53-scores.pdf", sep = ""), width = 7, height = 5)
print(ggviolin(new.noost, x = "status", y = "tp53_score", fill = "status",
               palette = colores,alpha = 0.5,
               add = "boxplot", 
               add.params = list(fill = "white"))+
        stat_compare_means(comparisons = list(c("wild-type", "altered")), label = "p.signif")+ # Add significance levels
        stat_compare_means(label.y = 1.5) + ##global p
        theme_Publication() + 
        xlab("Mutation Status") + ylab('TP53 Classifier Score'))
dev.off()