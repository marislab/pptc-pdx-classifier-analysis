###### Calculate TMB and create plots for Figures 4 and S4########
#
#     Author: Jo Lynne Rokita
#     Updated 2019-09-22
################################################################

#read in tmb data
bp <- read.delim(paste0(dataDir, "Breakpoints_rawdata.txt"), sep = "\t",
                 header = T, as.is = T)
tmb <- read.delim(paste0(dataDir, "mutations-per-model.txt"), sep = "\t",
                  header = T, as.is = T)
#read in signatures
sigs <- read.delim(paste0(dataDir, "allpdx-signature-weights.txt"),
                   sep = "\t", header = T, as.is = T)
rownames(sigs) <- sigs$Model
sigs$Model <- NULL
##remove sigs with no values and unknown sig
Summat <-  colSums(sigs,na.rm = TRUE)
nonzero.sigs <- sigs[, Summat!=0]
nonzero.sigs$Unknown <- NULL
nonzero.sigs$Model <- rownames(nonzero.sigs)

##merge breakpoints and scores
bp.scores <- merge(bp[,c("Model","n.breakpoints")], scores, all = T)
summary(lm(bp.scores$n.breakpoints ~ bp.scores$tp53_score))
##add tmb to the merge
bp.scor.tmb <- merge(bp.scores, tmb[,c("Model", "MutperMB")], all = T)
##log tmb
bp.scor.tmb$tmb <- log(bp.scor.tmb$MutperMB, 2)
###add signatures to merge
all.mer <- merge(bp.scor.tmb, nonzero.sigs, all = T)
##move model to rownames for correlation plots
rownames(bp.scor.tmb) <- bp.scor.tmb$Model

###create first cor matrix 
first <- bp.scor.tmb[,c("tp53_score", "n.breakpoints", "tp53_shuffle", "tmb")]
###remove ost models to see if correlation still holds
no.ost.models <- var.score.noost$Model
first.no.ost <- first[no.ost.models,]

cor.test(first$tp53_score, first$n.breakpoints, method = "pearson")
cor.test(first.no.ost$tp53_score, first.no.ost$n.breakpoints, method = "pearson") ##yes
cor.test(first$tmb, first$n.breakpoints, method = "pearson")
#remove model column 
cor.df <- all.mer[,2:ncol(all.mer)]

####create correlation matrix
cormat <- round(cor(first, use = "pairwise", method = "pearson"),2)
source(paste0(script.folder, "cormat.R"))

# Create a ggheatmap for figure 4E
pdf(paste0(results.folder, Sys.Date(), "-tp53-breakpoint-corplot.pdf"), height = 6, width = 6)
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "#005789", high = "#D55E00", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_Publication()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
# Print the heatmap with labels
print(ggheatmap + 
        geom_text(aes(Var2, Var1, label = value), color = "black", size = 6) +
        theme_Publication()+
        theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                         size = 12, hjust = 1))+
        labs(y ="", x = ""))
dev.off()


#######second cor matrix: - sig 9 not in cor.df anymore
second <- cor.df[,c("tp53_score", "n.breakpoints", "tp53_shuffle", "nf1_score", "nf1_shuffle", 
                    "Signature.1", "Signature.13", "Signature.2", "Signature.15", "tmb",
                    "Signature.3", "Signature.6", "Signature.20", "Signature.26", "Signature.9",
                    "Signature.10", "Signature.21")]
cor.test(second$Signature.1, second$Signature.6, method = "pearson")
cor.test(second$Signature.1, second$Signature.3, method = "pearson")
cor.test(second$Signature.2, second$Signature.13, method = "pearson")

cormat <- round(cor(second, use = "pairwise", method = "pearson"),2)
source(paste0(script.folder, "cormat.R"))

###figure S4E
pdf(paste0(results.folder, Sys.Date(), "-sig-breakpoint-corplot.pdf"), height = 10, width = 10)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "#005789", high = "#D55E00", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_Publication()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
# Print the heatmap
print(ggheatmap + 
        geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
        theme_Publication()+
        theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                         size = 12, hjust = 1))+
        labs(y ="", x = ""))
dev.off()

