# APOC1

# rm(list=ls())

# chr19:45,417,504-45,422,606
# (GRCh37/hg19 by Ensembl) 

library(data.table)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(ieugwasr)
library(locuszoomr)
library(MendelianRandomization)
library(MRPRESSO)
library(patchwork)

# if(0){
genechr <- 19
genestart <- 45417504
genestop <- 45422606
genewindow <- 1000

d_tc <- fread("data/TC_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz",
               check.names = T)
d_tc <- d_tc[CHROM %in% genechr &
                 POS_b37 > genestart - genewindow*1000 &
                 POS_b37 < genestop + genewindow*1000,]

d_tc[,"p" := as.numeric(pvalue)]

d_tc <- d_tc[!(is.na(rsID)),]

loc <- locus(data = d_tc,
             gene = "APOC1", ens_db = "EnsDb.Hsapiens.v75",
             yvar = "pvalue_neg_log10",
             flank = 0)
locus_plot(loc)

loc_20kb <- locus(data = d_tc,
             gene = "APOC1", ens_db = "EnsDb.Hsapiens.v75",
             yvar = "pvalue_neg_log10",
             flank = 2e4)
locus_plot(loc_20kb)


d_ad <- fread("data/AD_sumstats_Jansenetal_2019sept.txt.gz",
              check.names = T)
d_ad <- d_ad[CHR %in% genechr &
               BP > genestart - genewindow*1000 &
               BP < genestop + genewindow*1000,]

d_ad[,"minuslog10p" := -(log(2) + pnorm(abs(Z), lower.tail = FALSE, log = TRUE))/log(10)]

d_ad[,"pval" := ifelse(P == 0, .Machine$double.xmin*2, P)]
d_ad[,"pval" := ifelse(pval < .Machine$double.xmin*2, .Machine$double.xmin*2, pval)]


loc2 <- locus(data = d_ad,
              gene = "APOC1", ens_db = "EnsDb.Hsapiens.v75",
              pos = "BP",
              yvar = "minuslog10p",
              flank = 0)
locus_plot(loc2)

loc2_20kb <- locus(data = d_ad,
              gene = "APOC1", ens_db = "EnsDb.Hsapiens.v75",
              pos = "BP",
              yvar = "minuslog10p",
              flank = 2e4)
locus_plot(loc2_20kb)




d_tc0 <- d_tc[rsID %in% loc$data$rsID,]

d_tc0[,"p1" := ifelse(p == 0, .Machine$double.xmin*2, p)]
d_tc0[,"p1" := ifelse(p1 < .Machine$double.xmin*2, .Machine$double.xmin*2, p1)]

d_top <- d_tc0[pvalue_neg_log10 > -log10(5e-8),]
# d_top <- d_ldl_0[rsid %in% d_cad$rsid,]

d_clumped2 <- ld_clump(dat = data.frame(rsid = d_top$rsID, pval = d_top$p1),
                       clump_r2 = 0.15,
                       bfile = "data/EUR", # downloaded from http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz
                       plink_bin = plinkbinr::get_plink_exe())

d_exposure <- d_tc[rsID %in% d_clumped2$rsid,]
d_outcome <- d_ad[SNP %in% d_clumped2$rsid,]

intersect(names(d_exposure), names(d_outcome))

setnames(d_outcome, "SE", "SE_outcome")

d_exposure <- d_exposure[,c("rsID", "ALT", "REF", "EFFECT_SIZE", "SE", "p")]
d_outcome <- d_outcome[,c("SNP", "A1", "A2", "BETA", "SE_outcome", "P")]
setnames(d_outcome, c("rsid", "EA_outcome", "NEA_outcome", "BETA_outcome", "SE_outcome", "P_outcome"))

d <- d_exposure[d_outcome, on = c("rsID" = "rsid")]

all.equal(d$ALT, d$EA_outcome)

# save.image(file = "APOC1_example.RData")
# }

# load("APOC1_example.RData")

d[,"bx" := abs(EFFECT_SIZE)]
d[,"by" := sign(EFFECT_SIZE)*BETA_outcome]



input_data <- mr_input(bx = d$bx, bxse = d$SE, by = d$by, byse = d$SE_outcome)

res_mr <- mr_ivw(input_data)
print(res_mr)

loo0 <- mr_loo(input_data)
loo0$data$tmp <- ifelse(loo0$data$snps == "IVW estimate", 1, 0)

p_loo <- ggplot(data = loo0$data, aes(x = estimates, y = snps, shape = factor(tmp))) +
  geom_point(size = 1) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.005) +
  geom_vline(xintercept = 0, col = "grey", linetype = "dashed") +
  scale_shape_manual(values = c(20, 18)) +
  scale_x_continuous(name = "Leave-one-out\nMR estimate",
                     limits = c(-0.1, 1.6),
                     breaks = log(c(1, 2.25, 4.5)), labels = c(1, 2.25, 4.5)) +
  scale_y_discrete(name = "Excluded variant", labels = rev(c(d$rsID, "IVW estimate"))) +
  guides(shape = "none") +
  theme_classic(base_size = 7.5)

forest0 <- mr_forest(input_data)
forest0$data$tmp <- ifelse(forest0$data$snps == "IVW estimate", 1, 0)

p_forest <- ggplot(data = forest0$data, aes(x = estimates, y = snps, shape = factor(tmp))) +
  geom_point(size = 1) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.005) +
  geom_vline(xintercept = 0, col = "grey", linetype = "dashed") +
  scale_shape_manual(values = c(20, 18)) +
  scale_x_continuous(name = "MR estimate",
                     breaks = log(c(1, 3, 10, 30, 60)), labels = c(1, 3, 10, 30, 60)) +
  scale_y_discrete(name = "Variant", labels = rev(c(d$rsID, "IVW estimate"))) +
  guides(shape = "none") +
  theme_classic(base_size = 7.5)

p_scatter <- ggplot(d, aes(x = bx, y = by)) +
  geom_hline(yintercept = 0, col = "grey", linetype = "dashed") +
  geom_vline(xintercept = 0, col = "grey", linetype = "dashed") +
  geom_abline(intercept = 0, slope = res_mr@Estimate, col = "darkblue") +
  geom_point(size = 0.5) +
  geom_errorbar(aes(ymin = by - 1.96*SE_outcome, ymax = by + 1.96*SE_outcome), width = 0.005*(2.5/8)) +
  geom_errorbarh(aes(xmin = bx - 1.96*SE, xmax = bx + 1.96*SE), height = 0.005) +
  theme_classic() +
  scale_x_continuous(name = bquote(beta[GX]), limits = c(0, 0.25)) +
  scale_y_continuous(name = bquote(beta[GY]), limits = c(-0.2, 0.6))



res_median <- mr_median(input_data)
print(res_median)

res_mode <- mr_mbe(input_data)
print(res_mode)

exp(res_mr@Estimate + c(0, -1, 1)*1.96*res_mr@StdError)
exp(res_median@Estimate + c(0, -1, 1)*1.96*res_median@StdError)
exp(res_mode@Estimate + c(0, -1, 1)*1.96*res_mode@StdError)

(res_mr@Estimate - res_median@Estimate)/sqrt(res_mr@StdError^2 + res_median@StdError^2)
(res_mr@Estimate - res_mode@Estimate)/sqrt(res_mr@StdError^2 + res_mode@StdError^2)

res_presso <- mr_presso("BETA_outcome", "EFFECT_SIZE", "SE_outcome", "SE",
                        data = as.data.frame(d),
                        seed = as.numeric(paste0(utf8ToInt("APOC1"), collapse = "")) %% .Machine$integer.max)

res_egger <- mr_egger(input_data)
print(res_egger)



LDmat <- ld_matrix(c(d$rsID, "rs429358", "rs7412"),
                   bfile = "data/EUR",
                   plink_bin = plinkbinr::get_plink_exe())

LDmat[grep("rs429358", rownames(LDmat)),]^2
LDmat[grep("rs7412", rownames(LDmat)),]^2

### Figures

variants <- intersect(loc_20kb$data$rsID, loc2_20kb$data$SNP)

dt1 <- loc_20kb$data[loc_20kb$data$rsID %in% variants,]

ps1 <- ggplot(
  data = data.frame(x = dt1$POS_b37, y = dt1$pvalue_neg_log10),
  aes(x = x, y = y)) +
  geom_point(size = 1) +
  scale_y_continuous(name = bquote(-log[10](italic(p)))) +
  ggtitle("Total cholesterol levels") +
  theme_classic(base_size = 7.5) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

dt2 <- loc2_20kb$data[loc2_20kb$data$SNP %in% variants,]

ps2 <- ggplot(
  data = data.frame(x = dt2$BP, y = dt2$minuslog10p),
  aes(x = x, y = y)) +
  geom_point(size = 1) +
  scale_y_continuous(name = bquote(-log[10](italic(p)))) +
  ggtitle("Alzheimer's disease") +
  theme_classic(base_size = 7.5) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

pgenes <- gg_genetracks(loc_20kb, highlight = "APOC1",
                        cex.text = 0.5,
                        cex.axis = 0.6,
                        cex.lab = 0.7)

(p_locus <- ps1 / ps2 / pgenes +
    plot_layout(heights = c(1, 1, 1/3)))

p_scatter <- p_scatter + theme_classic(base_size = 7.5)

p_full <- (p_locus | (free(p_scatter) / (p_forest + p_loo))) +
  plot_layout(widths = c(2/(1+sqrt(5)), 1)) +
  plot_annotation(tag_levels = list(c("A", "", "", "B", "C", "D")))
p_full

ggsave("Figure3.pdf", p_full, width = 170, height = 170/2, units = "mm", dpi = 600)
