
## @knitr setup, include=FALSE, cache=FALSE
#options(replace.assign=TRUE,width=60)
opts_chunk$set(fig.path='graphics/Cockell-', cache.path='cache/graphics-', fig.align='center', fig.height=3, tidy.opts=list(width.cutoff=60))

old_par = par(mfrow=c(1,1))
options(width=60)
  ##Formating for the table. Not used yet
  sn = function(x,digits) {
      if (x==0) return("0")
      ord = floor(log(abs(x),10))
      x = x / 10^ord
      if (!missing(digits)) x = format(x,digits=digits)
      if (ord==0) return(as.character(x))
      paste(x,"\\times 10^{",ord,"}",sep="")
    }

##Tables should be format on the decimal place
##Change 
format_dp = function(x, digits) {
    x = format(x, digits=digits)
    gsub("\\.", "&", x)
  }

simple_cap = function(x) {
      s = strsplit(x, " ")[[1]]
      s[1] = paste(toupper(substring(s[1], 1, 1)), substring(s[1], 2),
                               sep = "", collapse = " ")
      paste(s, collapse=" ")
  }
knit_theme$set("bw.css")





## @knitr eval=FALSE
## source("http://bioconductor.org/biocLite.R")


## @knitr eval=FALSE
## biocLite()


## @knitr eval=FALSE, tidy=FALSE
## ## From Bioconductor
## biocLite(c("ArrayExpress", "arrayQualityMetrics", "GOstats",
##            "GEOquery", "lumi", "lumiHumanAll.db",
##            "lumiHumanIDMapping", "RamiGO"))
## ## From CRAN
## install.packages("gplots")


## @knitr echo=FALSE, warning=FALSE,message=FALSE, cache=TRUE, results='hide', eval=TRUE
library("GEOquery")
library("lumi")
raw_data = lumiR(fileName="data/raw_data.txt")


## @knitr message=FALSE, warning=FALSE, cache=TRUE
library("ArrayExpress")
library("lumi")


## @knitr eval=FALSE
## ae = getAE("E-MTAB-1593", type="full")
## raw_data = lumiR(ae)


## @knitr results='hide', message=FALSE, warning=FALSE, cache=TRUE
str(raw_data)


## @knitr figure1, fig.keep='none', cache=TRUE, message=FALSE, tidy=FALSE
## Substitute "density" for one of the other plot types:
## "boxplot", "pair", "MA" or "cv"
plot(raw_data, what="density")


## @knitr message=FALSE, warning=FALSE, results='hide', cache=TRUE
vst_data = lumiT(raw_data)


## @knitr message=FALSE, warning=FALSE, results='hide', cache=TRUE
##Use robust spline normalisation (rsn)
rsn_data = lumiN(vst_data, method="rsn")


## @knitr message=FALSE, warning=FALSE, results='hide', cache=TRUE
qc_data = lumiQ(rsn_data)


## @knitr fig.keep='none', cache=TRUE, message=FALSE, tidy=FALSE
plot(qc_data)
plot(raw_data)


## @knitr SC_analysis_QC1, message=FALSE
library("arrayQualityMetrics")


## @knitr SC_analysis_QC2, eval=FALSE, tidy=FALSE
## arrayQualityMetrics(expressionset=qc_data, outdir="qc")


## @knitr SC_analysis_reload, warning=FALSE, results='hide', cache=TRUE
raw_data_post = raw_data[ ,-23]
vst_data_post = lumiT(raw_data_post)
rsn_data_post = lumiN(vst_data_post, method="rsn")
qc_data_post = lumiQ(rsn_data_post)


## @knitr SC_analysis_sampleQC1, cache=TRUE, tidy=FALSE
exprs_data = exprs(qc_data_post)
treatments = c("Ctrl", "TF1", "TF2", "TF3", "TF4", 
               "TF1TF4", "TF2TF4", "TF3TF4")
array_names = rep(treatments, each=3)[1:23]
colnames(exprs_data) = array_names


## @knitr SC_analysis_sampleQC2, cache=TRUE, tidy=FALSE
present_count = detectionCall(raw_data_post)
select_data = exprs_data[present_count > 0, ]


## @knitr cache=TRUE
nrow(select_data)/nrow(exprs_data)




## @knitr echo=FALSE,warning=FALSE, message=FALSE,
library(lumi)


## @knitr SC_analysis_IDmapping, message=FALSE
library("lumiHumanAll.db")
library("annotate")
probe_list = rownames(select_data)
nuIDs = probeID2nuID(probe_list)[ ,"nuID"]
symbol = getSYMBOL(nuIDs, "lumiHumanAll.db")
name = unlist(lookUp(nuIDs, "lumiHumanAll.db", "GENENAME"))


## @knitr 
anno_df = data.frame(ID= nuIDs, probe_list, symbol, name)


## @knitr SC_analysis_design, tidy=FALSE
library("limma")
design = model.matrix(~0 + factor(array_names, levels=treatments))
colnames(design) = treatments
num_parameters = ncol(design)
fit = lmFit(select_data, design)


## @knitr SC_analysis_contrasts, tidy=FALSE
cont_mat = makeContrasts(TF1-Ctrl, TF2-Ctrl, TF3-Ctrl, TF4-Ctrl,
                  TF1TF4-TF1-TF4+Ctrl, TF2TF4-TF2-TF4+Ctrl,
                  TF3TF4-TF3-TF4+Ctrl, levels=treatments)
fit2  = contrasts.fit(fit, contrasts=cont_mat)


## @knitr fit_contrasts
fit2 = eBayes(fit2)
fit2$genes = anno_df


## @knitr SC_analysis_DEG, results='hide', tidy=FALSE, cache=TRUE
## Filter by fold change (1.5x) and p-value (0.05) cutoffs
## Adjusted using Benjimini-Hochberg false discovery rate
topTable(fit2, coef="TF1 - Ctrl", p.value=0.05, lfc=log2(1.5))


## @knitr echo=FALSE, cache=TRUE
tt = topTable(fit2, coef="TF1 - Ctrl", p.value=0.05, lfc=log(1.5, 2))
tt_gs = tt$symbol
tt_logFC = format_dp(tt$logFC, 3)
tt_adj = sapply(tt$adj.P.Val, sn, 2)
tt_expr = format_dp(tt$AveExpr, 3)




## @knitr cache=TRUE, tidy=FALSE
gene_list = topTable(fit2, coef="TF1 - Ctrl", number=nrow(anno_df))


## @knitr figure2, fig.keep='none', cache=TRUE, tidy=FALSE
plot(gene_list$logFC, -log10(gene_list$adj.P.Val), 
  col=1+(abs(gene_list$logFC) > 1 & gene_list$adj.P.Val < 0.05))


## @knitr cache=TRUE
## p-values are adjusted for multiple testing
results = classifyTestsP(fit2, p.value=0.01, method="fdr")


## @knitr cache=TRUE, fig.keep='none', tidy=FALSE
vennDiagram(results[, c("TF1 - Ctrl", "TF4 - Ctrl", "TF1TF4 - TF1 - TF4 + Ctrl")])


## @knitr cache=TRUE, tidy=FALSE
tf1_table = topTable(fit2, coef="TF1 - Ctrl", 
                     n=length(probe_list), sort.by="logFC")
tf2_table = topTable(fit2, coef="TF2 - Ctrl", 
                     n=length(probe_list), sort.by="logFC")
tf3_table = topTable(fit2, coef="TF3 - Ctrl", 
                     n=length(probe_list), sort.by="logFC")
tf4_table = topTable(fit2, coef="TF4 - Ctrl", 
                     n=length(probe_list), sort.by="logFC")


## @knitr cache=TRUE, tidy=FALSE
all_signames = unique(c(rownames(tf1_table[1:10,]), 
                        rownames(tf2_table[1:10,]), 
                        rownames(tf3_table[1:10,]), 
                        rownames(tf4_table[1:10,])))


## @knitr cache=TRUE, tidy=FALSE
full = data.frame(tf1_table[all_signames,]$logFC, 
                  tf2_table[all_signames,]$logFC,
                  tf3_table[all_signames,]$logFC, 
                  tf4_table[all_signames,]$logFC,
                  row.names=all_signames)
colnames(full) = paste0("TF", 1:4)


## @knitr cache=TRUE
row_names = as.character(tf1_table[rownames(full),]$symbol)
## Filter out unannotated genes
hm_data = full[!is.na(row_names), ]
row_names = row_names[!is.na(row_names)]


## @knitr echo=FALSE
par(old_par)


## @knitr figure5, fig.keep='none', cache=TRUE, message=FALSE, fig.height=8, tidy=FALSE
library("gplots")
breaks = c(seq(min(full), 0, length.out=128), seq(0, max(full), length.out=128))
heatmap.2(as.matrix(hm_data), dendrogram="row", Colv=FALSE,
              col=bluered(255), key=TRUE, labRow=row_names,
              breaks=breaks, symkey=FALSE, density.info="none", 
              trace="none", cexRow=0.5, cexCol=0.75)




## @knitr message=FALSE, warning=FALSE, tidy=FALSE, cache=TRUE
library("GOstats")


## @knitr message=FALSE, warning=FALSE, tidy=FALSE, cache=TRUE
sig_values = tf1_table[tf1_table$adj.P.Val < 0.05 & 
                         abs(tf1_table$logFC) > log2(1.5), ]
sig_probes = as.character(sig_values$probe_list)


## @knitr message=FALSE, warning=FALSE, tidy=FALSE, cache=TRUE
## Map probe id to Entrez gene identifiers
entrez = unique(unlist(lookUp(nuIDs[sig_probes], 
                               "lumiHumanAll.db", "ENTREZID")))
entrez = as.character(entrez[!is.na(entrez)])


## @knitr message=FALSE, warning=FALSE, tidy=FALSE, cache=TRUE
## Determine the universe of possible entrez ids
entrez_universe = unique(unlist(
              lookUp(nuIDs, "lumiHumanAll.db", "ENTREZID")))
entrez_universe = as.character(entrez_universe[!is.na(entrez_universe)])


## @knitr message=FALSE, warning=FALSE, tidy=FALSE, cache=TRUE
params = new("GOHyperGParams",
    geneIds=entrez,
    universeGeneIds=entrez_universe,
    annotation="lumiHumanAll.db",
    ontology="BP",
    pvalueCutoff= 0.01,
    conditional=FALSE,
    testDirection="over")


## @knitr message=FALSE, warning=FALSE, tidy=FALSE, cache=TRUE
hyperg_result = hyperGTest(params)


## @knitr cache=TRUE
print(hyperg_result)


## @knitr message=FALSE, warning=FALSE, tidy=FALSE, cache=TRUE
## Adjust the p-values for multiple testing (FDR)
go_fdr = p.adjust(pvalues(hyperg_result), method="fdr")


## @knitr message=FALSE, warning=FALSE, tidy=FALSE, cache=TRUE
## Select the Go terms with adjusted p-value less than 0.01
sig_go_id = names(go_fdr[go_fdr < 0.01])
## Retrieve significant GO terms for BP (Biological process)
sig_go_term = getGOTerm(sig_go_id)[["BP"]]


## @knitr echo=FALSE, cache=TRUE
## Preparing data for the table
col1 = sig_go_id[1:5]
## Capitialise first letter
col2 = as.vector(sapply(sig_go_term[sig_go_id[1:5]], simple_cap))

## Put univerise in brackets
col3 = summary(hyperg_result)[1:5,]$Count
col3 = paste(col3, "(")
col3 = paste0(col3, summary(hyperg_result)[1:5,]$Size)
col3 = paste0(col3, ")")
##Format scientific 
col4 = sapply(go_fdr[sig_go_id[1:5]], sn, 2)


## @knitr eval=FALSE, cache=TRUE
## ## Visualize the enriched GO categories
## library("RamiGO")
## amigo_tree = getAmigoTree(sig_go_id)




## @knitr echo=FALSE, cache=FALSE, message=FALSE
source("http://bioconductor.org/biocLite.R")



