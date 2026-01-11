# DESeq2 contrast logic example_Additional
#
#   DESeq2 contrast logic example (illustrative)
#   This script is a minimal, *illustrative* snippet showing the DESeq2 model design and
#   contrast definitions used to generate the key differential-expression comparisons for
#   Figure 4 and the accompanying DE output tables (main effects and interaction).

#   - Defines an interaction model: ~ pollution + temp + pollution:temp
#   - Demonstrates how contrasts were extracted for:
#     (i)  warming in non-polluted conditions (black vs white)
#     (ii) pollution under ambient temperature (polluted vs non-polluted in white)
#     (iii) the interaction term (effect modification)
#     (iv) the combined treatment vs control (black+polluted vs white+non-polluted)

#   For most users, the recommended starting point is:
#   - Raw reads: ArrayExpress (Accession: E-MTAB-16541)
#   - Derived summary tables: 'Dataframes/Annotated DEG list and interaction classes.xlsx'
#     and 'Dataframes/GO enrichment results.xlsx'

####################

library("DESeq2")
library("tximport")

# import Salmon counts via tximport
dds_int <- DESeqDataSetFromTximport(genes.salmon, coldata_int, ~pollution + temp + pollution:temp)

# remove low counts (as per default pipeline)
dds_int <- dds_int[rowSums(counts(dds_int)) > 10, ]

dds_int <- DESeq(dds_int)
resultsNames(dds_int)

# extract comparisons to identify significantly differentially expressed genes (Padj<0.05)
# Warming only, i.e.comparing black vs white tiles in non-polluted areas
res_BNPvWNP <- results(dds_int, contrast=c("temp","black","white"))


# Pollution only, polluted vs non-polluted in white panels
res_WPvWNP <- results(dds_int, contrast=c("pollution","polluted","non_polluted"))

# identify significantly interactive genes - Padj < 0.05
res_interaction <- results(dds_int, name="pollutionpolluted.tempblack")


# also get log2FC for combined effect vs control, i.e. black-polluted vs white-non-polluted.
res_BPvWNP <- results(dds_int, list(c("temp_black_vs_white", "pollution_polluted_vs_non_polluted", "pollutionpolluted.tempblack")))