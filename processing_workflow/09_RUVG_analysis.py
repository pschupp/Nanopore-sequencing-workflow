import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.stats import zscore

# f(x) input
project_name = "nano_24_sec_seq"
expr_counts = pd.read_csv("nano_24_sec_min.tsv", sep="\t")
expr_counts = expr_counts.drop(index=range(5682, 5690))
ercc_counts = pd.read_csv("~/Documents/24_sec_seq/root/06_alignment_to_genome/minimap2/GRCh37.p7_alignment/ERCC_alignment/ercc.counts.tsv", sep="\t")
oldham_genes = pd.read_csv("Oldham_fidelity_4_celltypes_genes.csv", sep=",")
oldham_genes = oldham_genes.iloc[:, 1:]

# Processing with only numbers in actual dataframe for easier numerical processing later on. rownames are genes/ercc, colnames are samples
expr_counts.set_index(expr_counts.columns[0], inplace=True)
genes = expr_counts.index
expr_counts = expr_counts.apply(pd.to_numeric)

expr_col_sums = expr_counts.sum(axis=0)
expr_norm = expr_counts.div(expr_col_sums, axis=1)
expr_norm_l2 = np.log2(expr_norm + 1)
expr_norm_l2.columns = expr_counts.columns
expr_norm_l2.index = genes

ercc_counts.set_index(ercc_counts.columns[0], inplace=True)
ercc_counts = ercc_counts.iloc[:, 2:]
ercc_counts = ercc_counts.drop(columns=[ercc_counts.columns[7], ercc_counts.columns[15]])
ercc_counts.columns = expr_counts.columns
ercc_counts = ercc_counts[(ercc_counts != 0).sum(axis=1) <= 12]  # remove ercc synthetic transcripts for which there is a dearth of count data

# Using RUVg we are measuring two metrics to validate our approach.
# First we will be using PCA analysis to look at how clustering improves through RUVg normalization.
# Second we will use the log2 ratio of gene count / median gene count over all sections.

# First we will plot these plots with baseline data, log2 normalized and not log2 normalized.
# Note: RUVg function needs to be defined or imported from a library

# expr_ruvg = RUVg(x=expr_counts.values, cIdx=ercc_counts.values, k=1, isLog=False)
# expr_ruvg_df = pd.DataFrame(expr_ruvg['normalizedCounts'])
# expr_ruvg_pca = PCA()
# expr_ruvg_pca.fit(expr_ruvg_df.T)
# expr_ruvg_pca_eig = expr_ruvg_pca.explained_variance_
# expr_ruvg_pca_pvar = expr_ruvg_pca_eigenvalues / np.sum(expr_ruvg_pca_eigenvalues)
# xlab_char = f"PC1 {np.round(expr_ruvg_pca_pvar[0] * 100, 3)} % var. explained"
# ylab_char = f"PC2 {np.round(expr_ruvg_pca_pvar[1] * 100, 3)} % var. explained"

# plt.figure()
# plt.title("Plot of First 2 PCs of Log2-Norm Data")
# plt.xlabel(xlab_char)
# plt.ylabel(ylab_char)
# plt.scatter(expr_ruvg_pca.transform(expr_ruvg_df.T)[:, 0], expr_ruvg_pca.transform(expr_ruvg_df.T)[:, 1])
# plt.savefig(f"{project_name}_pca_count_k0.pdf")
# plt.close()

# Similar plotting for expr.ruvg.pca.vec_count_k0

expr_rmedian = expr_counts.median(axis=1).values
expr_counts_r = pd.DataFrame(index=expr_counts.index, columns=expr_counts.columns)

for row in range(expr_counts.shape[0]):
    expr_counts_r.iloc[row, :] = expr_counts.iloc[row, :] + 1 / expr_rmedian[row]

expr_counts_r_l2 = np.log2(expr_counts_r)
expr_counts_r_0 = expr_counts_r.copy()
expr_counts_r[expr_counts_r == 0] = np.nan
expr_counts_r_l2[expr_counts_r_l2 == 0] = np.nan

plt.figure(figsize=(18, 6))
plt.boxplot(expr_counts_r_l2)
plt.title("RLE of Unnormalized Data")
plt.ylabel("RLE log2(count/median count)")
plt.xlabel("Samples")
plt.savefig(f"{project_name}_rle_count_k0.pdf")
plt.close()

# Similar plotting for RE of Unnormalized Data and RLE of Unnormalized Data with incremented zeros

k_ruv = 0
expr_counts.to_csv(f"{project_name}_expression_matrix_log2_norm_k{k_ruv}.tsv", sep="\t", index=True)

# Now doing iterative RUVg with various k's
for k_ruv in range(1, ercc_counts.shape[0] + 1):
    # expr_ruvg = RUVg(x=expr_counts.values, cIdx=ercc_counts.values, k=k_ruv, isLog=False)
    # expr_ruvg_df = pd.DataFrame(expr_ruvg['normalizedCounts'])
    # expr_ruvg_df = np.log2(expr_ruvg_df + 1)

    # expr_ruvg_df.to_csv(f"{project_name}_expression_matrix_log2_norm_k{k_ruv}.tsv", sep="\t", index=True)

    vars = expr_ruvg_df.var(axis=1)
    expr_ruvg_df = expr_ruvg_df[vars != 0]
    if expr_ruvg_df.shape[0] > 0:
        expr_ruvg_pca = PCA()
        expr_ruvg_pca.fit(expr_ruvg_df.T)
        expr_ruvg_pca_eig = expr_ruvg_pca.explained_variance_
        expr_ruvg_pca_pvar = expr_ruvg_pca_eigenvalues / np.sum(expr_ruvg_pca_eigenvalues)
        xlab_char = f"PC1 {np.round(expr_ruvg_pca_pvar[0] * 100, 3)} % var. explained"
        ylab_char = f"PC2 {np.round(expr_ruvg_pca_pvar[1] * 100, 3)} % var. explained"

        plt.figure()
        plt.title("Plot of First 2 PCs of Log2-Norm Data")
        plt.xlabel(xlab_char)
        plt.ylabel(ylab_char)
        plt.scatter(expr_ruvg_pca.transform(expr_ruvg_df.T)[:, 0], expr_ruvg_pca.transform(expr_ruvg_df.T)[:, 1])
        plt.savefig(f"{project_name}_pca_count_k{k_ruv}.pdf")
        plt.close()

        # Similar plotting for expr.ruvg.pca.vec_count_k

        expr_rmedian = (expr_ruvg_df + 1).median(axis=1).values
        expr_counts_r = pd.DataFrame(index=expr_ruvg_df.index, columns=expr_ruvg_df.columns)

        for row in range(expr_ruvg_df.shape[0]):
            expr_counts_r.iloc[row, :] = expr_ruvg_df.iloc[row, :] + 1 / expr_rmedian[row]

        expr_counts_r_l2 = np.log2(expr_counts_r)
        expr_counts_r_l2[expr_counts_r_l2 == 0] = np.nan

        plt.figure(figsize=(18, 6))
        plt.boxplot(expr_counts_r_l2)
        plt.title("RLE of Unnormalized Data")
        plt.ylabel("RLE log2(count/median count)")
        plt.xlabel("Samples")
        plt.savefig(f"{project_name}_rle_count_k{k_ruv}.pdf")
        plt.close()

# Need to work automating this
k_files = [pd.read_csv(f"nano_24_sec_seq_expression_matrix_log2_norm_k{i}.tsv", sep="\t") for i in range(21)]

expr_counts = np.log2(expr_counts + 1)

for type in range(4):
    expr_ind = expr_counts.index[expr_counts.index.isin(oldham_genes.iloc[:, type])].tolist()
    expr_norm_ind = expr_norm_l2.index[expr_norm_l2.index.isin(oldham_genes.iloc[:, type])].tolist()
    k_inds = [k.index[k.index.isin(oldham_genes.iloc[:, type])].tolist() for k in k_files]

    expr_bicor = bicor(expr_counts.loc[expr_ind])
    expr_bicor = expr_bicor[np.triu_indices(len(expr_bicor), k=1)]
    expr_norm_bicor = bicor(expr_norm_l2.loc[expr_norm_ind])
    expr_norm_bicor = expr_norm_bicor[np.triu_indices(len(expr_norm_bicor), k=1)]

    k_bicors = [bicor(k.loc[k_ind]) for k, k_ind in zip(k_files, k_inds)]
    k_bicors = [k[np.triu_indices(len(k), k=1)] for k in k_bicors]

    type_c = oldham_genes.columns[type]
    type_cor = pd.DataFrame(np.column_stack([expr_bicor, expr_norm_bicor] + k_bicors))
    type_cor.columns = ["Unnormalized", "Global-Scaling"] + [f"RUVg, k={i+1}" for i in range(len(k_bicors))]

    plt.figure(figsize=(15, 6))
    plt.boxplot(type_cor)
    plt.ylabel("Biweight Midcorrelation")
    plt.title(f"Biweight Midcor of Top 50 {type_c} genes")
    plt.savefig(f"bwmcor_top50_{type_c}_genes.pdf")
    plt.close()

