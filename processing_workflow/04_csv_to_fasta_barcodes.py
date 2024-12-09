import pandas as pd
import numpy as np

path = "/Users/Patrick/Downloads/"  # don't forget final "/"
file = "cas_guides.csv"
wpath = path + file
barcodes = pd.read_csv(wpath, sep=",", header=None)

b_mat = np.empty((barcodes.shape[0] * 2, 1), dtype=object)

even = np.arange(2, barcodes.shape[0] * 2 + 1, 2)
odd = np.arange(1, barcodes.shape[0] * 2 + 1, 2)

names = barcodes.iloc[:, 0].apply(lambda x: f">{x}").tolist()
b_mat[odd - 1] = names  # Adjust for 0-based indexing

seq = barcodes.iloc[:, 1].astype(str).tolist()
b_mat[even - 1] = seq  # Adjust for 0-based indexing

output_path = path + "singe_cell_barcodes.fasta"
np.savetxt(output_path, b_mat, fmt='%s', delimiter=',', header='', comments='')


