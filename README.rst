coo_utils works with numpy as scipy.sparse (2D only) to store
compressed multidimensional arrays both in memory and on disks. This can
DRASTICALLY reduce the space required.

In reality, what is stored is a list-of-lists containing sparse matrix
nodes in the scipy.sparse.coo_matrix format. The storage format can
either be normal or differenced (helpful for same-value blocks data).

Some terms used in the code:

CooHD: list (or deeply nested list(s)) of scipy.sparse.coo_matrix's (HD means higher-dimensional)

nnzs: An array with the same shape as the nested lists in the cooHD (aka without the 2 coo_matrix dimensions) that gives the length of each coo_matrix (nnz) at each node of the tree

RCD: A sparse array represented as a flat row-column-data matrix. Multi-dimensional shape data is recoverable using the nnzs array

The disk storage format uses the RCD fomat for maximum simplicity:

\*_rcd.npy    -  The rcd matrix stored in the .npy format (int32)

\*_nnzs.npy   -  The nnzs array store in .npy format (int)

\*_shape.txt  -  a simple text file with the full shape of the array in comma-separated plain text
