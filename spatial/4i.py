import scanpy as sc
import squidpy as sq

# load
adata = sq.datasets.four_i()


# plot
sc.pl.spatial(adata, color="cluster", spot_size=1)