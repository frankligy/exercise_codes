## How to set up UCSC cell browser and cellxgene

This tutorial applies to both local machine and server (cchmc)

- See [cellxgene](#cellxgene).
- See [UCSC cell browser](#ucsc-cell-browser).



## cellxgene

First enter to the environment, install the package, may also need to install `scanpy` and `anndata`.

```
pip install cellxgene   # if on local machine
./bin/python3.6 -m pip install cellxgene  # on cchmc server
```

Then start to run,

```shell
cellxgene launch /path/to/file.h5ad --open # on local machine
cellxgene launch /path/to/file.h5ad --open --host server_id  # cchmc
```

`server_id` can be accessed after launch an interactive `bsub`, then type `hostname`, paste the output to the `--host` argument.

Open the browser to see the result. Remember, first use `adata.raw.to_data().write()` to make sure all the genes are included in `adata`.




## UCSC cell browser

First enter to the environment, install the package, may also need to install `scanpy` and `anndata`.

```shell
pip install cellbrowser   # if on local machine
./bin/python3.6 -m pip install cellbrowser  # on cchmc server
```

Then start to run,

```shell
# make sure the adata.obs only contains one column named celltype
cbImportScanpy -i /path/to/file.h5ad -o /a/folder/hold/scanpy/data
cd /a/folder/hold/scanpy/data
cbBuild -o /a/folder/hold/web/files -p 8888
# then open the page, either localhost:8888 or server_id:8888
```

`server_id` is the same means to get as [cellxgene](#cellxgene).