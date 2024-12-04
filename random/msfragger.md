# How to run MSFragger?

First please refer to this issue page which detailed how I installed and configured things from scratch (https://github.com/Nesvilab/FragPipe/issues/1914), but
it is more for me not for users.

Now let's talk about how to run MSFragger through the FragPipe. You have raw files and fasta file.

```bash
1. /gpfs/data/yarmarkovichlab/MSFragger/test_thermo_data/20240110_E_OdinLC_IC_PDX_HD_19.raw (raw file)
2. /gpfs/data/yarmarkovichlab/MSFragger/test_db/ensembl_protein_add_decoy.fasta (20k canonical proteome with reverted as decoy)
```

We need to add the fasta to the workflow file (`/gpfs/data/yarmarkovichlab/MSFragger/fragpipe/workflows/Nonspecific-HLA-customized.workflow`), see below

```bash
# Workflow: Nonspecific-HLA

database.db-path=/gpfs/data/yarmarkovichlab/MSFragger/test_db/ensembl_protein_add_decoy.fasta

crystalc.run-crystalc=false
```

We then need to create a `manifest` file to indicate the raw file path (`/gpfs/data/yarmarkovichlab/MSFragger/test_manifest/test_manifest.fp-manifest`), should be self-explainable.

Now we ready to run:

```bash

module load jdk/21.0.1
module load python/cpu/3.9.13
module load mono/5.20.1

/gpfs/data/yarmarkovichlab/MSFragger/fragpipe/bin/fragpipe --headless \
    --workflow /gpfs/data/yarmarkovichlab/MSFragger/fragpipe/workflows/Nonspecific-HLA-customized.workflow \
    --manifest /gpfs/data/yarmarkovichlab/MSFragger/test_manifest/test_manifest.fp-manifest \
    --workdir /gpfs/data/yarmarkovichlab/MSFragger/test_thermo_data \
    --config-tools-folder /gpfs/data/yarmarkovichlab/MSFragger/fragpipe/tools 
    --config-diann /gpfs/data/yarmarkovichlab/MSFragger/fragpipe/tools/diann/1.8.2_beta_8/linux/diann-1.8.1.8 
    --config-python /gpfs/share/apps/python/cpu/3.9.13/bin/python
```

The results should be in here (`/gpfs/data/yarmarkovichlab/MSFragger/test_thermo_data`).






