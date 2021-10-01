def make_sure_adata_writable(adata,delete=False):
    '''
    maks sure the adata is able to write to disk, since h5 file is stricted typed, so on mixed dtype is allowd.
    this function basically is to detect the column of obs/var that are of mixed types, and delete them.
    :param adata: Anndata
    :param delete: boolean, False will just print out what columns are mixed type, True will automatically delete those columns
    :return: Anndata
    Examples::
        from sctriangulate.preprocessing import make_sure_adata_writable
        make_sure_adata_writable(adata,delete=True)
    '''
    # check index, can not have name
    var_names = adata.var_names
    obs_names = adata.obs_names
    var_names.name = None
    obs_names.name = None
    adata.var_names = var_names
    adata.obs_names = obs_names

    # make sure column is pure type, basically, if mixed tyep, delete the column, and print out the delete columns
    # go to: https://github.com/theislab/scanpy/issues/1866

    var = adata.var
    obs = adata.obs
    
    for col in var.columns:
        if var[col].dtypes == 'O':
            all_type = np.array([type(item) for item in var[col]])
            first = all_type[0]
            if (first==all_type).all() and first == str:   # object, but every item is str
                continue
            else:   # mixed type
                print('column {} in var will be deleted, because mixed types'.format(col))
                if delete:
                    adata.var.drop(columns=[col],inplace=True)
    
    for col in obs.columns:
        if obs[col].dtypes == 'O':
            all_type = np.array([type(item) for item in obs[col]])
            first = all_type[0]
            if (first==all_type).all() and first == str:   # object, but every item is str
                continue
            else:   # mixed type
                print('column {} in obs will be deleted, because mixed types'.format(col))
                if delete:
                    adata.obs.drop(columns=[col],inplace=True)

    return adata


















