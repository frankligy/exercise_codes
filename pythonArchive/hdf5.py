import numpy
import h5py

'''
Two type of object:
1. group  (file is also a group)
2. dataset

* each object will have a name, which serve as key to look for them
* each dataset will have a attrs dictionary storing associated metadata, using dict(dataset.attrs) to display
* group.visititems(print) to see all its members (groups and datasets)

f
|_____ dataset1
|_____ subgroup1
            |________ dataset2
            |________ dataset3
'''

with h5py.File('/Users/ligk2e/Desktop/test.hdf5','w') as f:
    f.create_dataset('dataset1',(100,),dtype='<i8')
    g = f.create_group('subgroup1')
    g.create_dataset('dataset2',(100,),dtype='<i8')
    f.create_dataset('subgroup1/dataset3',(100,),dtype='<i8')

f = h5py.File('/Users/ligk2e/Desktop/test.hdf5','r')


# if wanna know h5ad, see below:
# https://anndata.readthedocs.io/en/latest/fileformat-prose.html




