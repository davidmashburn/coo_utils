'''coo_utils works with numpy as scipy.sparse (2D only) to store
compressed multidimensional arrays both in memory and on disks. This can
DRASTICALLY reduce the space required.

In reality, what is stored is a list-of-lists containing sparse matrix
nodes in the scipy.sparse.coo_matrix format. The storage format can
either be normal or differenced (helpful for same-value blocks data).

Some terms used in the code:

CooHD: list (or deeply nested list(s)) of scipy.sparse.coo_matrix's
       (HD means higher-dimensional)
nnzs: An array with the same shape as the nested lists in the cooHD
      (aka without the 2 coo_matrix dimensions) that gives the length
      of each coo_matrix (nnz) at each node of the tree
RCD: A sparse array represented as a flat row-column-data matrix
    Multi-dimensional shape data is recoverable using the nnzs array

The disk storage format uses the RCD fomat for maximum simplicity:

*_rcd.npy    -  The rcd matrix stored in the .npy format (int32)
*_nnzs.npy   -  The nnzs array store in .npy format (int)
*_shape.txt  -  a simple text file with the full shape of the array in
                comma-separated plain text
'''

import numpy as np
import scipy.sparse

def VerifyCooHDShape(cooHD,shape):
    if len(shape)>2:                      # I'm a tree
        if len(cooHD)!=shape[0]:
            return False                 # I'm a wrong tree
        for i in cooHD:
            if not VerifyCooHDShape(i,shape[1:]):
                return False             # I'm a tree with a bad branch
    elif cooHD.shape!=tuple(shape):       # I'm a wrong leaf
        return False
    return True                          # All checks out!

def GetCooHDnnzs(cooHD,shape):
    if len(shape)>2:
        return np.array([ GetCooHDnnzs(i,shape[1:]) for i in cooHD ])
    else:
        return np.array(cooHD.nnz)

# rcdIndex here is a list purely so that the value can be pass-by object (global updates...)
def ConvertRCDToCooHD(rcdMatrix,nnzs,shape,rcdIndex=None,tolil=False):
    if rcdIndex==None:
        rcdIndex=[0]
    if len(shape)>2:
        return [ ConvertRCDToCooHD(rcdMatrix,nnzs[i],shape[1:],rcdIndex,tolil=tolil) for i in range(shape[0]) ]
    else:
        i = rcdIndex[0]
        row  = rcdMatrix[0,i:i+nnzs]
        col  = rcdMatrix[1,i:i+nnzs]
        data = rcdMatrix[2,i:i+nnzs]
        rcdIndex[0]+=nnzs
        coo = scipy.sparse.coo_matrix((data,(row,col)), shape=shape, dtype=np.int32)
        if tolil:  return coo.tolil()
        else:      return coo

# Find the nnzs using nnzs = GetCooHDnnzs(cooHD,shape)
def ConvertCooHDToRCD(cooHD,nnzs,rcdMatrix=None,rcdIndex=None,fromlil=False):
    if rcdMatrix==None:
        rcdMatrix = np.zeros([3,np.sum(nnzs)],dtype=np.int32) # row,col,data
        rcdIndex=[0]
    if len(nnzs.shape)>0:
        for i in range(nnzs.shape[0]):
            ConvertCooHDToRCD(cooHD[i],nnzs[i],rcdMatrix,rcdIndex,fromlil=fromlil)
    else:
        if fromlil:
            cooHD = cooHD.tocoo()
        i = rcdIndex[0]
        rcdMatrix[0,i:i+nnzs] = cooHD.row
        rcdMatrix[1,i:i+nnzs] = cooHD.col
        rcdMatrix[2,i:i+nnzs] = cooHD.data
        rcdIndex[0]+=nnzs
    return rcdMatrix

def SaveCooHDToRCDFile(cooHD,shape,baseName,fromlil=False):
    if not VerifyCooHDShape(cooHD,shape):
        print 'Shape does not match or cooHD is ill-formed!'
        return
    nnzs = GetCooHDnnzs(cooHD,shape)
    np.save( baseName+'_nnzs.npy', nnzs )
    rcdMatrix = ConvertCooHDToRCD(cooHD,nnzs,fromlil=fromlil)
    np.save( baseName+'_rcd.npy', rcdMatrix )
    shapeStr = ','.join(map(str,shape))
    open(baseName+'_shape.txt','w').write(shapeStr)

def GetShapeFromFile(baseName):
    shapeStr = open( baseName+'_shape.txt' , 'r' ).read()
    shape = tuple(map(int,shapeStr.split(',')))
    return shape

def LoadRCDFile(baseName):
    rcdMatrix = np.load( baseName+'_rcd.npy' )
    nnzs = np.load( baseName+'_nnzs.npy' )
    shape = GetShapeFromFile(baseName)
    return rcdMatrix,nnzs,shape
def LoadRCDFileToCooHD(baseName,tolil=False):
    rcdMatrix,nnzs,shape = LoadRCDFile(baseName)
    if not VerifyCooHDShape(nnzs,shape[:-2]):
        print 'Shape does not match nnzs!'
        return
    return shape,ConvertRCDToCooHD(rcdMatrix,nnzs,shape,tolil=tolil)

def ArrayToCooHD(arr,dtype=np.int32,tolil=True):
    if arr.ndim==2:
        arr=arr.astype(dtype)
        if tolil:
            return scipy.sparse.lil_matrix(arr,dtype=dtype)
        else:
            return scipy.sparse.coo_matrix(arr,dtype=dtype)
    elif arr.ndim>2:
        return [ArrayToCooHD(i) for i in arr]
    else:
        print 'ERROR! array has too few dimensions!'

def CooHDToArray(cooHD,dtype=np.uint16):
    if hasattr(cooHD,'shape'): # seems to be a defining characteristic of all scipy.sparse
        return cooHD.toarray()
    else: # so this must be a list...
        return np.array([CooHDToArray(i) for i in cooHD],dtype=dtype)

# Super-simple compression for integer block matrices
def ArrayToCooDiff(arr,dtype=np.int32,tolil=True):
    '''Watersheds are usually blocks, so by striping the array, we can store it more compactly (diff on axis 0); recover the array with np.cumsum'''
    if arr.ndim==2:
        cooDiff = arr.astype(dtype)
        cooDiff[1:]=np.diff(cooDiff,axis=0)
        if tolil:
            return scipy.sparse.lil_matrix(cooDiff,dtype=dtype)
        else:
            return scipy.sparse.coo_matrix(cooDiff,dtype=dtype)
        return cooDiff
    elif arr.ndim>2:
        return [ArrayToCooDiff(i) for i in arr]
    else:
        print 'ERROR! array has too few dimensions!'
# And uncompression...
def CooDiffToArray(cooDiff,dtype=np.uint16):
    if hasattr(cooDiff,'shape'): # seems to be a defining characteristic of all scipy.sparse
        return np.cumsum(cooDiff.toarray(),axis=0) # that was easy...
    else: # so this must be a list...
        return np.array([CooDiffToArray(i) for i in cooDiff],dtype=dtype)

