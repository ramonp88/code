import pandas as pd
import numpy as np
from sklearn.cross_validation import KFold

#filename = '/ramon/MEGAsync/prodgraph/code/RandomWalk/Debug/small.csv'
filename = '/ramon/MEGAsync/prodgraph/movielens/ml-latest-small/ratings.csv'

############ CLEANING DATA SET

def clean_data(ds, key, thres):
    v = ds.groupby(key).size()
    res = [i for i in v.keys() if v[i] >= thres]
    idx = ds[key].isin(res)
    da = ds[idx]
    return da

def reduce_dataset(ds, key, nkeep):
    res = ds.groupby(key).size().sort_values(ascending=False)
    idxs = ds[key].isin(res.keys()[:nkeep])
    da = ds[idxs]
    return da

def do_reduction(filename):
    ds = pd.read_csv(filename)
    n = ds['movieId'].unique().size
    reduced = reduce_dataset(ds, 'movieId', n/2)
    reduced = clean_data(reduced, 'movieId', 5)
    return reduced
    
'''
reduced = do_reduction(filename)
reduced.to_csv('reduced.csv', index=False, header=True)
'''


############################ ALGORITHM
    
n_folds = 5
alpha = 1e-6 #alpha distribuicao dirichlet
p_rw = 0.8 #probabilidade de sucesso da distribuicao geometrica para tamanho random walk

ds = pd.read_csv(filename)

#aux = pd.read_csv('/home/ramon/MEGAsync/prodgraph/code/RandomWalk/Debug/small.csv')

kf = KFold(ds.shape[0], n_folds, shuffle=True)
for train, test in kf:
    df = ds.iloc[train,:]
    aux = ds.iloc[train,:]

    ## substituindo notas por escala inteira de 1 - max
    c = 0
    for i in sorted(aux['rating'].unique()):
        df.loc[aux['rating'] == i,'rating'] = c
        c += 1 

    ## substituindo movieId por 1 - max
    c = 0
    for i in sorted(aux['movieId'].unique()):
        df.loc[aux['movieId'] == i,'movieId'] = c
        c += 1

    nratings = len(df['rating'].unique())
    nmovies = len(df['movieId'].unique())

    ## calculando vetor g
    g = [0]*nratings 
    for i,v in df['rating'].value_counts(normalize=True).iteritems():
        g[int(i)] = v


    ## contando quantidade de notas por faixa
    grouped = df.groupby('movieId')['rating']
    d = dict()
    for name,group in grouped:
        d[name] = [0]*nratings
        
        for i,v in group.value_counts(normalize=True).iteritems():
            d[name][int(i)] = v


    ## criando lista de vizinhanca
    import time
    start_time = time.time()
    m = [[] for i in range(nmovies)]
    for i in range(nmovies):
        i_idx = df.loc[df['movieId'] == i, 'userId']
        for j in range(i+1,nmovies):
            j_idx = df.loc[df['movieId'] == j, 'userId']
            if(len(set(i_idx).intersection(set(j_idx))) > 0):
                m[i].append(j)
                m[j].append(i)

    elapsed_time = time.time()-start_time
    print(elapsed_time)
    break

    u = 1
    u_prod = df.loc[df['userId'] == u]

    g *= alpha

    n_max = 100

    walks = []
    for i in range(n_max):
        for i in u_prod:
            rw = [i]
            len_rw =  np.random.geometric(p_rw, 1)[0]
            while(len(rw) < len_rw):   
                p = rw[-1]
                wmax = -1
                for j in m[p]:
                    theta = np.random.dirichlet(alpha + d[j], 1)
                    y = np.flatnonzero(np.random.multinomial(1, theta, size=1))
                    if(y > wmax):
                        wmax = y
                        p = j
                rw.append(p)
            walks.append(rw)
'''
grouped = ds.groupby('movieId')['rating']
#grouped.size()
for name,group in grouped:
    print(name)
    print(group)
    group.size()
    break
'''