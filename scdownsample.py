import numpy as np


class DownSampler(object):
    """
    downsample scRNA-seq datasets
    use cellID as key for downsampling
    """

    def __init__(self):
        pass
    
    def DE(x):
        import scanpy as sc
        import anndata as ann

        adata = ann.Anndata(x)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        adata = adata[:, adata.var['highly_variable']]
        return adata

    def downsample(
        self,
        cellID :np.ndarray = None,
        cell_types :np.ndarray = None,
        TPM :np.ndarray = None,
        method = 'cluster',
        cellNumber :aim downsample cell number = None,
        *args,
        **kwargs,
    ) -> np.ndarray:
        
        """
        user cellID as key for downsampling
        input original cellID and other accompanying info and returns selected(downsampled) cellID
        """
        
        original_cellID = cellID


        if method == 'cluster':
            # algorithms here
            
            dsrate = cellNumber / len(cellID)
            sampled_ds_1 = []
            for i in set(celltype):
                indexes = np.where(celltype == i)
                sampling_number = int(sum(celltype == i) * dsrate)
                #print(sampling_number)
                ds_1 = np.random.choice(indexes[0], sampling_number, replace=False)
                sampled_ds_1 = np.append(sampled_ds_1, ds_1, axis=0)
            #print(len(sampled_ds_1))
            #print(max(sampled_ds_1).max())
            selected_cellID = celltype[sampled_ds_1]
            pass

        if method == '2d-tSNE': 

            from openTSNE import TSNE
            from openTSNE.callbacks import ErrorLogger

            #######calculate 2d-tsne
            tsne = TSNE(perplexity=30,
            metric="cosine",
            callbacks=ErrorLogger(),
            initialization='pca',
            random_state=42,
            early_exaggeration_iter=50,
            n_iter=250,
            n_jobs=-1)

            TPM = self.DE(TPM)
            embedding_train = tsne.fit(TPM)
            df_tsne_cor = pd.DataFrame(embedding_train, columns=['TSNE-1','TSNE-2'], index=original_cellID)
            raw_cor = df_tsne_cor.to_numpy()
            floor_cor = np.floor(raw_cor)

            grids = np.array([str(floor_cor[i][0]) + ' * ' + str(floor_cor[i][1]) for i in range(0, floor_cor.shape[0])])
            df_tsne_cor['grid'] = grids
            counts = Counter(grids)

            sampled_tsne_1 = np.array([]).reshape(-1, 1)
            number_in_grid = []
            for grid in set(grids):
                # print('%sth itetration' %i)
                indexes = np.where(grids == grid)
                grid_cor = raw_cor[indexes]

                # sampling number please see here.
                sampling_number = np.floor(np.log(counts[grid] + 1)).astype(int)
                _tsne_1 = np.random.choice(grid_cor[:, 0], sampling_number, replace=False)
                number_in_grid.append(_tsne_1.shape[0])
                sampled_tsne_1 = np.append(sampled_tsne_1, _tsne_1.reshape(-1,1), axis=0)

            sampled_tsne_1=sampled_tsne_1.reshape(len(sampled_tsne_1))
            selected_cellID = df_tsne_cor[df_tsne_cor['TSNE-1'].isin(sampled_tsne_1)].index
            pass

        elif method == '1d-tSNE':

            from openTSNE import TSNE
            from openTSNE.callbacks import ErrorLogger

            # algorithms here
            #######calculate 1d-tsne
            tsne = TSNE(perplexity=30,
            metric="cosine",
            callbacks=ErrorLogger(),
            initialization='pca',
            random_state=42,
            early_exaggeration_iter=50,
            n_iter=250,
            n_components=1,
            n_jobs=-1)
    
            TPM = self.DE(TPM)
            embedding_train = tsne.fit(TPM)

            df_tsne_cor=pd.DataFrame(embedding_train, columns=['TSNE-1'],index=original_cellID)
            ######print("cell number before downsample:",len(df_tsne_cor))
            raw_cor = df_tsne_cor.to_numpy()
            floor_cor =  np.floor(raw_cor)

            intervels = np.array([str(floor_cor[i]) for i in range(0, floor_cor.shape[0])])
            counts = Counter(intervels)
            df_tsne_cor['interval'] = floor_cor

            sampled_tsne_1 = np.array([]).reshape(-1,1)
            number_in_grid = []
            for i in set(intervels):
                # print('%sth itetration' %i)
                indexes = np.where(intervels == i)
                interval_cor = raw_cor[indexes]

                ##sampling number please see here.
                sampling_number = np.floor(np.sqrt(counts[i])).astype(int)
                _tsne_1 = np.random.choice(interval_cor[:, 0], sampling_number, replace=False)
                number_in_grid.append(_tsne_1.shape[0])
                sampled_tsne_1 = np.append(sampled_tsne_1, _tsne_1.reshape(-1,1), axis=0)

            sampled_tsne_1=sampled_tsne_1.reshape(len(sampled_tsne_1))
            selected_cellID = df_tsne_cor[df_tsne_cor['TSNE-1'].isin(sampled_tsne_1)].index
            pass

        return selected_cellID