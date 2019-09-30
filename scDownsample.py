class Downsample:
    def __init__(self, data, tsneData, index, dim=2, filterCell=True, leftThreshold=0.7, rightThreshold=0.7,
                 calculateCorrelation=True):

        self.data = data
        self.tsneData = tsneData
        self.index = index
        self.dim = dim
        self.filterCell = filterCell
        self.leftThreshold = leftThreshold
        self.rightThreshold = rightThreshold
        self.calculateCorrelation = calculateCorrelation

        '''
        data: anndata
        tsnedata: np.array of tsne
        index: cell name
        '''
        
    def Mix(self):    
        
        def downsample():

            if self.dim == 2:
                df_tsne_cor=pd.DataFrame(self.tsneData, columns=['TSNE-1','TSNE-2'],index=self.index)
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
                ######print(sampled_tsne_1)
                df_tsne_cor_sampled = df_tsne_cor[df_tsne_cor['TSNE-1'].isin(sampled_tsne_1)]
            
            if self.dim == 1:
                df_tsne_cor=pd.DataFrame(self.tsneData, columns=['TSNE-1'],index=self.index)
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
                df_tsne_cor_sampled = df_tsne_cor[df_tsne_cor['TSNE-1'].isin(sampled_tsne_1)]

            #######print(df_tsne_cor_sampled)
            print("cell number before downsample:",len(df_tsne_cor))    
            print("cell number after downsample:",len(df_tsne_cor_sampled))
            print("filter rate:", np.around(1-len(df_tsne_cor_sampled)/len(df_tsne_cor),2))
            self.ds_data = self.data[df_tsne_cor_sampled.index,:]
            #######print(self.ds_data )
            #######print(self.data )
            sc.pp.neighbors(self.ds_data, n_neighbors=20, n_pcs=40)
            sc.tl.umap(self.ds_data)
            try:
                sc.tl.louvain(self.ds_data)
            except:
                sc.tl.leiden(self.ds_data)
            self.ds_data.obs['raw_louvain']=self.data[self.ds_data.obs.index].obs.iloc[:,0]
            #print(self.ds_data)
                       
            return self.data, self.ds_data
        downsample()
        
        if self.filterCell:
            def filterCell():
                # Purify cell population          

                ad=self.ds_data.obs

                a={}
                for i in range(len(set(ad.iloc[:,0]))):
                    a[i]=Counter(ad[ad.iloc[:,0]=='%s' %i].loc[:,'raw_louvain'])

                b={}
                for i in range(len(set(self.data.obs.iloc[:,0]))):
                    b[i]=Counter(ad[ad['raw_louvain']=='%s' %i].iloc[:,0])

                cell_index=[]
                for i in range(len(a)):
                    for j in a[i]:
                        p1 = a[i][j] / sum(a[i].values())

                        if p1 > self.rightThreshold:
                            cell_index.append(ad[(ad['raw_louvain']==j)&(ad.iloc[:,0]=='%s' %i)].index)

                        else:
                            p2 = b[int(j)][i] / sum(b[int(j)].values())
                            if p2 > self.leftThreshold:
                                cell_index.append(ad[(ad['raw_louvain']==j)&(ad.iloc[:,0]=='%s' %i)].index)

                h=[]
                for i in range(len(cell_index)):
                    h += cell_index[i].values.tolist()

                self.ds_data=self.ds_data[h]
                return self.data, self.ds_data
            filterCell()
            #print(self.ds_data,"filter")
        
        if self.calculateCorrelation:
            def clusterCorrelation():
                ds_corr=[]
                for i in set(self.ds_data.obs.iloc[:,0]):
                    ds_filter=[]
                    test={}
                    raw_ds=[]
                    ####print(self.ds_data,"a")
                    ds_filter += [self.ds_data[self.ds_data.obs.iloc[:,0]=='%s' %i].X.mean(axis=0)]
                    test[i]=set(self.ds_data[self.ds_data.obs.iloc[:,0]=='%s' %i].obs.loc[:,'raw_louvain'])
                    #print(ds_filter,"b")
                    for j in test[i]:
                        #print(self.data.obs.iloc[:,0]==j,"c")
                        raw_ds = self.data[self.data.obs.iloc[:,0]==j].X.mean(axis=0)
                        #print(raw_ds.shape,"d")
                        #print(np.array(ds_filter).shape)
                        test=pd.DataFrame({"ds_filter": np.array(ds_filter).reshape(np.array(ds_filter).shape[1]), 
                            "raw_ds": np.array(raw_ds)})
                        #print(ds_filter,"e")
                        #calculate the correlation
                        ds_corr.append(np.triu(test.corr(method="spearman"))[0][1])
                        #print(ds_corr,"f")
                self.ds_corr=ds_corr
                return self.ds_corr
            clusterCorrelation()
        
        return self.data,self.ds_data, self.ds_corr