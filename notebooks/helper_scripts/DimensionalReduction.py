from sklearn.decomposition import PCA, SparsePCA
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import warnings

import seaborn as sns

import altair as alt


def label_points_(data,ax,max_labels=10):

    assert data.shape[1]==2, 'Expect data n x 2'

    sample_names= data.index

    N_samples = len(sample_names)
    if N_samples < max_labels:

        for i in range(N_samples):
            ax.annotate(s= sample_names[i] ,
                        xy=(data.iloc[i,0],data.iloc[i,1]),
                         xytext=(0,10),textcoords='offset points')


class DimensionalReduction:
    def __init__(self, data, decomposition= PCA , transformation= None, n_components=None, **kargs):



        if n_components is None:
            n_components = data.shape[0]


        if data.shape[0] > data.shape[1]:
            print("you don't need to reduce dimensionality or your dataset is transposed.")

        self.decomposition=  decomposition(n_components=n_components, **kargs)


        self.rawdata = data

        #self.variable_names = data.columns

        #self.sample_names = data.index

        if transformation is None:
            self.data_ = self.rawdata

        else:

            self.data_ = data.applymap(transformation)

        Xt = self.decomposition.fit_transform(self.data_)

        self.transformed_data= pd.DataFrame(Xt[:,:(n_components+1)] , index= data.index, columns= np.arange(n_components)+1 )


        name_components = ["components_"]

        for name in name_components:
            if hasattr(self.decomposition,name):
                self.components = pd.DataFrame(getattr(self.decomposition,name), index= np.arange(n_components)+1, columns= data.columns)

        if not hasattr(self,"components"):
            warnings.warn("Couldn't define components, wil not be able to plot loadings")


    def set_axes_labels_(self,ax,components):
        if hasattr(self.decomposition,'explained_variance_ratio_') :

            ax.set_xlabel("PC {} [{:.1f} %]".format(components[0],self.decomposition.explained_variance_ratio_[components[0]-1]*100))
            ax.set_ylabel("PC {} [{:.1f} %]".format(components[1],self.decomposition.explained_variance_ratio_[components[1]-1]*100))

        else:

            ax.set_xlabel("Component {} ".format(components[0]))
            ax.set_ylabel("Component {} ".format(components[1]))

    def plot_explained_variance_ratio(self,**kwargs):

        explained_variance_ratio = self.decomposition.explained_variance_ratio_

        return plt.bar(np.arange(len(explained_variance_ratio)),
                explained_variance_ratio,
                **kwargs
                )


    def plot_Components_2D(self,components=(1,2),ax= None,**scatter_args):

        components = list(components)
        assert len(components)==2, "expect two components"

        if ax is None:
            ax= plt.subplot(111)

        sns.scatterplot(components[0],components[1],data=self.transformed_data,ax=ax,**scatter_args)

        ax.axis('equal')
        self.set_axes_labels_(ax,components)
        label_points_(self.transformed_data[components],ax)

        return ax



    def plot_Loadings_2D(self,components=(1,2),ax= None,**scatter_args):


        if ax is None:
            ax= plt.subplot(111)



        x,y = self.components.loc[components[0]], self.components.loc[components[1]]


        sns.scatterplot(x,y,ax=ax,**scatter_args)
        ax.axis('equal')
        self.set_axes_labels_(ax,components)




        return ax



import altair as alt

import warnings


import altair as alt
def altair_plot2D(data,variables=None, **kws):

    if variables is None:
        variables = data.columns[:2]


    plot_data = data[variables].copy()

    if plot_data.shape[0] > 5000:

        warnings.warn("The number of rows in your dataset is greater than the maximum allowed (5000). I subsample. For information on how to plot larger datasets in Altair, see the documentation")


        plot_data = plot_data.loc[plot_data.abs().max(1).sort_values(ascending=False).index[:4998]]


    plot_data=plot_data.reset_index()
    plot_data.columns= plot_data.columns.astype(str)

    chart= alt.Chart(plot_data).mark_point().encode(
        x=str(variables[0]),
        y=str(variables[1]),
        tooltip= plot_data.columns[0]
    ).interactive()

    return chart
