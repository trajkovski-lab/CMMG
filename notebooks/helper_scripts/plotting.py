from __future__ import print_function, absolute_import
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from . import Significance
from itertools import combinations
import warnings
import matplotlib
from scipy import stats


try:
    import seaborn as sns

    #rc={'legend.fontsize': 16, 'axes.labelsize': 16, 'xtick.labelsize': 12, 'ytick.labelsize': 12,'axes.titlesize':20,}

    sns.set_context('paper',font_scale=1.5)

except :
    print("Couldn't load seaborn")
    #malkes plot R like
    #plt.style.use('ggplot')

def plotting_Setup(sns_contex='paper',font_scale=1.5, save_dpi=150):
    import matplotlib
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set_context(sns_contex,font_scale=font_scale)
    matplotlib.rcParams['savefig.dpi']=save_dpi
    matplotlib.rcParams['pdf.fonttype']=42

def transform2inch(*tupl, initial_unit='mm'):
    if initial_unit=='cm':
        inch = 2.54
    elif initial_unit=='mm':
        inch = 25.4

    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)


def saveplot(name,figurefolder='../Figures/',formats=['.pdf','.png'],SAVEPLOT=True,print_Latex=False):
    matplotlib.rcParams['pdf.fonttype']=42 # to save as vector format
    if SAVEPLOT:
        if not os.path.exists(figurefolder):
            os.makedirs(figurefolder)
        if len(formats)==0:
            raise Exception('no format specified')

        for format in formats:
            plt.savefig(os.path.join(figurefolder,name+format),bbox_inches='tight')


        if print_Latex:
            print('''\\begin{figure}
        \includegraphics[width=\linewidth]{{0}.pdf}
        \caption{}
        \label{fig:{0}}
    \end{figure}
    '''.format(name))






def format_p_value(p,Stars=False):
    if not Stars:
        return 'P = {:.1g}'.format(p).replace('0.','.')
    else:
        if p<0.001:
            return '***'
        elif p<0.01:
            return '**'
        if p < 0.05:
            return '*'
        else:
            return 'n.s.'

## private method
from matplotlib.font_manager import FontProperties
def label_sig_(x1,x2,y,text='*',ax=None,fontoffset=0,linewidth=1.5,font_size=None):

    if font_size is None:
        font_size='xx-small'

    font=FontProperties(size=font_size)
    if ax is None:
        ax=plt.gca()


    ax.annotate(text, xy=((x1+x2)/2.,y),textcoords='offset points',xytext=(0,fontoffset),
                horizontalalignment='center',verticalalignment='bottom',font_properties=font)
    ax.plot((x1,x2),(y,y),marker=3,markeredgewidth=linewidth, linewidth=linewidth,color='k')

def label_sig(order,group1,group2,y,text,hue_offset=0,ax=None, fontoffset=0, linewidth=1.5,font_size=None):
    ''' wraper for putting a sig label on a plot'''
    group_id=dict(zip(order,range(len(order))))
    label_sig_(group_id[group1]+hue_offset,group_id[group2]+hue_offset,y,text=text,
              ax=ax, fontoffset=fontoffset, linewidth=linewidth,font_size=font_size)



def plot_all_sig_labels(order,K,y0=None,deltay=None,x_offset=0,ax=None,sig_Stars=False,plot_only_significant=True,**kws):
    if ax is None:
        ax=plt.gca()

    Lim=ax.get_ylim()
    if deltay is None:
        deltay=(Lim[1]-Lim[0])/10
    if y0 is None:
        y=Lim[1]+(Lim[1]-Lim[0])/10
    else: y=y0

    if plot_only_significant:
        K_sig=K.loc[K<0.05]
    else:
        K_sig=K

    S=K_sig.apply(format_p_value, Stars=sig_Stars)

    for s in S.iteritems():
        (group1,group2)=s[0]
        label_sig(order,group1,group2,y,text=s[1],hue_offset=x_offset,ax=ax,**kws)
        y+=deltay


def count_sample_for_legend(variable, Grouping_Variable, order):
    '''
    generate legend such as ['E.coli\nn=7', 'Other\nn=7', 'negative\nn=15', 'rotaEAEC\nn=16', 'rotavirus\nn=26']
    '''
    variable_ =pd.Series(variable)
    count_samples= variable_.groupby(Grouping_Variable).count().loc[order]
    return ['{}\nn={}'.format(*label) for label in count_samples.iteritems()]



import string

def annotate_subplot(ax,Letter):
    ax.text(-0.1, 1.1, Letter, transform=ax.transAxes,
        size="medium", weight='bold')


def set_titles(ax,nrows=2,ncols=2):
    Labels=np.array(list(string.ascii_uppercase[:nrows*ncols])).reshape((nrows,ncols))
    ax2=np.atleast_2d(ax.reshape((nrows,ncols)))
    Labels2=np.atleast_2d(Labels)
    for i in range(nrows):
        for j in range(ncols):
            annotate_subplot(ax2[i,j],Labels2[i,j])



def Boxplot(variable,grouping_variable, order=None, ax=None,swarmplot=False,
            grouping_variable2=None, order2=None, Header_Description=None,
            box_params=None,swarm_params=None):

    variable=pd.Series(variable)


    if order is None:
        order=np.unique(grouping_variable)

    sns.set_context('paper',font_scale=1.8)

    if ax is None:
        ax=plt.subplot(111)

    ##default params
    box_params_={}
    #if grouping_variable2 is None:
    #    box_params_= dict(color='darkgray',medianprops={'color':'r'},linewidth=1.5)
    #else:
    #    box_params_ =dict(linewidth=1.5)

    if box_params is not None:
        box_params_.update(box_params)



    sns.boxplot(y= variable, x=grouping_variable,order=order,ax=ax,
                hue=grouping_variable2,hue_order=order2,**box_params_)

    if swarmplot:
        handles, labels = ax.get_legend_handles_labels() #keep boxplot legend
        #default params
        if grouping_variable2 is not None:
            order2_=np.unique(grouping_variable2)
            swarm_params_=dict(palette=dict(zip(order2_,['k']*len(order2_))), size=4, split=True)
        else:
            swarm_params_=dict(color='k')
        if swarm_params is not None:
            swarm_params_.update(swarm_params)


        sns.swarmplot(y= variable, x=grouping_variable,order=order,ax=ax,
                      hue=grouping_variable2,hue_order=order2,**swarm_params_)

        ## only use legend of boxplot
        if grouping_variable2 is not None:

            ax.legend(handles, labels)


    ax.set_xticklabels(count_sample_for_legend(variable, grouping_variable, order))
    ax.set_xlabel('')
    if Header_Description is not None:
        ax.set_ylabel(Header_Description)

    return ax



def N(x): return np.sum(~np.isnan(x))
def Q1(x): return np.percentile(x,25)
def Q3(x): return np.percentile(x,75)

def Sumarize_Table(Data,grouping_variable, order=None ,aggregation_functions=[N,'mean','std',Q1,'median', Q3 ],Pairwise_Sig=None, calculate_diff_for=['median'] ):
    """ Sumarize a datatable grouped by a grouping variable using the functions in aggregation function

        eg. Row 1 mean std Q1 median Q2

        If a table with Pairwise p-values is given it is added to the table
        and the differences for the variables in the list calculate_diff_for are calculated

    """

    if order is None:
        order=np.unique(grouping_variable)

    Summary=Data.groupby(grouping_variable).aggregate(aggregation_functions).T.unstack()

    Summary= Summary[order]

    for measure in calculate_diff_for:
        assert measure in aggregation_functions, "Need {measure} values for {measure} difference calculations".format(measure=measure)

        for (G1,G2) in combinations(order,2):
            comparison='{} <-> {}'.format(G1,G2)
            Summary[(measure+'_diff',comparison)] = Summary[(G2,measure)]-Summary[(G1,measure)]

    if Pairwise_Sig is not None:
    ## corrected P values
        pv=Pairwise_Sig.copy()
        pv.index=[['P_values']*pv.shape[0],
                  pv.index.get_level_values(0).values+' <-> '+pv.index.get_level_values(1).values]
        Summary= pd.concat((Summary,pv.T),axis=1)

    else:
        warnings.warn('Pairwise_Sig is not defined no P_values')
    return Summary

def correct_SumTable_for_multiple_testing(ST,correction_method='Bonferroni'):
    """Adds corected p values to summarize table ST"""

    for comparison in ST.P_values.columns:
        [G1,G2]=comparison.split(' <-> ')
        ST[("P "+correction_method,comparison)]=Significance.correct_pvalues_for_multiple_testing(ST[('P_values',comparison)],correction_method)






def reshape_Data(data, grouping_variable,value_name='value',var_name='variable',data_columns=None):
    """ reshapes data from ID x Gene to a matrix where all the values are in one column"""

    D=data.copy()

    D['ID']=D.index
    D=D.join(grouping_variable)

    return pd.melt(D,id_vars=['ID',grouping_variable.name], value_name=value_name,value_vars=data_columns,var_name=var_name)





## Correlation

def plot_correlation(Data,grouping_variable,order,Correlations,Pvalues,N,g1,g2,corr_function='spearmanr'):
    """# Example
    #V = Viewpoint
    po.plotting_Setup(font_scale=0.75)
    sns.set_palette(['green','darkred'])
    g1,g2='Escherichia','Escherichia phage'

    Correlations,Pvalues,N= Significance.correlation_byGroup(V.Data,V.grouping_variable,V.order)

    plot_correlation(V.Data,V.grouping_variable,V.order,Correlations,Pvalues,N,g1,g2,
    corr_function='spearmanr')
    """
    d_=Data.join(grouping_variable)

    sns.pairplot(d_,hue=grouping_variable.name,hue_order=order,
                x_vars=[g1],y_vars=[g2],size=4)


    statistics=pd.DataFrame()

    statistics[corr_function]= Correlations.loc[:,g1,g2]
    statistics['P values']=Pvalues.loc[:,g1,g2]
    statistics['N']=N

    return statistics


def load_Viewpoint_from_pickle(pickle_filename):
    import pickle

    with open(pickle_filename,'rb') as infile:
        V = pickle.load(infile)
    return V

## Viewpoint

class Viewpoint:
    def __init__(self, Data,grouping_variable,order=None,grouping_variable2=None,order2=None,ref_group=None,
                 Pairwise_Sig="mannwhitneyu", Name='Rawdata',Selection_columns=None,Selection_rows=None,Header_Description=None):
        self.Name=Name

        Data=pd.DataFrame(Data)

        # Selection of subdataset
        if Selection_columns is None: Selection_columns=Data.columns
        if Selection_rows is None: Selection_rows=Data.index
        self.Data=Data.loc[Selection_rows,Selection_columns]

        if type(Selection_rows) == pd.Series:
            Selection_rows=Selection_rows.loc[Selection_rows].index.values # select index values

        self.grouping_variable=grouping_variable.loc[Selection_rows]


        if order is None:
            self.order=np.unique(self.grouping_variable)
        else:

            for group in order:
                if not any(grouping_variable==group):
                    raise NameError("{} not in grouping_variable".format(group))

            self.order=order


        if grouping_variable2 is not None:
            self.grouping_variable2=grouping_variable2.loc[Selection_rows]

            if order2 is None:
                self.order2=np.unique(self.grouping_variable2)
            else:
                self.order2=order2
                for group in order2:
                    if not any(grouping_variable2==group):
                        raise NameError("{} not in grouping_variable2".format(group))
        else:
            self.grouping_variable2 =None
            self.order2=None

        self.ref_group=ref_group

        if Header_Description is None:
            self.Header_Description=pd.Series(self.Data.columns, index=self.Data.columns,name=("Description","Description"))
        else:
            self.Header_Description= pd.Series(Header_Description.copy())
            self.Header_Description.name=name=("Description","Description")

            assert self.Data.columns.isin(self.Header_Description.index).all(), "Not all headers have a description!"

        if type(Pairwise_Sig)==str:
            self.Pairwise_Sig,self.Statistic_Values=Significance.Paired_StatsTest(self.Data,
                              self.grouping_variable,
                              test=Pairwise_Sig,
                              ref_group=self.ref_group,
                              return_statistic=True)
        else:
            assert type(Pairwise_Sig)==pd.core.frame.DataFrame, "Pairwise_Sig should be a string for a statistical test or the output of Significance.Paired_StatsTest"
            self.Pairwise_Sig=Pairwise_Sig


    def to_pickle(self, pickle_filename):
        import pickle

        out_folder= os.path.dirname(pickle_filename)

        if not os.path.exists(out_folder): os.makedirs(out_folder)

        with open(pickle_filename,'wb') as outf:
            pickle.dump(self,outf)


    def Boxplot(self,variable_to_plot, ax=None,swarmplot=True,Header_Description=None,
            box_params=None,swarm_params=None,plot_sig_labels=True,sig_labels_params=None):

        if ax is None:
            ax=plt.subplot(111)


        Boxplot(self.Data[variable_to_plot],order=self.order,
                grouping_variable=self.grouping_variable,
                ax=ax,swarmplot=swarmplot,Header_Description=self.Header_Description[variable_to_plot],
                box_params=box_params,swarm_params=swarm_params,
                grouping_variable2=self.grouping_variable2, order2=self.order2)


        if plot_sig_labels:
            if sig_labels_params is None: sig_labels_params={}


            self.plot_all_sig_labels(variable_to_plot,ax=ax,**sig_labels_params)
        return ax


    def plot_all_sig_labels(self,variable_to_plot,y0=None,deltay=None,x_offset=0,ax=None,sig_Stars=False,**kws):
        plot_all_sig_labels(self.order,self.Pairwise_Sig[variable_to_plot], y0=y0,deltay=deltay,x_offset=x_offset,ax=ax,sig_Stars=sig_Stars,**kws)


    def Sumarize_Table(self, aggregation_functions=[N,'mean','std',Q1,'median', Q3 ],correction=None ,calculate_diff_for=['median']):
        """correction one of "Bonferroni", "Bonferroni-Holm", "Benjamini-Hochberg" """
        ST=Sumarize_Table(self.Data,self.grouping_variable, self.order, aggregation_functions=aggregation_functions,
                             Pairwise_Sig=self.Pairwise_Sig, calculate_diff_for=calculate_diff_for)


        if correction is not None:
            correct_SumTable_for_multiple_testing(ST,correction_method=correction)

        ST[('Description','Description')]=self.Header_Description

        return ST



    def Boxplot_all_data(self,ax=None,plot_params=None,
                      sns_plot_function=sns.boxplot,x_tick_param={"rotation":45, "horizontalalignment":"right"},
                        y_label='value',legend_title=None,data_columns=None):
        """
        Box plot of all data grouped according to grouping variable

        Parameters
        ----------
        data_columns: list-like
            list of parameter to be plotted

        ax:
            axis to plot
        sns_plot_function: function
            seaborn plotting function default boxplot
        plot_params: dict



        Optional parameters
        -------------------

        x_tick_param



        """
        if legend_title is None:
            var_name=''
        else:
            var_name=legend_title

        if ax is None:
            ax= plt.subplot(111)
        if plot_params is None: plot_params= {}


        data=reshape_Data(self.Data,self.grouping_variable,value_name=y_label,var_name=var_name,data_columns=data_columns)

        data[var_name]= data[var_name].map(self.Header_Description)

        hue_order= self.Header_Description.loc[data_columns]

        sns_plot_function(hue=var_name,x= self.grouping_variable.name,
                          order= self.order, y=y_label,data=data,hue_order= hue_order,ax=ax, **plot_params)

        _= ax.set_xticklabels(ax.get_xticklabels(),**x_tick_param)


        return ax






        return Boxplot_all_data(self.Data,self.grouping_variable,order=self.order,ax=ax,plot_params=plot_params,
                  sns_plot_function=sns_plot_function,data_columns=data_columns, Header_Description= self.Header_Description, **param)
## 2x2 plot

#f, ax=plt.subplots(2,2,sharex=True,figsize=(10,7))

#set_titles(ax)

#metaplot('WAZ',ax=ax[0,0])
#plot_all_sig_labels(order,K.loc['WAZ'],y0=4,deltay=1.1,ax=ax[0,0])
#metaplot('HAZ',ax=ax[0,1])
#plot_all_sig_labels(order,K.loc['HAZ'],y0=4,deltay=1.1,ax=ax[0,1])
#metaplot('MUAC',ax=ax[1,0])
#plot_all_sig_labels(order,K.loc['MUAC'],y0=19,deltay=2,ax=ax[1,0])
#metaplot('BMI',ax=ax[1,1])
#plot_all_sig_labels(order,K.loc['BMI'],deltay=1.8,y0=22,ax=ax[1,1])
#plt.tight_layout()

#saveplot('MetadataPlot_'+'BiometricVariables',figurefolder=os.path.join(out_folder,"FinalPlots"),SAVEPLOT=SAVEPLOT)
#plt.show()


def boxplot_with_connections(V,varaible_2_plot,SampleID,ax=None):
    """plots data from viewpoint V grouped by grouping_variable 2 and grouping_variable and
    connects all samples with the same SampleID, to better see the evolution within the different groups (grouping_variable)

    assumes there are only one data point with each sampleID in eaach group of grouping_variable
    """
    assert V.grouping_variable2 is not None

    if ax is None:
        ax= plt.subplot(111)


    Boxplot(V.Data[varaible_2_plot],V.grouping_variable2,order=V.order2,
               grouping_variable2=V.grouping_variable,order2=V.order,ax=ax)

    G=V.Data[ varaible_2_plot].groupby(V.grouping_variable2)

    for i,group_name in enumerate(V.order2):
        for id,d in G.get_group(group_name).groupby(SampleID):
            ax.plot(np.array([0,1])*0.4-0.2+i,d.values,'o-k')


def Boxplot2groupingvaraibles(Data,variable2plot,grouping_variable,order,grouping_variable2,order2,y0=None):
    """Experimental: polts data in seperate subplots for the grouping groupingvariable and statistics for groupingvariable2 """
    N=len(order)

    f,axe= plt.subplots(1,N,sharey=True)

    Viewpoints=[Viewpoint(Data[variable2plot],
                             grouping_variable=grouping_variable2,
                             order=order2,
                             Selection_rows=(grouping_variable==order[i])) for i in range(N)]

    for i in range(N):
        Viewpoints[i].Boxplot(variable2plot,ax=axe[i])
        Viewpoints[i].plot_all_sig_labels(variable2plot,ax=axe[i],y0=y0)

        axe[i].set_title(order[i])
    axe[1].set_ylabel('')
    f.suptitle(grouping_variable.name)

def Matrix_Plot(V,dim,y0=None,fig_name=None,y_label=None,vars2plot=None,plot_all_sig_labels=True,Boxplot_param=None,
    figsize=None,sharey=True):
    """plot several boxplots for vars2plot in a matrix defined by dim
    y0 is used for significanc elabels or as max value

    V Viewpoint
    """


    if vars2plot is None:
        vars2plot=V.Data.columns
    if sharey:
        if y0 is None:
            y0 = V.Data.max().max()

        deltay=0.2*y0
    else:
        y0=None
        deltay=None

    N= len(vars2plot)

    nrows,ncols=dim
    assert N ==nrows*ncols, 'Missmatch between dimensions and len(vars2plot): ({} x {}) != {N}'.format(*dim,N=N)

    # by default no swarmplot is displayed when two grouping variables are used
    if Boxplot_param is None:
        Boxplot_param ={}
    if V.grouping_variable2 is not None:
        Boxplot_param_,Boxplot_param= Boxplot_param,dict(swarmplot=False)
        Boxplot_param.update(Boxplot_param_)


    f,axe= plt.subplots(nrows,ncols,sharey=sharey,sharex=True,figsize=figsize)

    for i in range(nrows):
        for j in range(ncols):
            var=vars2plot[i*ncols+j]
            ax=axe[i,j]
            V.Boxplot(var,axe[i,j],plot_sig_labels=False, **Boxplot_param)
            if plot_all_sig_labels:
                V.plot_all_sig_labels(var,ax=axe[i,j],y0=y0,deltay=deltay)

            ax.set_title(ax.get_ylabel(),fontdict={'fontsize':'xx-small'})
            ax.set_ylabel('')
            if ax.legend_ is not None:
                ax.legend_.set_visible(False)

    if sharey:

        ax.set_ylim(bottom=V.Data.min().min()- y0*0.1)

        if not plot_all_sig_labels:
            ax.set_ylim(top=y0)
    if y_label is not None:
        axe[0,0].set_ylabel(y_label)
        axe[1,0].set_ylabel(y_label)

    if fig_name is not None:

        f.suptitle(fig_name)

    #legend
    ax= axe[0,0]
    if ax.legend_ is not None:
        handles,labels=ax.get_legend_handles_labels()

        ax.legend(handles[:3],labels[:3],
                        bbox_to_anchor=(0.07,1),loc='upper left',bbox_transform=f.transFigure,ncol=3,fontsize='xx-small')



def multiplot(V,n_rows,n_cols,subset=None,sharey=True,figsize=None,sig_labels_params=None,**kws):
    """
    Plots multiple variables from the Viewpont as a multi pannel figure, define subset, n_col and n_rows
    """

    N=  n_cols*n_rows

    if subset is None:
        subset= V.Data.columns[:N]

    f,axe = plt.subplots(n_rows,n_cols,sharey=sharey,figsize=figsize)
    if sig_labels_params is None:
        sig_labels_params={}
        if sharey:
            max_value= V.Data[subset].max().max()
            sig_labels_params=dict(y0=max_value,deltay=max_value/10.)

    axe= np.ravel(axe)

    for i in range(N):
        V.Boxplot(subset[i],ax= axe[i],sig_labels_params=sig_labels_params,**kws)
    import string

    for i,ax in enumerate(axe):
        annotate_subplot(ax,string.ascii_uppercase[i])
    return axe




def plot_trace_(data, color=None, alpha=0.1, what='MeanSem',label=None,ax= None, **kws):

    if ax is None:
        ax = plt.subplot(111)

    if what.startswith('Mean'):
        data.mean().plot(marker='o',color=color,label=label,ax=ax,**kws)
        color=ax.lines[-1].get_color()

        if what=='MeanStd':


            ax.fill_between(data.columns,data.mean()+data.std(ddof=1),data.mean()-data.std(ddof=1),
                     alpha=alpha,color=color)

        elif what=='MeanSem':

            ax.fill_between(data.columns,data.mean()+data.sem(),data.mean()-data.sem(),
                     alpha=alpha,color=color)

        elif what=='MeanCI':

            ci=stats.t.ppf((1.95) / 2, data.shape[0] - 1)*data.sem()

            ax.fill_between(data.columns,data.mean()+ci,data.mean()-ci,
                     alpha=alpha,color=color)
        elif not what=='Mean':
            raise Exception("'what' must be one of 'MeanStd','Mean','MeanCI', 'MeanSem','MedianIQR'")


    elif what=='MedianIQR':
        data.median().plot(marker='o',color=color,label=label,ax=ax,**kws)
        color=ax.lines[-1].get_color()
        ax.fill_between(data.columns,data.quantile(q=0.25),data.quantile(q=0.75),
                     alpha=alpha,color=color)
    else:
        raise Exception("'what' must be one of 'MeanStd','MeanCI','Mean', 'MeanSem','MedianIQR'")

def plot_traces(data,grouping_variable ,order=None,colors=None, what='MedianIQR',ax=None,**kws):
    if ax is None:
        ax = plt.subplot(111)

    if order is None:
        order= grouping_variable.unique()

    if colors is None:
        colors = dict(zip(order, sns.color_palette('deep',len(order))))

    G= data.groupby(grouping_variable)
    for g in order:
        plot_trace_(G.get_group(g),label=g,color= colors[g],ax=ax,what=what, **kws)

    ax.legend()

    return ax



from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms

from scipy.stats import chi2


def confidence_ellipse(x, y, ax,ci=0.95, color='red',facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    **kwargs
        Forwarded to `~matplotlib.patches.Ellipse`

    Returns
    -------
    matplotlib.patches.Ellipse
    """

    if ax is None:
        ax=plt.gca()

    if len(x) < 4:
        raise InputError("need more than 3 data points")

    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor,edgecolor=color ,**kwargs)


    s= chi2.ppf(ci,2)

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0] * s)
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1] * s)
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)
