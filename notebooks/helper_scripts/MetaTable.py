import pandas as pd
from pandas import DataFrame
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns


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

def correct_pvalues_for_multiple_testing(pvalues, correction_type = "Benjamini-Hochberg"):
    """
    correction_type: one of "Bonferroni", "Bonferroni-Holm", "Benjamini-Hochberg"
    consistent with R - print correct_pvalues_for_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1])
    """
    from numpy import array, empty
    pvalues = array(pvalues)

    #sort p vlaues and prepare unsort index
    sort_index = np.argsort(pvalues)
    pvalues= pvalues[sort_index]
    unsort_index = np.argsort(sort_index)

    n = pvalues.shape[0]
    new_pvalues = empty(n)
    n= float(n)


    if correction_type == "Bonferroni":
        new_pvalues = n * pvalues
    elif correction_type == "Bonferroni-Holm":
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
        values.sort()
        for rank, vals in enumerate(values):
            pvalue, i = vals
            new_pvalues[i] = (n-rank) * pvalue
    elif correction_type == "Benjamini-Hochberg":
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
        values.sort()
        values.reverse()
        new_values = []
        for i, vals in enumerate(values):
            rank = n - i
            pvalue, index = vals
            new_values.append((n/rank) * pvalue)
        for i in range(0, int(n)-1):
            if new_values[i] < new_values[i+1]:
                new_values[i+1] = new_values[i]
        for i, vals in enumerate(values):
            pvalue, index = vals
            new_pvalues[index] = new_values[i]


    return new_pvalues[unsort_index]


from matplotlib.font_manager import FontProperties
def label_sig_(x1,x2,y,text='*',ax=None,fontoffset=0,linewidth=1.5,font_size=None):

    if font_size is None:
        font_size='small'

    font=FontProperties(size=font_size)
    if ax is None:
        ax=plt.gca()


    ax.annotate(text, xy=((x1+x2)/2.,y),textcoords='offset points',xytext=(0,fontoffset),
                horizontalalignment='center',verticalalignment='bottom',font_properties=font)
    ax.plot((x1,x2),(y,y),marker=3,markeredgewidth=linewidth, linewidth=linewidth,color='k')


def stats_test_all_on_once_(values1,values2,test,**test_kws):
        ResultsDB= pd.DataFrame(index=values1.columns)
        res= test(values1,values2,**test_kws)
        ResultsDB['Statistic']= res.statistic
        ResultsDB['Pvalue']= res.pvalue
        return ResultsDB

def stats_test_per_value_(values1,values2,
                   test,**test_kws):



    Pairwise_comp=pd.DataFrame(columns=['Statistic','Pvalue'],index=values1.columns)

    for variable in values1.columns:
        try:
            Pairwise_comp.loc[variable] = test(values1[variable],values2[variable],**test_kws)
        except ValueError:
            Pairwise_comp.loc[variable] = np.nan,np.nan

    return Pairwise_comp

def Paired_StatsTest(data,
                        grouping_variable,
                        Comparisons=None,
                        ref_group=None,
                        test='welch', test_kws=None,
                        correct_for_multiple_testing=True):

    """test: a parwise statistical test found in scipy e.g ['mannwhitneyu','ttest_ind']
        or a function wich takes two argumens. Additional keyword arguments can be specified by test_kws"""

    if test_kws is None:
        test_kws={}

    if type(test)==str:

        if test=='welch':
            test_kws.update(dict(equal_var=False))
            test='ttest_ind'

        assert hasattr(stats,test), "test: {} is not found in scipy".format(test)
        scipy_test= getattr(stats,test)

        if test in ['ttest_ind']:

            Test= lambda values1,values2: stats_test_all_on_once_(values1,values2,scipy_test,**test_kws)

        else:
            Test= lambda values1,values2: stats_test_per_value_(values1,values2,scipy_test,**test_kws)

    elif callable(test):
        Test=test
    else:
        Exception("Test should be a string or a callable function got {}".format(type(test)))


    data= pd.DataFrame(data)

    Groups=np.unique(grouping_variable)




    if ref_group is not None:
        Combinations=[(ref_group,other) for other in Groups[Groups!=ref_group]]

    elif Comparisons is not None:
        Combinations=Comparisons
        #TODO: check if ok
    else:
        from itertools import combinations
        Combinations= list(combinations(Groups,2))



    Results={}
    for group1,group2 in Combinations:


        values1=data.loc[ grouping_variable==group1,:]
        values2=data.loc[ grouping_variable==group2,:]



        Pairwise_comp= Test(values1,values2)


        Pairwise_comp['median_diff'] = values2.median()-values1.median()

        if correct_for_multiple_testing:
            Pairwise_comp['pBH']= Pairwise_comp[['Pvalue']].dropna().apply(correct_pvalues_for_multiple_testing, axis=0,correction_type='Benjamini-Hochberg')

        Results[group2+'_vs_'+group1]=Pairwise_comp

    Results= pd.concat(Results,axis=1)
    Results.columns = Results.columns.swaplevel(0,-1)
    Results.sort_index(axis=1,inplace=True)

    return Results


def _plot_sig_labels_hue(P_values,order,y0,deltay,ax,x_offset,show_ns=True,width=0.8,
                        Stars=True,labelkws=None):

    if not show_ns:
        P_values=P_values.loc[P_values<0.05]

    P_values=P_values.apply(format_p_value, Stars=Stars)

    if labelkws is None:
        labelkws={}

    y=y0
    for idx,text in P_values.iteritems():

        def calculate_hue_offset(group,order):
            return (order.index(group) - len(order)*0.5+0.5)/len(order)*width

        x1 = calculate_hue_offset(idx.split('_vs_')[0],order)+ x_offset
        x2 = calculate_hue_offset(idx.split('_vs_')[1],order)+x_offset



        label_sig_(x1,x2,y,text,ax=ax,**labelkws)

        y+=deltay


def _plot_sig_labels_xaxis(P_values,order,y0,deltay,
                        ax,show_ns=True,width=0.8,
                        Stars=True,labelkws=None):



    if not show_ns:
        P_values=P_values.loc[P_values<0.05]

    P_values=P_values.apply(format_p_value, Stars=Stars)

    if labelkws is None:
        labelkws={}

    y=y0
    for idx,text in P_values.iteritems():

        def calculate_x_offset(group,order):
            return order.index(group)

        x1 = calculate_x_offset(idx.split('_vs_')[0],order)
        x2 = calculate_x_offset(idx.split('_vs_')[1],order)



        label_sig_(x1,x2,y,text,ax=ax,**labelkws)

        y+=deltay



def _plot_all_sig_labels(P_values,order,orderG=None,show_ns=True,
                        y0='auto',deltay='auto',
                        ax=None,**kws):


    if ax is None:
        ax=plt.gca()

    Lim=ax.get_ylim()
    if deltay == 'auto':
        deltay=(Lim[1]-Lim[0])/10
    if y0 == 'auto':
        y0=Lim[1]+(Lim[1]-Lim[0])/10


    if orderG is None:
        _plot_sig_labels_xaxis(P_values,order,ax=ax,
                               deltay=deltay, y0=y0,show_ns=show_ns,
                            **kws)
    else:

        for cat,P in  P_values.groupby(level= list(range(len(P_values.index.levels)-1))):

                x_offset= orderG.index(cat)

                _plot_sig_labels_hue(P_values.loc[cat],order,ax=ax,
                                    x_offset=x_offset,deltay=deltay, y0=y0,show_ns=show_ns,
                                    **kws)



def clr(Data, min_multiply_factor=0.65,log= np.log2):

    min_nz_value= Data[Data>0].min().min()

    data= Data.astype(float)
    data[data==0.] = min_nz_value* 0.65

    data= log(data)
    data = (data.T-data.mean(1)).T

    return data
#data= pd.DataFrame(composition.clr(composition.multiplicative_replacement(data)),\n",
#                index=data.index,columns= data.columns)\n",


class Viewpoint():
    def __init__(self, data,sample_data=None,
                 feature_data=None,
                     test_variable=None,
                     order_grouping=None,
                     grouping_variables=None,
                     order_test=None,
                     colors=None
                ):


        self.data=DataFrame(data)

        if sample_data is None:
            self.samples = DataFrame(index=data.index,data=data.index)
        else:
            self.samples = DataFrame(sample_data.loc[data.index])
        if feature_data is None:
            self.features = pd.Series(self.data.columns,self.data.columns)
        else:
            self.features = pd.Series(feature_data.loc[data.columns])

        if test_variable is None:
            self.test_variable = self.samples.index
        elif (sample_data is not None) and (test_variable in self.samples.columns):
            self.test_variable =  self.samples[test_variable]
        else:
            self.test_variable =  test_variable

        if grouping_variables is not None:

            self.grouping_variables = self.samples[grouping_variables]
        else:
            self.grouping_variables = None

        self.stats =None

        if order_grouping is None:
            if self.grouping_variables is not None:
                self.order_grouping= np.unique(self.grouping_variables)
            else:
                self.order_grouping=None
        else:
            self.order_grouping= order_grouping
        if order_test is None:
            self.order_test= np.unique(self.test_variable)
        else:
            self.order_test= order_test

        if colors is None:
            self.colors= sns.color_palette('deep',len(order_test))
        else:
            assert len(colors)== len(self.order_test)
            self.colors=colors

    def __repr__(self):

        return f"MetaTable with {self.data.shape[0]} samples x {self.data.shape[1]} features\n"+\
            repr(self.data.iloc[:5,:5])+"\n\nSample annotations:\n"+\
            repr(self.samples.head())+"\n\nFeature annotations:\n"+\
            repr(self.features.head())

    def _apply_to_subsets(self,function,**kws):

        assert self.grouping_variables is not None

        results={}

        for subset, subset_data in self.data.groupby(self.grouping_variables):
            results[subset] = function(subset_data,**kws)
        return results

    def calculate_stats(self,comparisons=None,ref_group=None,test='welch',**kws):

        if comparisons is not None:
            for c in comparisons:
                assert c[0] in self.order_test, f"{c[0]} is not in the test variable or order_test"
                assert c[1] in self.order_test, f"{c[1]} is not in the test variabl or order_test"


        function= Paired_StatsTest

        kws= dict( grouping_variable= self.test_variable,
                  ref_group=ref_group,
                  Comparisons=comparisons,
                  test=test, **kws)


        if self.grouping_variables is None:
            results= function(self.data,**kws)
            #raise NotImplementedError()
        else:
            results= self._apply_to_subsets(function,**kws)

            results= pd.concat(results,axis=1)
            results.columns=results.columns.swaplevel(0,1)
            results.sort_index(axis=1,inplace=True)
        self.stats= results


    def boxplot(self,variable,distance_between_sig_labels='auto',box_params=None,swarm_params=None,ax=None,**labelkws):

        if ax is None:
            ax= plt.subplot(111)

        params= dict(y=self.data[variable],ax=ax)

        if self.grouping_variables is None:
            params.update(dict(x=self.test_variable,order=self.order_test))
        else:
            params.update(dict(x=self.grouping_variables,hue= self.test_variable,
                     hue_order=self.order_test,order=self.order_grouping))

        if box_params is None: box_params={}
        if swarm_params is None: swarm_params={}

        sns.boxplot(palette=self.colors,**params,**box_params)

        legend= ax.get_legend_handles_labels()

        sns.swarmplot(**params,color='k',dodge=True,**swarm_params)
        if self.grouping_variables is not None:
            ax.legend(*legend,bbox_to_anchor=(1,1))

        ax.set_ylabel(self.features.loc[variable])


        if (self.stats is not None) :
            P_values= self.stats.Pvalue.loc[variable].T
            _plot_all_sig_labels(P_values,self.order_test,self.order_grouping,
                                    deltay=distance_between_sig_labels,
                                    ax=ax,**labelkws)
