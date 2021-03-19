from __future__ import print_function
from scipy import stats
from itertools import combinations
import pandas as pd
import numpy as np






from scipy import stats

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

from scipy import stats
import warnings

def Significance_test(G,Grouping_variable,test=stats.kruskal,test_kws={}):

    Kr=pd.DataFrame(1,index=G.index,columns=['P_value'])

    for g in G.index:
        try:
            S=G.loc[g].dropna()
            groups=(gr.values for name,gr in S.groupby(Grouping_variable.loc[S.index]))
            Kr.loc[g]=test(*groups,**test_kws).pvalue
        except ValueError:
            pass
        except TypeError:
            warnings.warn("I have only one group for variable: {}\n".format(g))
            pass


    Kr['Q_value']=correct_pvalues_for_multiple_testing(Kr.P_value)
    Kr['Log_P']=-np.log10(Kr.P_value)
    Kr.sort_values('Q_value',inplace=True)
    return Kr



def test_manwitney(data,grouping_variable,Groupes1=None,Groupes2=None,Alternative='two-sided'):
    ''' gives P values for parwaise comparison
    Index: group1, group2 as all parwaise combination of Groupes1 with Groupes2
    Headers: Variables tested for'''

    if Groupes1 is None: Groupes1=np.unique(grouping_variable)
    if Groupes2 is None: Groupes2=np.unique(grouping_variable)

    Results=[]
    Group_pairs=[[],[]]

    for group1 in Groupes1:
        for group2 in Groupes2:
            if not group1==group2:
                data_=data.loc[(grouping_variable==group1)|(grouping_variable==group2)]


                Grouped=data_.groupby(grouping_variable)

                if not Alternative=='two-sided' and not group1<group2:
                    if Alternative=='greater':
                        alternative_='less'
                    elif Alternative=='less':
                        alternative_=='greater'
                    else:
                        raise Exception('Alternative not understood')
                else:
                    alternative_=Alternative

                Kr=Significance_test(data_.T,grouping_variable,
                                          test=stats.mannwhitneyu,test_kws=dict(alternative=alternative_))
                Results.append(Kr['P_value'])
                Group_pairs[0].append(group1)
                Group_pairs[1].append(group2)

    K=pd.DataFrame(Results,index=Group_pairs)
    K.index.names=['Group1','Group2']

    return K


def Paired_StatsTest(data,
                        grouping_variable,
                        Groups=None,
                        ref_group=None,
                        test='ttest_ind',
                        test_kws=None,
                        return_statistic=False):

    """test: a parwise statistical test found in scipy e.g ['mannwhitneyu','ttest_ind']
        or a function wich takes two argumens. Additional keyword arguments can be specified by test_kws"""

    if type(test)==str:
        assert hasattr(stats,test), "test: {} is not found in scipy".format(test)
        Test= getattr(stats,test)

    elif callable(test):
        Test=test
    else:
        Exception("Test should be a string or a callable function got {}".format(type(test)))


    if test_kws is None:
        test_kws={}

    data= pd.DataFrame(data)

    if Groups is None:
        Groups=np.unique(grouping_variable)
    Groups= pd.Series(Groups)



    if ref_group is None:
        from itertools import combinations
        Combinations= list(combinations(Groups,2))
    else:
        Combinations=[(ref_group,other) for other in Groups.loc[Groups!=ref_group]]

    Results=[]
    for group1,group2 in Combinations:

        Pairwise_comp=pd.DataFrame(columns=[['Statistic','P_Value'],
                                            [(group1,group2)]*2],
                                   index=data.columns)
        #Results.items.names=['Group1','Group2']

        for variable in data.columns:

            values1=data.loc[ grouping_variable==group1,variable].dropna()
            values2=data.loc[ grouping_variable==group2,variable].dropna()

            try:
                Pairwise_comp.loc[variable] = Test(values1,values2,**test_kws)
            except ValueError:
                Pairwise_comp.loc[variable] = np.nan,np.nan



        Results.append(Pairwise_comp)


    Results= pd.concat(Results,axis=1)


    if return_statistic:
        return Results.P_Value[Combinations].T ,Results.Statistic[Combinations].T
    else:
        return Results.P_Value[Combinations].T



def perm_test(x, y, N_perm=1000):
    x,y=np.array(x),np.array(y)
    n, k = len(x), 0
    diff = np.abs(np.mean(x) - np.mean(y))
    z = np.concatenate([x, y])
    for j in range(N_perm):
        np.random.shuffle(z)
        k += diff < np.abs(np.mean(z[:n]) - np.mean(z[n:]))

    if k==0:
        k=1

    return k / float(N_perm)

## Categorical variables

def get_count_matrix(variable,grouping_variable):
    variable= pd.Series(variable)
    return variable.groupby([variable,grouping_variable]).size().unstack().fillna(0)


def Chi2_test(variable,grouping_variable):

    Count_matrix=get_count_matrix(variable,grouping_variable)
    chi2,p_value,doff,expected=Results=stats.chi2_contingency(Count_matrix)

    return chi2,p_value


def Chi2_test_for_data(data,grouping_variable,cathegorical_colums,order=None):

    if order is None:
        order = list(grouping_variable.unique())

    if cathegorical_colums is None:
        cathegorical_colums= data.columns

    def concat(values):
        return ':'.join(values.astype('str'))

    Ratios= pd.DataFrame(index=cathegorical_colums,columns=order+['Description'])
    Statistics=pd.DataFrame(index=cathegorical_colums,columns=['Chi2','P_value'])

    for c in cathegorical_colums:
        Count_matrix=get_count_matrix(data[c],grouping_variable)
        Ratios.loc[c,'Description']=concat(Count_matrix.index.values)
        Ratios.loc[c,order]=Count_matrix.astype(int).apply(concat).loc[order].T

        chi2,p_value,doff,expected=stats.chi2_contingency(Count_matrix)
        Statistics.loc[c]= chi2,p_value

    Out=Ratios.join(Statistics)
    Out.columns=[['Ratios']*(len(order)+1)+['Statistics']*2,
                 Out.columns]

    return Out

def Odds_Ratio(variable,grouping_variable):
    Count_matrix=get_count_matrix(variable,grouping_variable)
    assert Count_matrix.shape==(2,2)
    C=Count_matrix.T.as_matrix()
    OR= float(C[0,0]*C[1,1])/(C[1,0]*C[0,1])
    return OR

def Relative_Risk(variable,grouping_variable):
    Count_matrix=get_count_matrix(variable,grouping_variable)
    assert Count_matrix.shape==(2,2)
    C=Count_matrix.T.as_matrix()
    RR= (1.* C[0,1]/C[0,1]) /(1.*C[0,:].sum()/C[1,:].sum())
    return RR




## Correlations

def calculate_corr(Data,corr_function='spearmanr'):
    corr,pv= getattr(stats,corr_function)(Data)
    corr=pd.DataFrame(corr,index=Data.columns, columns= Data.columns)
    pv=pd.DataFrame(pv,index=Data.columns, columns= Data.columns)

    return corr,pv

def correlation_byGroup(Data,grouping_variable,order=None,corr_function='spearmanr'):
    if order is None:
        order=grouping_variable.unique()

    Correlations,Pvalues,N={},{},{}

    corr,pv = calculate_corr(Data,corr_function)
    Correlations['All']=corr
    Pvalues['All']=pv
    N['All']=Data.shape[0]

    for group in order:
        d=Data.loc[grouping_variable==group]
        corr,pv = calculate_corr(d)
        Correlations[group]=corr
        Pvalues[group]=pv
        N[group]=d.shape[0]

    Correlations=pd.Panel.from_dict(Correlations).loc[['All']+order]
    Pvalues=pd.Panel.from_dict(Pvalues).loc[['All']+order]
    N=pd.Series(N)

    return Correlations,Pvalues,N
