from __future__ import print_function, division
import pandas as pd
import matplotlib.pylab as plt
import numpy as np
import warnings


from scipy import stats

def normalize_row(d):
    return d/d.sum()

def _validate(counts, suppress_cast=False):
    """Validate and convert input to an acceptable counts vector type.
    Note: may not always return a copy of `counts`!
    """
    counts = np.asarray(counts)

    if not suppress_cast:
        counts = counts.astype(int, casting='safe', copy=False)

    if counts.ndim != 1:
        raise ValueError("Only 1-D vectors are supported.")
    elif (counts < 0).any():
        raise ValueError("Counts vector cannot contain negative values.")

    return counts



def shannon(counts, base=2):
    """Calculate Shannon entropy of counts (H), default in bits.
    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.
    base : scalar, optional
        Logarithm base to use in the calculations.
    Returns
    -------
    double
        Shannon diversity index H.
    Notes
    -----
    The implementation here is based on the description given in the SDR-IV
    online manual [1]_, except that the default logarithm base used here is 2
    instead of :math:`e`.
    References
    ----------
    .. [1] http://www.pisces-conservation.com/sdrhelp/index.html
    """
    counts = _validate(counts,suppress_cast=True)
    freqs = counts / counts.sum()
    nonzero_freqs = freqs[freqs.nonzero()]
    return -(nonzero_freqs * np.log(nonzero_freqs)).sum() / np.log(base)


import os

def load_metaphlan_data(merged_table,Subset=None, normalize_after_selection=False, Level='All'):
    """ load metaphlan data from a merged table which can be created withe the script metaphlan merge table
        Relative abundance sums to 100 within each taxonomic Level.

        Parameters:

        merged_table:               path to merged table
        Subset:                     Select a subset based on a string found in all indexes e.g. k_Bacteria, k_Viruses, g_Bifidobacterium
                                    If you do a selection take yare if you want to normalize or not.
        normalize_after_selection:  Boolean if data is normalized to 100 after Subset selection.
        Level:                      one of ["Strain","Species","Genus","Family","All"] default "All"
                                    Does a subselection based of the taxonomic level"""


    M=pd.read_table(merged_table,index_col=0,comment='#')


    #remove .metagenome.txt from column name
    M.columns=M.columns.str.replace('.metagenome','')


    if Subset is not None:
        M=M.loc[M.index.str.contains(Subset)]


    # taxonomic level
    if Level=='All':
        Mg=M
    else:
        if Level=='Genus':
            Mg=M.loc[M.index.str.contains('g__')&~M.index.str.contains('s__')]
        elif Level=='Species':
            Mg=M.loc[M.index.str.contains('s__')&~M.index.str.contains('t__')]
        elif Level == 'Family':
            Mg=M.loc[M.index.str.contains('f__')&~M.index.str.contains('g__')]
        elif Level=='Strain':
            Mg=M.loc[M.index.str.contains('t__')]


        Mg.index=[e.split('|')[-1][3:] for e in Mg.index]

    if normalize_after_selection:
        Mg= Mg/Mg.sum()*100.

    return Mg


def read_b3_table(B3_folder,Level='L7',subfolder='b3_all_taxa_levels_group_significance_*/raw_data'):

    for root, dirs, files in os.walk(os.path.join(B3_folder,subfolder)):
        for name in files:
            if Level in name and 'txt' in name:
                OTU_file=os.path.join(root, name)
    if OTU_file is None:
        print('OTU_file not found')

    OTU_table=pd.read_table(OTU_file,header=0,index_col=0).T
    return OTU_table



def load_alpha_div(alphaD_table_file):
    """alphaD_table_file e.g. {B3_folder}/arare_max1000/alpha_div_collated/{alhaD_measure}.txt  """
    M=pd.read_table(alphaD_table_file,index_col=0,header=0)
    #take part with highest rarefaction
    M=M.loc[M.iloc[:,0]==M.iloc[:,0].max()].iloc[:,2:]
    D=M.median()
    D.name=os.path.splitext(os.path.basename(alphaD_table_file))[0]
    del M
    return D

def read_raw_OTU_table(OTU_file):
    """ returns taxonomy, and OTU_table
    """

    OTU_table=pd.read_table(OTU_file,header=1,index_col=0).T
    taxonomy=OTU_table.loc['taxonomy']
    OTU_table.drop('taxonomy',inplace=True)
    return taxonomy,OTU_table

def minfilter(data,value=15,mode='fixed',combine_as_Other=True):
    '''  sample ID X bacteria_taxa
        mode: ['fixed','atleast1','mean','proportion','fixedmean']
    '''


    if data.shape[0]> data.shape[1]:
        warnings.warn("May be your data is not in the right orientation sample x feature")


    if mode=='atleast1':
        comparison=data.T>value
        criteria=comparison.any(1)
    elif mode=='mean':
        criteria=data.mean()>value
    elif mode=='fixed':
        assert type(value) == int, "for this mode I expect an int: number of features selected"
        value=min(value,data.shape[1])-2
        max_values=data.max()
        criteria= (max_values>=max_values.sort_values(ascending=False).iloc[value])
    elif mode=='fixedmean':
        assert type(value) == int, "for this mode I expect an int: number of features selected"
        value=min(value,data.shape[1])-2
        max_values=data.mean()
        criteria= (max_values>=max_values.sort_values(ascending=False).iloc[value])
    elif mode=='proportion':
        assert type(value) == float, "for this mode I expect an float: proportion of cumulative sum"
        means=data.mean().sort_values(ascending=False)
        means/=means.sum()
        criteria=means.cumsum()<value
    else:
        raise Exception("Didn't understand mode: {}".format(mode))


    refined=data.loc[:,criteria].copy()

    if combine_as_Other:
        refined['Other']=data.loc[:,~criteria].sum(1)
    return refined

def core_filter(data, proportion=0.25, min_value=0, combine_as_Other=False):
    ''' sample ID x bacteria_taxa
    '''

    if data.shape[0]> data.shape[1]:
        warnings.warn("May be your data is not in the right orientation sample x feature")

    presence=(data>min_value).sum(0)/float(data.shape[0])
    criteria= presence > proportion

    refined=data.loc[:,criteria].copy()
    if combine_as_Other:
        refined['Other']=data.loc[:,~criteria].sum(1)
    return refined




def sort_by_most_aboundant_species(data):
    '''sort table by most aboundant index
    '''
    most_aboundant_species=data.mean().sort_values(ascending=False).index[0]
    return data.sort_values(most_aboundant_species).T


def Get_Spec(Taxonomy,split_character=None,Find_always_first=False):
    '''Get specification of genus family from Tax string
    e.g 'k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacteriales; f__Enterobacteriaceae; g__Escherichia; s__' '''




    def remove_taxid(tax_str):
        for tax in "kpcofgst":
            tax_str=tax_str.replace(tax+'__','')
        return tax_str



    def is_informative_taxname(tax_name):

        if tax_name=='Other' or tax_name=='other' or len(tax_name)<=4: #if smaler than 4 the name is only g__ or s__
            return False, 'other'

        elif tax_name=='sp' or tax_name=='sp.':
            return False, 'sp.'

        elif '[' in tax_name:
            return False, tax_name

        else:
            return True,tax_name

    def search_first(Taxionomy_):
        first=None
        for e in reversed(Taxionomy_[:-1]):
            is_informative,alternative_tax= is_informative_taxname(e)
            if not is_informative: continue
            else:
                first=e
                break;
        if first is None:
            raise(Exception("Didn't find first taxonomy for string: {}".format(Taxionomy_)))
        return remove_taxid(first)




    def get_Spec_(Taxstr,split_character):

        if split_character is None:
            if '|' in Taxstr: split_character= '|'
            elif ';' in Taxstr: split_character= ';'
            else: split_character=' '

        Taxionomy_=Taxstr.split(split_character)

        last=remove_taxid(Taxionomy_[-1])


        is_informative, alternative_tax= is_informative_taxname(last)

        if is_informative and not Find_always_first:
            tax_out= alternative_tax
        else:
            tax_out= search_first(Taxionomy_) + ' ' + alternative_tax



        return tax_out

    if hasattr(Taxonomy, '__iter__'):
        return [get_Spec_(i,split_character) for i in Taxonomy]
    elif type(Taxonomy) == str:
        return get_Spec_(Taxonomy,split_character)
    else :
        raise TypeError('Need string or iterable of strings got: ',type(Taxonomy))






import pandas as pd


def read_fasta_file(fasta_file):
    """reads a fasta file and puts all sequences in a pandas Series"""
    from Bio import SeqIO
    rep_seqs= SeqIO.parse(fasta_file,'fasta')
    return pd.Series(dict((s.id,str(s.seq)) for s in rep_seqs))

def table(S):
    return S.groupby(S).size()


def gen_names_for_range(N,prefix=''):

    n_leading_zeros= len(str(N))

    format_int=prefix+'{:0'+str(n_leading_zeros)+'d}'

    return [format_int.format(i) for i in range(N)]
