from skbio.stats import composition
import numpy as np
import pandas as pd


def clr(counts_data,log= np.log2):


    #TODO: check if count data

    # remove columns with all
    data= counts_data.loc[:,~(counts_data<=1).all()]

    #dataframe with replace zeros
    data= pd.DataFrame( composition.multiplicative_replacement(data),
                       columns=data.columns,
                       index= data.index
                      )

    data= log(data)
    data = (data.T-data.mean(1)).T

    return data
