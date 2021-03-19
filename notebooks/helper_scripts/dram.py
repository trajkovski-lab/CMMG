import pandas as pd

try:
    from tqdm import tqdm
except ImportError:
    tqdm = list

def read_dram_annotation_summary(annotation_summary):
    
    reader= pd.ExcelFile(annotation_summary)

    descriptions= []
    annotations= []

    for sheet_name in tqdm(['carbon utilization',
     'Transporters',
     'MISC',
     'Energy',
     'Organic Nitrogen',
     'carbon utilization (Woodcroft)']):



        table= pd.read_excel(reader,sheet_name=sheet_name,index_col=0)

        description_table= table.loc[:,:'subheader']
        description_table['Category']= sheet_name


        descriptions.append( description_table)
        annotations.append(table.iloc[:,4:])

    A= pd.concat(annotations,copy=False)
    Descriptions=pd.concat(descriptions,copy=False)
    

    # make tables unique
    unique_ids= ~ A.index.duplicated(keep='first')
    A=A.loc[unique_ids]

    Descriptions= Descriptions.loc[unique_ids]

    A= A.loc[(A>0).any(1)]
    Descriptions=Descriptions.loc[A.index]

    A=A.T
    
    return A, Descriptions


