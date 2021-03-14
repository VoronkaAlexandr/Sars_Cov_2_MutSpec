import pandas as pd

id_table = pd.read_csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/ideal_table.csv')
df = pd.read_csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data/modernized_nextstrain.csv')

id_table['MutExp'] = 1
id_table['MutObsNC'][id_table['MutObsNC'].isna()] = 0
id_table['Time'] = None
id_table['Country'] = None
id_table['Continent'] = None


for index, row in df.iterrows():
    a = id_table.loc[(id_table.Pos == df['Pos'][index]) & (id_table['RefNuc'] == df['RefNuc'][index]) & (id_table['AltNuc'] == df['AltNuc'][index])]
    a = id_table['Time'][a.index.tolist()]
    if id_table['Time'][a.index.tolist()].isna:
        id_table['Time'][a.index.tolist()] = df['Time'][index]
    else:
        id_table['Time'][a.index.tolist()] = id_table['Time'][a.index.tolist()] + df['Time'][index]
    if id_table['Country'][a.index.tolist()].isna:
        id_table['Country'][a.index.tolist()] = df['Country'][index]
    else:
        id_table['Country'][a.index.tolist()] = id_table['Country'][a.index.tolist()] + df['Country'][index]
    if id_table['Continent'][a.index.tolist()].isna:
        id_table['Continent'][a.index.tolist()] = df['Continent'][index]
    else:
        id_table['Continent'][a.index.tolist()] = id_table['Continent'][a.index.tolist()] + df['Continent'][index]

id_table.to_csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/ideal_table+nextstrain.csv')

