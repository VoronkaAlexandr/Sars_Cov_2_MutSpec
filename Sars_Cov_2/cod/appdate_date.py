import pandas as pd
from datetime import datetime, date, timedelta

df = pd.read_csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data/modernized_nextstrain.csv')

time_format = []
years = []
months = []
for i in range(len(df['Time'])):
    time = df['Time'][i]
    year = int(time)
    part_of_year = time-year
    day = int(part_of_year*365)
    d = date(year, 1, 1)
    time_format.append(d+timedelta(day))
df['Time'] = time_format
df['MutCount'] = 1
for i in df.Time:
    years.append(i.year)
    months.append(i.month)

df['Year'] = years
df['Month'] = months
df.to_csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data/norm_data_modernized_nextstrain.csv')