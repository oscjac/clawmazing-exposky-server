import pandas as pd
import numpy as np

# Read the CSV file and replace any all-empty columns with 0
df = pd.read_csv('calculations.csv').fillna(0)

# Perform the groupby operation and calculate the mean
merged_df = df.groupby(['pl_name', 'hostname']).agg(
    pl_rade=('pl_rade', lambda x: np.nanmean(x)), 
    pl_masse=('pl_masse', lambda x: np.nanmean(x)),
    pl_orbincl=('pl_orbincl', lambda x: np.nanmean(x)),
    st_teff=('st_teff', lambda x: np.nanmean(x))
).reset_index()

for index, row in merged_df.iterrows():
    row_dict = row.to_dict()  # Convert the row to a dictionary
    print(f"Index: {index}, Row data: {row_dict}")
