import numpy as np 
import pandas as pd

DB = pd.read_csv("DB.csv")
print(DB.shape)

DB_small = pd.read_csv("processed_output.csv")
print(DB_small.shape)
