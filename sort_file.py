import sys
import pandas as pd

pd.read_csv(sys.argv[1],sep='\t',header=None).sort_values(by=[1]).to_csv(sys.argv[2], header=None, index=None, sep='\t')
#print(pd.read_csv(sys.argv[1],sep='\t', header=None).columns)



