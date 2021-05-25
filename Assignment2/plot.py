#!/usr/bin/env python
# coding: utf-8

# In[47]:


import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt


df = pd.read_csv('data.txt',sep="\t",header=0)
arr = df.to_numpy()

sns.set()

#print(arr)

Bcast = []
Gather = []
Reduce = []
Alltoallv = []

for d in arr:
    if d[0]=="Bcast":
    	Db = []
    	for dataB in d[1:]:
    	    Db.append(dataB)
    	Bcast.append(Db)
    
    if d[0]=="Gather":
    	Dg = []
    	for dataG in d[1:]:
    	    Dg.append(dataG)
    	Gather.append(Dg)
    
    if d[0]=="Reduce":
    	Dr = []
    	for dataR in d[1:]:
    	    Dr.append(dataR)
    	Reduce.append(Dr)
    
    if d[0]=="Alltoallv":
    	Da = []
    	for dataA in d[1:]:
    	    Da.append(dataA)
    	Alltoallv.append(Da)
 
allData = [Bcast,Gather,Reduce,Alltoallv]
#print(Gather)



for id,Type in enumerate(allData):
    demo_input_format = pd.DataFrame.from_dict({
    "D": [],
    "P": [],
    "ppn": [],
    "mode": [],  # 1 --> optimized, 0 --> standard
    "time": [],
    })
    for Arr in Type:
    	demo_input_format = demo_input_format.append({"D": Arr[0], "P": Arr[1], "ppn": Arr[2], "mode": Arr[3], "time": Arr[4]
                }, ignore_index=True)
        
    demo_input_format["(P, ppn)"] = list(map(lambda x, y: ("(" + x + ", " + y + ")"), map(str, demo_input_format["P"]), map(str,       
    demo_input_format["ppn"])))
    
    sns.catplot(x="(P, ppn)", y="time", data=demo_input_format, kind="box", col="D", hue="mode")
    name=""
    if id==0:
        name="Bcast"
    if id==1:
        name="Gather"
    if id==2:
        name="Reduce"
    if id==3:
        name="Alltoallv"
    
    #plt.title("BoxPlot for Bcast")
    plt.legend()
    plt.savefig("plot_"+name)
    #plt.show()

