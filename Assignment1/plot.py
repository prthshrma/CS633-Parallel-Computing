import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv('data.txt',sep="\t",header=0)
arr = df.to_numpy()

processes = [16,36,49,64]

multiple = [[],[],[],[],[],[],[]]
packed = [[],[],[],[],[],[],[]]
derived = [[],[],[],[],[],[],[]]


for p in processes:
    for id, A in enumerate(arr):
        if A[3] == p:
            if A[1] == 16:
                if id % 3 == 0:
                    multiple[0].append(np.log2(A[2]))
                if id % 3 == 1:
                    packed[0].append(np.log2(A[2]))
                if id % 3 == 2:
                    derived[0].append(np.log2(A[2]))
            if A[1] == 32:
                if id % 3 == 0:
                    multiple[1].append(np.log2(A[2]))
                if id % 3 == 1:
                    packed[1].append(np.log2(A[2]))
                if id % 3 == 2:
                    derived[1].append(np.log2(A[2]))
            if A[1] == 64:
                if id % 3 == 0:
                    multiple[2].append(np.log2(A[2]))
                if id % 3 == 1:
                    packed[2].append(np.log2(A[2]))
                if id % 3 == 2:
                    derived[2].append(np.log2(A[2]))
            if A[1] == 128:
                if id % 3 == 0:
                    multiple[3].append(np.log2(A[2]))
                if id % 3 == 1:
                    packed[3].append(np.log2(A[2]))
                if id % 3 == 2:
                    derived[3].append(np.log2(A[2]))
            if A[1] == 256:
                if id % 3 == 0:
                    multiple[4].append(np.log2(A[2]))
                if id % 3 == 1:
                    packed[4].append(np.log2(A[2]))
                if id % 3 == 2:
                    derived[4].append(np.log2(A[2]))
            if A[1] == 512:
                if id % 3 == 0:
                    multiple[5].append(np.log2(A[2]))
                if id % 3 == 1:
                    packed[5].append(np.log2(A[2]))
                if id % 3 == 2:
                    derived[5].append(np.log2(A[2]))
            if A[1] == 1024:
                if id % 3 == 0:
                    multiple[6].append(np.log2(A[2]))
                if id % 3 == 1:
                    packed[6].append(np.log2(A[2]))
                if id % 3 == 2:
                    derived[6].append(np.log2(A[2]))

    plt.figure()

    M_median = []
    P_median = []
    D_median = []

    for m in multiple:
        M_median.append(np.median(m))

    for m in packed:
        P_median.append(np.median(m))

    for m in derived:
        D_median.append(np.median(m))

    plt.boxplot(multiple, labels=['16^2', '32^2', '64^2', '128^2', '256^2', '512^2', '1024^2'])
    plt.boxplot(packed, labels=['16^2', '32^2', '64^2', '128^2', '256^2', '512^2', '1024^2'])
    plt.boxplot(derived, labels=['16^2', '32^2', '64^2', '128^2', '256^2', '512^2', '1024^2'])

    x = [1, 2, 3, 4, 5, 6, 7]
    y = M_median
    plt.plot(x, y, linewidth=1, color='y', label='Multiple MPI_Send')

    y = P_median
    plt.plot(x, y, linewidth=1, color='g', label="MPI_Packed")

    y = D_median
    plt.plot(x, y, linewidth=1, color='r', label="MPI_Derived")

    plt.ylabel("Time in seconds [log2(time)]")
    plt.xlabel("N (data points per process)")

    plt.title("BoxPlot for processes "+str(p))
    plt.legend()
    plt.savefig("plot"+str(p))
    #plt.show()















