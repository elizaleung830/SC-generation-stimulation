import csv
from io import StringIO
import numpy as np
array = np.load("data_h_0.25_w_0.7_no.npz")
a = np.vstack((array["wls"], array["neff_list_te"]))
a = a.T

with open('TE_no', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(a)