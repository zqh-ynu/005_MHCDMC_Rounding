import numpy as np
import math
import matplotlib.pyplot as plt

shapes = len(np.linspace(-3,3,12))
A_loc = [
    [0, 0],
    [0, 1]
]
X, Y = np.meshgrid(np.linspace(-3,3,12), np.linspace(-3,3,12))
X = X.flatten()
Y = Y.flatten()
CLO_A = np.zeros(X.shape,dtype=np.int32)
D = np.zeros((len(A_loc), shapes * shapes))
for i in range(len(A_loc)):
    a = A_loc[i]
    x = a[0]
    y = a[1]
    D[i] = np.sqrt(np.power(X - x, 2) + np.power(Y - y, 2))
D = D.T
print(D.shape)
for i in range(X.shape[0]):
    CLO_A[i] = np.where(D[i] == np.max(D[i]))[0]

CLO_A = CLO_A.reshape(shapes, shapes)
X = X.reshape(shapes, shapes)
Y = Y.reshape(shapes, shapes)
# 每个栅格点距离最近的服务器距离
CLO_D = np.zeros((shapes, shapes))

for i in range(shapes):
    for j in range(shapes):
        aid = CLO_A[i][j]
        a = A_loc[aid]
        xa = a[0]
        ya = a[1]

        x = X[i][j]
        y = X[i][j]

        CLO_D[i][j] = math.sqrt(math.pow(xa - x, 2) + math.pow(ya - y, 2))
print(CLO_D)
print(CLO_D.shape)

# print(type(X))
# Z = (1 - X/2 + X**5 + Y**3) * np.exp(-X**2 - Y**2)
# levels = np.linspace(np.min(Z), np.max(Z), 10)
#
# fig, ax = plt.subplots()
#
# ax.contour(X, Y, Z, levels=levels)
# plt.show()