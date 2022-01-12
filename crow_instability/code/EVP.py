from cylinderModule import *

# shift value
sigma = 2 + 0.1j
# number of eigenmodes to compute
k = 10
kk = 1.4

print(cpp.la.linear_algebra_backends())
cylinder = LIAproblem2D(100,"./Results/DNS/DNS/")
# cylinder.steady_state()
print('teste')
#cylinder.eigenvalues(sigma=sigma, k=k, kk=kk)
cylinder.eigenvalues_k(sigma=sigma, k=k, kk_list=np.arange(0,1,0.2))
