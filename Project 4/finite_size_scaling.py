import numpy as np

L = np.array([20,40,60,80])				# L values
T_C = np.array([2.34,2.31,2.32,2.27])	# corresponding T_c(L) values

a = np.zeros(6)

# different ways to compute a by subtracting equations
a[0] = (T_C[3]-T_C[2])/(L[3]-L[2])
a[1] = (T_C[3]-T_C[1])/(L[3]-L[1])
a[2] = (T_C[3]-T_C[0])/(L[3]-L[0])
a[3] = (T_C[2]-T_C[1])/(L[2]-L[1])
a[4] = (T_C[2]-T_C[0])/(L[2]-L[0])
a[5] = (T_C[1]-T_C[0])/(L[1]-L[0])

exact = 2.269

print 'Different approximations of a: '
print a

a = sum(a)/6

print 'Mean value of a: ', a

print 'Finite size scaling approximation of T_c: ', T_C[3]-a*L[3]**-1
print 'Relative error: ', abs(exact-(T_C[3]-a*L[3]**-1))/exact