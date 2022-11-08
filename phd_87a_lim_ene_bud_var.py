from sympy.stats import P, E, variance, Die, Normal
from sympy import Eq, simplify

X, Y = Die('X', 6), Die('Y', 6) # Define two six sided dice
Z = Normal('Z', 0, 1) # Declare a Normal random variable with mean 0, std 1

P(X>3) # Probability X is greater than 3
E(X+Y) # Expectation of the sum of two dice
variance(X+Y) # Variance of the sum of two dice
simplify(P(Z>1)) # Probability of Z being greater than 1

ti_in = Normal('ti_in', 287, 50)
heat = Normal('heat', 0.7, 0.1)
dust = Normal('dust', 0.65, 0.2)
ir = Normal('ir', 220, 20)

co = ir - (ti_in*heat + ti_in*(1-heat)*dust)
Eco = float(E(co))
VARco = variance(co)
VARco = np.sqrt(float(VARco))

# Numerical verification in case on does not trust SymPy.
num_ti_in = np.random.normal(287, 50, 1000000)
num_heat = np.random.normal(0.7, 0.1, 1000000)
num_dust = np.random.normal(0.65, 0.2, 1000000)
num_ir = np.random.normal(220, 20, 1000000)
num_co = num_ir - (num_ti_in*num_heat + num_ti_in*(1-num_heat)*num_dust)
