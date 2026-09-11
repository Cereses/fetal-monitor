import numpy as np, itertools
u = np.loadtxt('captures/contractions_20260909_114206_uc4hz.csv', delimiter=',', skiprows=1)
t, v = u[:,0]-u[0,0], u[:,1]
z = v == 0.0
print(f'zeros {z.sum()} of {v.size} ({100*z.mean():.1f}%)')
i = 0
for k, g in itertools.groupby(enumerate(z), key=lambda p: p[1]):
    g = list(g); n = len(g)
    if k: print(f'  zero run {t[g[0][0]]:7.1f} to {t[g[-1][0]]:7.1f} s  ({n/4:.1f} s)')
print(f'cycle 7 press ends at 1115 s. Non-zero median {np.median(v[~z]):.0f}')