

ode = ODE(1, x -> x, [1.])
ode = PODE(1, (x, y) -> x, (x, y) -> 2y, [1.], [1.])
dae = DAE(2, 1, x -> x, (x, λ) -> [λ -λ], x -> x[2]-x[1], [1., 1.], [0.])
dae = PDAE(1, 1, (x, y) -> x, (x, y) -> 2y, (x, y, λ) -> λ, (x, y, λ) -> -λ, (x, y) -> x[1]-y[1], [1.], [1.], [0.])
