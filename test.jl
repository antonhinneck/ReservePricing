using JuMP

set = [i for i in 1:10]

m = Model()

@variable(m, Ap[set])
@variable(m, Am[set])
@variable(m, p)

@expression(m, tp, p .* Ap)
@expression(m, tm, p .* Am)
