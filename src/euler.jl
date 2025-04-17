using SciMLBase
using OrdinaryDiffEq


"""A representation of the following 1D ODE:
    x' = sin(x)
"""
function ode3310!(du, u, p, t)
    du[1] = sin(u[1])
end

"""A representation of the following 2D ODE:
    x' = 1.5x - xy
    y' = -3y + xy
"""
function ode3311!(du, u, p, t)
    du[1] = 1.5u[1] - u[1] * u[2]
    du[2] = -3 * u[2] + u[1] * u[2]
end

"""A representation of the following 1D ODE:
    y' = y + cos(y)
"""
function odeTMJets!(du, u, p, t)
    du[1] = u[1] * cos(u[1])
end

"""Use forward Euler to solve the given ODE.

    @param[in]
"""
function euler(vector_field, # ::ODEFunction
               time_horizon::Float64,
               time_step_size::Float64,
               initial_state::Vector{Float64})
    tspan = (0.0, time_horizon)
    prob = ODEProblem(vector_field, initial_state, tspan, dt=time_step_size)

    #sol::ODESolution = solve(prob, Euler())
    sol::ODESolution = solve(prob, AutoTsit5(Rosenbrock23()))
    tseries, vseries = sol.t, sol.u
    return [
        collect(zip(tseries, getindex.(vseries, var_idx)))
        for var_idx in eachindex(initial_state)
    ]
end
