using SciMLBase
using OrdinaryDiffEq

"""A representation of the following 2D ODE:
    x' = 1 + y
    y' = -x^2
"""
function ode332!(du, u, p, t)
    du[1] = 1.0 + u[2]
    du[2] = -(u[1]^2)
end

"""Use forward Euler to solve the given ODE.

    @param[in]
"""
function euler(vector_field, # ::ODEFunction
               time_horizon::Float64,
               time_step_size::Float64,
               initial_state::Vector{Float64})
               u0 = [0.0;0.0]
    tspan = (0.0, time_horizon)
    prob = ODEProblem(vector_field, initial_state, tspan, dt=time_step_size)

    sol::ODESolution = solve(prob, Euler())
    return sol.t, sol.u
end