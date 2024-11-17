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
	# FIXME: offset x by +5 and y by +2, since the TaylorModels
	# library requires 0 to be in the domains and such.
    du[1] = 1.5u[1] - u[1] * u[2] + 5
    du[2] = -3 * u[2] + u[1] * u[2] + 2
end

function trucktrailer!(du, u, p, t)
    # truck and trailer Parameters
    L0 = 0.3375  # truck length
    M0 = 0.1     # truck distance to center of mass
    L1 = 0.3     # trailer length
    M1 = 0.06    # trailer distance to center of mass

    # state variables: x1 (trailer x), y1 (trailer y), theta0 (truck yaw), theta1 (trailer yaw)
    # control inputs (v0, dtheta0) are now treated as additional states with constant derivatives (i.e., as controls)
    x1, y1, theta0, theta1, v0, dtheta0 = u

    # intermediate variables
    beta01 = theta0 - theta1  # angle between truck and trailer
    v1 = v0 * cos(beta01) + M0 * sin(beta01) * dtheta0  # trailer velocity

    # ODEs + zero(u[.]) for type stability
    du[1] = v1 * cos(theta1) + zero(u[1])
    du[2] = v1 * sin(theta1) + zero(u[2])
    du[3] = dtheta0 + zero(u[3])
    du[4] = (v0 / L1) * sin(beta01) - (M0 / L1) * cos(beta01) * dtheta0 + zero(u[4])

    du[5] = zero(u[5])  # v0 remains constant
    du[6] = zero(u[6])  # dtheta0 remains constant

    return du
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