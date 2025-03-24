""" Multivariate Taylor models in N independent variables with absolute
    remainder bounds.

    This module adopts a lot of code from the TaylorModels.jl library by
    Louis Benet and David P. Sanders. We do so for the following reasons.
        1. Currently, the TaylorModels.jl library has low test coverage.
        2. The TaylorModels.jl dependency prevents us from bumping the
           TaylorSeries and IntervalArithmetic libraries to the newest
           versions.
        3. We require only a subset of the TaylorModels.jl functionality,
           mostly related to the TaylorModelN type. We opt to adapt the
           subset of needed features and simplify where feasible to ensure
           we fully understand what operations and algorithms the code
           executes. This contributes to better understanding when writing
           the paper that accompanies this code.
        4. This allows full control to modify the code.
"""

using TaylorSeries, IntervalArithmetic

import Base: iterate, size, eachindex, firstindex, lastindex,
    length, getindex, show, promote

import Base: ==, +, -, *, /, ^, in, zero, one,
    inv, sqrt, exp, log, sin, cos, tan,
    asin, acos, atan, sinh, cosh, tanh

import TaylorSeries: pretty_print, get_order, evaluate, constant_term,
    linear_polynomial, nonlinear_polynomial


include("constructors.jl")
include("promotion.jl")
include("auxiliary.jl")
include("arithmetic.jl")
include("evaluate.jl")
include("bounds.jl")
include("rpa_functions.jl")
include("show.jl")
