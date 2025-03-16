""" Multivariate Taylor model displaying in N independent variables.

    This module heavily draws from the TaylorModels.jl library by
    Louis Benet and David P. Sanders.
"""

function show(io::IO, a::TaylorModelN)
    if TS._show_default[end]
        return Base.show_default(io, a)
    else
        return print(io, pretty_print(a))
    end
end

"""TaylorModelN pretty print.

    To disable/enable the use of big-O notation, see
    `TaylorSeries.displayBigO`.
"""
function pretty_print(a::TaylorModelN)
    return string(polynomial(a), " + ", a.rem)
end
