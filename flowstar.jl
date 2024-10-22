using TaylorSeries, IntervalArithmetic

# TODO: add the partial derivs wrt time?
function lie_derivative(f, g)
    TaylorSeries.jacobian(g) * f
end
