# Dependencies

## Important

Some of the packages require us to cite them:

1. IntervalArithmetic, citation required, [link](https://github.com/JuliaIntervals/IntervalArithmetic.jl/blob/master/CITATION.bib)
1. TaylorSeries, citation required, [link](https://github.com/JuliaDiff/TaylorSeries.jl/blob/master/CITATION.bib)
1. TaylorModels, **no citation required?**

## ToC - Recommended Changes

- [IntervalArithmetic](#recommended-changes-intervals)
- [TaylorSeries](#recommended-changes-taylor-series)
- [TaylorModels](#recommended-changes-taylor-models)

## General Julia Notes

### Evaluating Expressions

Every Julia expression is a nested structure of `Expr` objects that allows manipulation of the expression as if it were an abstract syntax tree.

Additionally, "_Every module has its own eval function that evaluates expressions in its global scope. Expressions passed to eval are not limited to returning values – they can also have side-effects that alter the state of the enclosing module's environment:_"


### Links

- Julia, Evaluating Expressions, [link](https://docs.julialang.org/en/v1/manual/metaprogramming/#Evaluating-expressions)



## IntervalArithmetic

The main datastructure is the `Interval`. Prefer `Interval` over `BareInterval`. By definition it is $[a, b] = \{x \in \mathbb{R} \mid a \le x \le b\}$. Any interval creation that violates $a \le b$, such as $[3, 1]$, is instead mapped to the empty set $\emptyset$.

As per the IntervalArithmetic docs: "_Basic arithmetic operations (+, -, *, /, ^) are defined for pairs of intervals in a standard way: the result is the smallest interval containing the result of operating with each element of each interval._"

### Recommended Changes (intervals)

- Add assertions to all code to require that `isguaranteed(i) == true` is maintained for all interval computations?
- Add assertions to all code to require that `decoration(i) == com` is maintained for all interval computations?
- Try to **bump the IntervalArithmetic version** up to the most recent version. Though this may be **blocked** by the TaylorModels package. There are differences between v0.21.x and v0.22.x:
  - `setformat(...)` was renamed to `setdisplay(...)` in v0.22
  - `isguaranteed(...)` was introduced in v0.22
  - v0.22 includes decoration functionality in the `Interval` struct, while versions < v0.22 may use a `DecoratedInterval` sub-type instead.
  - In v0.22.x you must explicitly import the shorthands `..` and `±` for the constructor `interval()` with `using IntervalArithmetic.Symbols`.

### Check the Following

- Does TaylorModels.jl make use of a recent version of IntervalArithmetic.jl? Does it consider guarantees and such?
- Which rounding "[flavor](https://github.com/JuliaIntervals/IntervalArithmetic.jl/blob/master/src/intervals/flavor.jl)" does TaylorModels.jl use?

### Rounding

Intervals work with fixed size floating point numbers, meaning rounding errors are a concern when manually constructing intervals. The IntervalArithmetic interval parsing interface tries to solve this.

```julia
# Julia suffers from parse-time rounding errors of fixed size floats.
julia> x = 0.1 # typeof(x) == Float64 
julia> x
0.1
julia> x > 1//10 # typeof(1//10) = Rational, this is a fraction without rounding error
true # x = 0.1 > 1/10! This indicates a rounding error!

# What if we construct an interval?
julia> interval(0.1) # The degenerate interval [0.1, 0.1]
Interval(0.1, 0.1)
julia> interval(0.1).lo > 1//10 # intervals suffer from parse-time rounding errors!
true # 0.1 > 1/10! This indicates a rounding error!

# We must use the parsing interface of IntervalArithmetic.
# The constructed interval takes rounding error into account.
julia> I"0.1" # The interval [a, b]
Interval(0.09999999999999999, 0.1)
julia> I"0.1".lo <= 0.1
true # a <= 0.1
julia> 0.1 >= I"0.1".hi
true # 0.1 <= b
```



### Decorations

Every `Interval` contains a decoration member that specifies some properties of the given interval. These properties have an absolute order: `com > dac > def > trv > ill` or `common > defined and continuous > defined > trivial > ill-formed`. This is useful to inspect the result of applying a function to an interval. See the [docs](https://juliaintervals.github.io/IntervalArithmetic.jl/stable/manual/construction/#Decorations).

```julia
#
# For IntervalArithmetic v0.22.23
#
julia> using IntervalArithmetic.Symbols # include the ".." interval constructor
julia> div(x) = 1 / x
julia> decoration(-1..1) # The constructor `interval()` is com by default
com::Decoration = 4
julia> decoration(div(-1..1)) # 1 / x is undefined for x=0, so the function is discontinuous at that point in the domain leading to the decoration "trivial".
trv::Decoration = 1


#
# For IntervalArithmetic < v0.22.x
#
julia> div(x) = 1 / x
julia> decoration(-1..1) # The constructor `interval()` is com by default
com::Decoration = 4
# In older IntervalArithmetic versions an Interval object is
# distinct from a DecoratedInterval object.
julia> decoration(div(-1..1)) # This does not yield "trv" as expected
dac::Decoration = 3
julia> decoration(dic(DecoratedInterval(-1, 1))) # This yields "trv" as expected
trv::Decoration = 1
```

Older versions of IntervalArithmetic may prefer the use of `DecoratedInterval` for using decorations, while the newer versions include that functionality in the primary `Interval` struct instead.


### Guarantees

Introduced in IntervalArithmetic v0.22.x.

As per the IntervalArithmetic docs: "_Specifically, an interval x constructed via interval satisfies `isguaranteed(x) == true`. However, if a call to `convert(::Type{<:Interval}, ::Real)` occurs, then the resulting interval x satisfies `isguaranteed(x) == false`, receiving the "NG" (not guaranteed) label. For instance, consider the following examples:_"

```julia
julia>  convert(Interval{Float64}, 1.) # considered "not guaranteed" as this call can be done implicitly
Interval{Float64}(1.0, 1.0, com)_NG
julia> interval(1) # considered "guaranteed" as the user explicitly constructed the interval
Interval{Float64}(1.0, 1.0, com)
```

The source code comments for `isguaranteed(...)` state: "_Test whether the interval is not guaranteed to encompass all possible numerical errors. This happens whenever an [`Interval`](@ref) is constructed using `convert(::Type{<:Interval}, ::Real)`, which may occur implicitly when mixing intervals and `Real` types._"

As per the IntervalArithmetic docs: **"_A user interested in validated numerics should always have a resulting interval for which isguaranteed is true._"**.

### Loose Evaluation

As per the IntervalArithmetic docs, "_The process of discarding parts of an interval that are not in the domain of a function is called loose evaluation._" For example, for domain $x \in [-2, 2]$ and the function sqrt, the evaluation $sqrt([-2, 2])$ discards the negative part $[-2, 0[$ of the interval and evaluates sqrt on the remaining section $[0, 2]$.

```jl
julia> sqrt(-2..2)
[0, 1.41422]
```

Taylor models exclusively contain polynomials, **so we should expect that every single one of our computations has `com` type decorations. But we should take care that this still holds w.r.t. the error accrued due to Taylorizing the dynamics/ODEs.**

Additionally, as per the IntervalArithmetic docs: "_The decoration trv is an indicator of information loss. Often this also reveals that something unexpected occured. Therefore, any interval marked by this decoration may not be trusted and the code may need to be revised._".


### Polynomial Interval Enclosure

By **section 2.2.2 _Interval evaluation for polynomial functions_** of Xin Chen's thesis, an interval extension $P$ of a polynomial $p$ can be easily obtained and is inclusion isotonic and Lipschitz (see Lemma 2.2.9 in Xin Chen's thesis). This means that an interval evaluation of a polynomial is always an over-approximation of the true range of that polynomial over the evaluated intervals.

More intuitively, note that interval arithmetic operations are overapproximate forms of their set operation counterparts; the true solution of the set operation is always contained inside the solution computed via interval arithmetic. For example, for $X = \{x \in \mathbb{R} \mid a_x \le x \le b_x\}, Y = \{y \in \mathbb{R} \mid a_y \le y \le b_y\}$ the elementwise set product is $X \cdot Y = \{x \cdot y \mid \forall x \in X, \forall y \in Y\}$ and the interval arithmetic counterpart is $X \odot Y = [a_x, b_x] \odot [a_y, b_y]$. Then $X \cdot Y \subseteq X \odot Y$.

Interval arithmetic defines such over-approximate formulas for all operations required by polynomial equations: addition, subtraction, multiplication, exponentiation and division. Polynomials are simply the composition of variables using these basic operations, so it stands to reason that **an interval evaluation of a polynomial**, i.e. replacing each variable $x$ by the corresponding interval domain $[a_x, b_x]$ where $a_x \le b_x$ and replacing each coefficient $c_x$ by the degenerate interval $[c_x, c_x]$ and finally reducing the expression via interval arithmetic, **produces an over-approximation of the range of the polynomial**, since the interval extension is inclusion isotonic and Lipschitz.

### Version Problems

I have IntervalArithmetic stuck at v0.20.9 instead of the most recent version v0.22.23 because using Guillermo's version of TaylorModels blocks this.

Check the file `NEWS.md` in various releases of the IntervalArithmetic github for these issues.

- The _Trivial_ section of the constructor docs says that `decoration(sqrt(-2..2)) == trv`, but I found this to result in decoration `dac` instead. To yield `trv`, you need ``decoration(sqrt(DecoratedInterval(-2, 2))) == trv``.

- In version v0.22 of IntervalArithmetic the function `setformat(...)` was renamed to `setdisplay(...)`. But I have to call `setformat(:full, decorations=true)` instead of `setdisplay(:full, decorations=true)`.


### Links

1. IntervalArithmetic, Github repo, [link](https://github.com/JuliaIntervals/IntervalArithmetic.jl)
1. IntervalArithmetic, main documentation, [link](https://juliaintervals.github.io/IntervalArithmetic.jl/stable)
1. IntervalArithmetic, Decorations docs, [link](https://juliaintervals.github.io/IntervalArithmetic.jl/stable/manual/construction/#Decorations)


## TaylorSeries

The documentation (which contains a User guide), the github repository and the accompanying paper `TaylorSeries.jl: Taylor expansions in one and several variables in Julia` are the main sources of information for how to use this package.

We will focus on the `TaylorN` struct as the main datastructure, since we require multivariate Taylor polynomials. The `Taylor1` and `HomogeneousPolynomial` types will be mostly ignored.


### Recommended Changes (taylor series)

- Make use of the function `taylor_expand` when shifting the Taylorized vector field?
- Replace the $t - t$ hack in Lie derivatives with the `TaylorN` constructor for a constant 0 of the desired (default) order?

### Check the Following

- Do we correctly make use of `TaylorN.order` VS `get_order()` as the order of constructed `TaylorN` objects? i.e. generally `get_order()` seems more of an upper bound for all `TaylorN` models, while a single object's `TaylorN.order` value only affects that object itself. Though the interplay of multiple `TaylorN.order` values during arithmetic must be taken into account; TaylorSeries always lower the order of multiple objects to the lowest order of all operands?
- Do all TaylorN polynomials part of the TaylorModelN objects have the exact same order? Do they have to? See the [evaluation](#evaluation) section and the test code, where we sometimes run into issues when evaluating if the `TaylorN` order of the inputs differ.
- Does the _mixtures_ TODO next to the definition of `power_by_squaring!` in `TaylorSeries.power.jl` prevent us from calling $evaluate(p, [v_1, ..., v_n])$ for $v_1, \dots, v_n$ being of type `TaylorN` and having different orders, i.e. $\exist v_i, v_j : i\ne j \land get\_order(v_i) \ne get\_order(v_j)$? I think it does, but I'm not sure. See the test code for an example that fails when evaluating using mixed orders for the inputs. **This could be a problem when we allow the polynomials to change order dynamically, individually, independently from the result returned by `get_order()`.**
- The TaylorSeries user guide states that "_To evaluate a Taylor series at a given point, Horner's rule is used via the function `evaluate(a, dt)`._" Is this equivalent with transforming a polynomial into Horner normal form pre evaluation?

### TaylorN Design

The datastructure of interest is `TaylorN`. It represents a multivariate Taylor polynomial expanded at the origin. So the expansion points are $c_1 = 0, \dots, c_n = 0$ for an $n$-variate polynomial. This is because the polynomial does not explicitly store the expansion point; **the expansion point must be manually added into the polynomial itself**. `TaylorN` contains two members:

- coeffs: A vector where the $i$-th element is a `HomogeneousPolynomial` of order $i+1$, i.e. the polynomial containing all (and only) terms of order $i+1$.
- order: The maximum order of the Taylor polynomial. i.e. All terms of order greater than this are effectively truncated from the polynomial.

**Note that the `order` member is essentially the truncation degree of the polynomial, it is *not* the degree of the currently largest term actually part of the polynomial.** For example, a `TaylorN` can have order $10$ and represent a polynomial $p = x \cdot y^2$ that has actual degree $3 \le 10$.

This difference matters for Taylor series arithmetic. Computing $q = p^3 = x^3 \cdot y^6$ results in a polynomial of actual degree $9 \le 10$, so no truncation happens. But $q = p^4 = x^4 \cdot y^8$ has actual degree $12 \gt 10$, so $q = 0$ is all that remains after truncation.


<details>
  <summary>TaylorN struct source code + doc.</summary>

```julia
"""
    TaylorN{T<:Number} <: AbstractSeries{T}

DataType for polynomial expansions in many (>1) independent variables.

**Fields:**

- `coeffs  :: Array{HomogeneousPolynomial{T},1}` Vector containing the
`HomogeneousPolynomial` entries. The ``i``-th component corresponds to the
homogeneous polynomial of degree ``i-1``.
- `order   :: Int`  maximum order of the polynomial expansion.

Note that `TaylorN` variables are callable. For more information, see
[`evaluate`](@ref).
"""
struct TaylorN{T<:Number} <: AbstractSeries{T}
    coeffs  :: Array{HomogeneousPolynomial{T},1}
    order   :: Int

    function TaylorN{T}(v::Array{HomogeneousPolynomial{T},1}, order::Int) where T<:Number
        coeffs = isempty(v) ? zeros(HomogeneousPolynomial{T}, order) : zeros(v[1], order)
        @inbounds for i in eachindex(v)
            ord = v[i].order
            if ord ≤ order
                coeffs[ord+1] += v[i]
            end
        end
        new{T}(coeffs, order)
    end
end
```

</details>


TaylorSeries internally uses a singleton (const, global) object, which can be accessed as `TaylorSeries._params_TaylorN_`, of type `ParamsTaylorN` to store information on how to construct objects of type `TaylorN`. The docs + definition for this struct follow.

<details>
  <summary>ParamsTaylorN struct source code + doc.</summary>

```julia
"""
    ParamsTaylorN

DataType holding the current parameters for `TaylorN` and
`HomogeneousPolynomial`.

**Fields:**

- `order            :: Int`  Order (degree) of the polynomials
- `num_vars         :: Int`  Number of variables
- `variable_names   :: Vector{String}`  Names of the variables
- `variable_symbols :: Vector{Symbol}`  Symbols of the variables

These parameters can be changed using [`set_variables`](@ref)
"""
mutable struct ParamsTaylorN
    order            :: Int
    num_vars         :: Int
    variable_names   :: Vector{String}
    variable_symbols :: Vector{Symbol}
end
```

</details>


The members of the singleton can be accessed using parameterless helper functions:
- order: `get_order()`
- num_vars: `get_numvars()`
- variable_names: `get_variable_names()`
- variable_symbols: `get_variable_symbols()`

The function `set_variables` exposes many parameter combinations to update the members of the singleton. It is the primary way to update the members.

**Note that TaylorSeries does not store the variables as TaylorN objects. It only constructs the actual TaylorN objects for each variable as specified by `TaylorSeries._params_TaylorN_` in response to a call of `get_variables` or `set_variables`. This means the objects returned by subsequent calls to these functions are distinct, and different, existing TaylorN objects are not implicitly modified after calls to those functions.**


### Stored Maximum Order

The TaylorSeries package internally stores a **maximum order**, which all TaylorN polynomials must adhere to. This is the `order` member of `TaylorSeries._params_TaylorN_`, which is the value returned by calling `get_order()`. It is used as the default order when constructing `TaylorN` objects, and as the truncation order for Taylor polynomial arithmetic.

This affects many operations.

- Given `TaylorN` objects $a$ and $b$, many arithmetic operations use `fixorder(a, b)` to force $a$ and $b$ to take on order $min(order(a), order(b))$ before applying the operation. In other words, operations that combine different orders propagate the lowest order. This results in truncation for at most one of the polynomials!
- You can only specify $order \le get\_order()$ for the function `get_variables(order::Int)`.

#### fixorder

We are mostly interested in the function definition of `fixorder(a, b)` for $a$ and $b$ being a `TaylorN` or `Taylor1`, which is [found in](https://github.com/JuliaDiff/TaylorSeries.jl/blob/0298820a6d1f903185e20849c8178ffb0c1cd503/src/auxiliary.jl#L259) `TaylorSeries/auxiliary.jl`.

It returns a tuple $(a', b')$ where $a'$ and $b'$ are copied objects that correspond to the `TaylorN` objects $a$ resp. $b$, but with their order restricted to $min(get\_order(a), get\_order(b))$. For example,

$$
\begin{align*}
  a  &= x + x^2 + x^3 & \text{with get\_order(a) = 4}\\
  b  &= y + y^2       & \text{with get\_order(b) = 2}\\
  a' &= x + x^2       & \text{with get\_order(a') = 2}\\
  b' &= y + y^2       & \text{with get\_order(b') = 2}\\
\end{align*}
$$


### Elementary Functions

TaylorSeries is able to dynamically generate the Taylor expansion of many elementary functions, such as exp, sin, etc. They dynamically generate these expansions. For example, the expansion of $sin(2x)$ is generated by first constructing a `TaylorN` expansion $t_1$ for $2x$ and calling $sin(t_1)$ to generate the final expansion $t_2$. They make use of function overloading to provide a seamless API.


### Evaluation

A $n$-variate polynomial $p$ of type `TaylorN`, which is implicitly expanded at the origin $c_1 = 0, \dots, c_n = 0$, can be evaluated at $x_1 = v_1, \dots, v_n = c_n$ using either of the following, equivalent function calls: $evaluate(p, [v_1, \dots, v_n])$ or $p([v_1, \dots, v_n])$. Thanks to Julia's Expr representation of statements, $v_i$ can be of any type as long as it defines the required arithmetic operators (+, -, etc.) and or elementary functions. To evaluate at $x_1 = 0, \dots, x_n = 0$ just call $evaluate(p)$ or $p()$.

**Note that the TaylorSeries user guide docs claim that "_To evaluate a Taylor series at a given point, Horner's rule is used via the function `evaluate(a, dt)`._" This is equivalent with transforming the polynomial into Horner normal form before evaluating it, to minimize dependency problem we will encounter w.r.t. Taylor models.**

**Note that evaluations are built on top of the TaylorSeries arithmetic implementation. This means that the arithmetic will still propagate the lowest polynomial order of all input (substitution) TaylorN values. But, during experimentation evaluation using different TaylorN orders sometimes produced an error in the `power_by_squaring!.power_by_squaring!` function.**

==> Does this have anything to do with TaylorN order [mixtures](https://juliadiff.org/TaylorSeries.jl/stable/userguide/#Mixtures)? Yes! At the time of writing, there is a TODO to implement `power_by_squaring!` for mixtures right next to the definition of the `power_by_squaring!` function definition. Here _mixtures_ refers to a `Taylor1` or `TaylorN` that makes use of variables that have different orders, i.e. a mixture of variable orders.

```julia
# in-place form of power_by_squaring
# this method assumes `y`, `x` and `aux` are of same order
# TODO: add power_by_squaring! method for HomogeneousPolynomial and mixtures
for T in (:Taylor1, :TaylorN)
    @eval function power_by_squaring!(y::$T, x::$T, aux::$T, p::Integer)
      ...
```


### Taylorization

The function `taylor_expand` is provided to generate a Taylor expansion of a given function at a specified expansion point $(c_1, \dots, c_n)$. The function `update!(p, a)` can be called to shift the expansion point of $p$ by $a$. See [the docs](https://juliadiff.org/TaylorSeries.jl/stable/api/#TaylorSeries.taylor_expand).

### Links

1. TaylorSeries, main documentation, [link](https://juliadiff.org/TaylorSeries.jl/stable/)
1. TaylorSeries, `fixorder(a, b)` Github source code, [link](https://github.com/JuliaDiff/TaylorSeries.jl/blob/0298820a6d1f903185e20849c8178ffb0c1cd503/src/auxiliary.jl#L259)
1. TaylorSeries, functions.jl Github source code, [link](https://github.com/JuliaDiff/TaylorSeries.jl/blob/master/src/functions.jl)
1. TaylorSeries, `taylor_expand` docs, [link](https://juliadiff.org/TaylorSeries.jl/stable/api/#TaylorSeries.taylor_expand)


## TaylorModels

### Recommended Changes (taylor models)

- Make use of the function `bound_remainder` to compute the Lagrange remainder explicitly?
- Look at `TaylorModels.shrink_wrapping!` as a reference? It seems to work similarly to preconditioning. Or did we not want to use library functions for this?
- Look at `TaylorModels.rpa` as a reference? It seems to compute the Taylorization of the dynamics. Or did we not want to use library functions for this?
- Look at `TaylorModels.initialize!` as a reference? It seems to compute the Taylorization of the dynamics. Or did we not want to use library functions for this?

### Check the Following

- Does [range bounding](#range-bounding) give a clue on computing the God-given domain?

### Range Bounding

The TaylorModel documentation provides an [example](https://juliaintervals.github.io/TaylorModels.jl/dev/range_bounding/) of range bounding, "_that is to find an interval $I \subseteq \mathbb{R}$ such that $f(x) \in I$ on a given domain $D$._" **Does this give a clue on actually computing the God-given domain instead?**

### fixorder (TMs)

Similar to TaylorSeries, the TaylorModels library defines `fixorder(a::TaylorModelN, b::TaylorModelN)` in `TaylorModels/auxiliary.jl`. This function truncates the polynomials of $a$ and $b$ and adds a bound (over-approximation) of the truncated terms (neglected part) to the corresponding Taylor models' remainders. This bound is computed by calling `TaylorModels.bound_truncation(::Type{TaylorModelN}, a::TaylorN, aux::IntervalBox, order::Int)` which essentially calls the function `TaylorSeries.evaluate` on the `TaylorN` containing only the terms of degree greater or equal to the truncation order. See the source code below.

<details>
<summary>fixorder and bound_truncation source code.</summary>

```julia
function fixorder(a::TaylorModelN, b::TaylorModelN)
    @assert tmdata(a) == tmdata(b)
    get_order(a) == get_order(b) && return a, b

    order = min(get_order(a), get_order(b))
    apol0, bpol0 = polynomial.((a, b))
    apol, bpol = TaylorSeries.fixorder(apol0, bpol0)

    # Bound for the neglected part of the polynomial
    dom = centered_dom(a)
    Δa = bound_truncation(TaylorModelN, apol0, dom, order) + remainder(a)
    Δb = bound_truncation(TaylorModelN, bpol0, dom, order) + remainder(b)

    return TaylorModelN(apol, Δa, expansion_point(a), domain(a)),
        TaylorModelN(bpol, Δb, expansion_point(b), domain(b))
end

function bound_truncation(::Type{TaylorModelN}, a::TaylorN, aux::IntervalBox,
        order::Int)
    order ≥ get_order(a) && return zero(aux[1])
    res = deepcopy(a)
    res[0:order] .= zero(res[0])
    return res(aux)
end
```

</details>

**Since this makes use of the TaylorSeries `evaluate()` function, we assume that the bounds of the truncated terms is computed using Horner's rule.**


### Links

1. TaylorModels, main documentation, [link](https://juliaintervals.github.io/TaylorModels.jl/dev/)
1. TaylorModels, Github repo, [link](https://github.com/JuliaIntervals/TaylorModels.jl)
1. TaylorModels, Juliacon 2018 introductory video of the library, [link](https://www.youtube.com/watch?v=o1h7BUW04NI)
1. TaylorModels, TMJets TM integration Github issue, [link](https://github.com/JuliaReach/ReachabilityAnalysis.jl/issues/588)
1. TaylorModels, `fixorder` source code, [link](https://github.com/JuliaIntervals/TaylorModels.jl/blob/bfc579131c4542b9f10b5455f2b0a26c1729f96f/src/auxiliary.jl#L92)


# Other References

1. JuliaDiff, main website including a list of Julia differentiation packages, "_JuliaDiff is an informal GitHub organization which aims to unify and document packages written in Julia for evaluating derivatives_", [link](https://juliadiff.org/)
1. JuliaDiff, Further reading about _Automatic Differentiation_, [link](https://juliadiff.org/ChainRulesCore.jl/stable/FAQ.html#Where-can-I-learn-more-about-AD-?)
