# TODO: reduce code duplication among the four different dispatches while maintaining
# TODO: performance (i.e. no / little allocation) for each
# TODO:  - Vectors x, y, and dx
# TODO:  - Vectors x, y; scalar dx (uniform mesh)
# TODO:  - Vectors x and dx; function f
# TODO:  - Vector x; scalar dx; function f

# TODO: One option could be to produce an integration weights vector which can just be
# TODO: dotted with the integrand to produce the integral. This has a few benefits:
# TODO:  - Single allocation of a vector which can be reused for multiple integrals, e.g.
# TODO:    for all `q`-point evaluations of a given quantity's Hankel transform
# TODO:  - GPU-friendly because integration becomes a matrix/vector operation
# TODO:  - Centralizes the unique part of each integration method, can reduce code
# TODO:    duplication

abstract type QuadratureMethod end

@doc raw"""
Trapezoidal rule quadrature implemented as described in the
[Wikipedia article](https://en.wikipedia.org/wiki/Trapezoidal_rule).

Non-uniform grid
```math
\int_a^b f(x) dx \approx
\sum_{i=0}^{N} \frac{f(x_{i-1}) + f(x_i)}{2} (x_i - x_{i-1})
```

Uniform grid
```math
\int_a^b f(x) dx \approx
\Delta x \left( \sum_{i=1}^{N-1} f(x_i) + \frac{f(x_N) + f(x_0)}{2} \right)
```
"""
struct Trapezoid <: QuadratureMethod end

@doc raw"""
Composite Simpson's rule quadrature implemented as described in the
[Wikipedia article](https://en.wikipedia.org/wiki/Simpson%27s_rule).

## Non-uniform grid

```math
\int_a^b f(x) dx \approx
\sum_{i=0}{N/2-1} \frac{\Delta x_{2i} + \Delta x_{2i+1}}{6} \left[
    \left( 2 - \frac{\Delta x_{2i+1}}{\Delta x_{2i}} \right) f_{2i} +
    \frac{
        \left( \Delta x_{2i} + \Delta x_{2i+1} \right)^2
        }{
        \Delta x_{2i} \Delta x_{2i+1}
    } f_{2i+1} +
    \left( 2 - \frac{\Delta x_{2i}}{\Delta x_{2i + 1}} \right) f_{2i + 2}
\right]
```
where
```math
f_{k} = f \left( a + \sum_{i=0}^{k-1} \Delta x_{i} \right)
```

In the case of an odd number of subintervals, the above formulae are used up to the second
to last interval, and the last interval is computed seperately as
```math
\frac{
    2 \Delta x_{N-1}^2 + 3 \Delta x_{N-1} \Delta x_{N-2}
    }{
    6\left( \Delta x_{N-1} + \Delta x{N-1} \right)
} f_{N} +
\frac{
    \Delta x_{N-1}^2 + 3 \Delta x_{N-1} \Delta x_{N-1}
    }{
    6 \Delta x_{N-1}
} f_{N-1} +
\frac{
    \Delta x_{N-1}^3
    }{
    6 \Delta x_{N-2} \left( \Delta x_{N-2} + \Delta x_{N-1} \right)
} f_{N-2}
```

## Uniform grid

```math
\int_a^b f(x) dx \approx
\frac{1}{3} \Delta x \left[
    f(x_0) + 4 \sum_{i=0}{N/2} f(x_{2i-1}) + 2 \sum_{i=0}{N/2-1} f(x_{2i}) + f(x_N)
\right]
```

In the case of an odd number of subintervals, the above formula is used for the first N-1
subintervals, and the last subinterval is handled with the Trapezoid method. This approach
is generally more well behaved in the specific case of pseudopotential quantities because
the function is generally close to zero in the last interval and the error made by
the Trapezoid rule w.r.t. Simpson's rule has a smaller effect on the value of the integral.
This approach is also equivalent to the approach taken in the non-uniform case.
"""
struct Simpson <: QuadratureMethod end

@doc raw"""
QuantumESPRESSO Simpson's (1/3) rule quadrature. The expression is equivalent to the
uniform grid case of `Simpson` for an even number of subintervals. However, the
derivative of the mesh `dx` is used in place of the finite difference between adjacent
mesh points `Δx`.
"""
struct QESimpson <: QuadratureMethod end

@doc raw"""
ABINIT quadrature -- a collection of composite closed Newton-Cotes rules. Behavior depends
on the number of available points. This method _only_ supports quantities on uniform
meshes!

```math
\int_a^b f(x) dx \approx \begin{cases}
    \frac{\Delta x}{72} \left[
        23.75 f_1 + 95.10 f_2 + 55.20 f_3 + 79.30 f_4 + 70.65 f_5 +
        72 \sum_{i=6}^{N-5} f_i +
        70.65 f_{N-4} 79.30 f_{N-3} + 55.20 f_{N-2} + 95.10 f_{N-1} + 23.75 f_{N}
    \right] & N >= 10 \\
    \frac{\Delta x}{48} \left[
        17 f_1 + 59 f_2 + 43 f_3 + 49 f_4 +
        48 f_5 +
        49 f_6 + 43 f_7 + 59 f_8 + 17 f_9
    \right] & N = 9 \\
    \frac{\Delta x}{48} \left[
        17 f_1 + 59 f_2 + 43 f_3 + 49 f_4 +
        49 f_5 + 43 f_6 + 59 f_7 + 17 f_8
    \right] & N = 8 \\
    \frac{\Delta x}{48} \left[
        17 f_1 + 59 f_2 + 43 f_3 + 50 f_4 + 43 f_5 + 59 f_6 + 17 f_7
    \right] & N = 7 \\
    \frac{\Delta x}{48} \left[
        17 f_1 + 59 f_2 + 44 f_3 + 44 f_4 + 59 f_5 + 17 f_6
    \right] & N = 6 \\
    \frac{\Delta x}{3} \left[
        f_1 + 4 f_2 + 2 f_3 + 4 f_4 + f_5
    \right] & N = 5 \\
    \frac{\Delta x}{8} \left[
        3 f_1 + 9 f_2 + 9 f_3 + 3 f_4
    \right] & N = 4 \\
    \frac{\Delta x}{3} \left[
        1 f_1 + 8 f_2 + 1 f_3
    \right] & N = 3 \\
    \frac{\Delta x}{2} \left[
        f_1 + f_2
    \right] & N = 2 \\
    \frac{\Delta x} \left[
        f_1
    \right] & N = 1
\end{cases}
```
"""
struct AbinitCorrectedTrapezoid <: QuadratureMethod end
