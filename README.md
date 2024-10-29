LambertW for Java
=================

Java implementation of the Lambert W(x) function.

The lambert W function is defined as the solution for: $x \cdot e^{x} = z$ , then $W(z) = x$

$$
W(x \cdot e^{x}) = x
$$

The implementation is based on the FORTRAN algorithm by Toshio Fukushima, and the [Boost C++ libraries](https://www.boost.org/doc/libs/develop/libs/math/doc/html/math_toolkit/lambert_w.html)

## Usage
In the interval $\left[0, \infty\right)$ there is a single real solution that can be computed by

```Java
    double w = Lambert.lambertW0(z);
```

In the interval $\left[ -\frac{1}{e},0 \right)$ there are two real solutions that can be computed by

```Java
    double w1 = Lambert.lambertW0(z);
    double w2 = Lambert.lambertWm1(z);
```

## Reference
Toshio Fukushima, “Precise and fast computation of Lambert W-functions without transcendental function evaluations”\
May 2013, Journal of Computational and Applied Mathematics 244:77-83\
[DOI:10.1016/j.cam.2012.11.021](https://doi.org/10.1016/j.cam.2012.11.021)

