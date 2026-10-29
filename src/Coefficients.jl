"""
    coefficient_isfinite(x)

Backend hook for finiteness validation. Defaults to isfinite(x). Exact fields
without infinities may specialize this to true; uncertain enclosures must not
be reported as certified finite/invertible by an adapter.
"""
coefficient_isfinite(x) = isfinite(x)

"""
    coefficient_imaginary_unit(x)

Return an imaginary unit in a complex coefficient domain compatible with x.
Backends with parent objects can specialize this without subtyping Number.
"""
coefficient_imaginary_unit(x) = complex(zero(real(x)),one(real(x)))

"""Construct x+i*y in a compatible coefficient domain, using the backend hook."""
coefficient_complex(x,y) = x+coefficient_imaginary_unit(x)*y
