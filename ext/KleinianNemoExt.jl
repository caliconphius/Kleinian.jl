module KleinianNemoExt

using Kleinian
import Nemo

# QQBar has exact real/complex algebraic elements, but no Base.isfinite in
# Nemo 0.56. It contains no infinities. No methods on foreign functions/types
# are added here: only Kleinian's coefficient protocol is extended.
Kleinian.coefficient_isfinite(::Nemo.QQBarFieldElem) = true
Kleinian.coefficient_imaginary_unit(x::Nemo.QQBarFieldElem) = sqrt(-one(x))

end
