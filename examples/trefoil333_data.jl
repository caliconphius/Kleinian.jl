using Kleinian, Poincare

"""
    trefoil333(; root3=sqrt(3.0), atol=nothing)

The (3,3,3) trefoil domain from the shared Compute Face Pairings construction:
two ideal HMR tetrahedra with their three exterior maps and shared-face gluing.
Pass an exact algebraic sqrt(3), e.g. Nemo QQBar, for exact matrix/shape checks.
Source: https://chatgpt.com/share/6ae1939a-e5a0-83ec-ad4b-cfb285e755f4
"""
function trefoil333(;root3=sqrt(3.0),atol=nothing)
    onevalue,zerovalue = one(root3),zero(root3)
    i = coefficient_imaginary_unit(root3)
    omega = (-onevalue+i*root3)/2
    u,d = (2+omega)/3,(1-omega)/3
    boundary(z) = IdealPoint(1+real(z*conj(z)),2real(z),2imag(z),real(z*conj(z))-1)
    points = (infinity=IdealPoint(onevalue,zerovalue,zerovalue,onevalue),
              zero=boundary(zerovalue),one=boundary(onevalue),u=boundary(u),d=boundary(d))
    plus = HyperbolicTetrahedron(points.infinity,points.zero,points.one,points.u;atol)
    minus = HyperbolicTetrahedron(points.infinity,points.zero,points.one,points.d;atol)
    # Source matrices act on z=x+iy; Kleinian's documented matrix spinor is
    # zeta=x-iy, so conjugate once at this boundary-coordinate interface.
    a = MobiusIsometry(conj.([-1-omega 1+omega; zerovalue omega]))
    b = MobiusIsometry(conj.([-omega omega; zerovalue 1+omega]))
    c = MobiusIsometry(conj.([-onevalue onevalue; omega-2 1-omega]))
    pairings = (FacePairing((2,2),(2,3),a;label="a"),
                FacePairing((1,2),(1,3),b;label="b"),
                FacePairing((2,1),(1,1),c;label="c"),
                FacePairing((1,4),(2,4),one(a);label="shared",internal=true))
    domain = TetrahedralComplex((plus,minus),pairings;atol)
    f,g = c*inv(a),a
    relations = (a^3,b^3,(b*inv(a))^3,c*inv(a)*inv(c)*inv(b)*c)
    (domain=domain,points=points,pairings=(a=a,b=b,c=c),generators=(f=f,g=g),relations=relations)
end
