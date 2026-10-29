# A parent-dependent exact field element that deliberately is not a Number,
# and deliberately has no zero(::Type) or one(::Type) implementation.
struct ParentCoefficient
    value::Rational{BigInt}
    parent::Symbol
end
parentcoefficient(x)=ParentCoefficient(Rational{BigInt}(x),:test_field)
Base.zero(x::ParentCoefficient)=ParentCoefficient(0//big(1),x.parent)
Base.one(x::ParentCoefficient)=ParentCoefficient(1//big(1),x.parent)
Base.real(x::ParentCoefficient)=x
Base.isreal(::ParentCoefficient)=true
Base.isfinite(::ParentCoefficient)=true
Base.iszero(x::ParentCoefficient)=iszero(x.value)
Base.isone(x::ParentCoefficient)=isone(x.value)
Base.:-(x::ParentCoefficient)=ParentCoefficient(-x.value,x.parent)
Base.:^(x::ParentCoefficient,n::Integer)=ParentCoefficient(x.value^n,x.parent)
for op in (:+,:-,:*,:/)
    @eval function Base.$op(a::ParentCoefficient,b::ParentCoefficient)
        a.parent==b.parent || throw(ArgumentError("different parents"))
        ParentCoefficient($op(a.value,b.value),a.parent)
    end
end
Base.:(==)(a::ParentCoefficient,b::ParentCoefficient)=a.parent==b.parent && a.value==b.value
Base.isless(a::ParentCoefficient,b::ParentCoefficient)=isless(a.value,b.value)
Base.:<(a::ParentCoefficient,b::ParentCoefficient)=a.value<b.value
Base.abs(a::ParentCoefficient)=ParentCoefficient(abs(a.value),a.parent)
Base.hash(a::ParentCoefficient,h::UInt)=hash((a.parent,a.value),h)

@testset "Value-based coefficient protocol" begin
    e0,e1,e2,e3=(basisvector(i,parentcoefficient(1)) for i in 0:3)
    o=HyperbolicPoint(e0)
    p=HyperbolicPoint(clvector(parentcoefficient(5),parentcoefficient(1),parentcoefficient(2),parentcoefficient(3)))
    L=Geodesic(e0∧e3)
    H=Horosphere(e0+e3)
    h=HalfTurn(L)
    P=Parabolic(H,e1)
    @test !(eltype(e0)<:Number)
    @test one(e0)==scalar(parentcoefficient(1))
    @test zero(e0)==scalar(parentcoefficient(0))
    @test undual(dual(e0))==e0
    @test parentcoefficient(2)*e0==e0+e0
    @test e0/parentcoefficient(2)==basisvector(0,parentcoefficient(1//2))
    @test h(h(p))==p
    @test inv(P)(P(p))==p
    @test P(H)==H
    @test P(o) in H
    @test rescale(p,parentcoefficient(2))==p
    @test eltype(versor(P))==ParentCoefficient
end
