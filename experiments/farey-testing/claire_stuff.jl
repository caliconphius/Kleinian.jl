Qb = algebraic_closure(QQ)
Ka, (x,y,z,a,b) = Qb[:x,:y, :z, :a, :b]

DD = (1-x^2 - y^2 - z^2 + a*y*y + b*z*z)
DD - 2x^2
dx1 = Nemo.derivative(DD*x,x)
dx1 = Nemo.derivative(DD(y,x,z,a,b)*y,y)
DD(x,1+x,x,a,b)

(DD(1+x,0,0,a,b) - DD(1,0,0,a,b))//x
(DD(1,x,y,a,b) - DD(1,0,0,a,b))//x

DD(x,x,1+x,a,b)
DD(x,y,1+z,a,b)
dx1(1,0,0,a,b)
DD(0,0,1,a,b)
Nemo.derivative(DD*x,y)
2*x*y*(a-1)
2*x*z*(b-1)
