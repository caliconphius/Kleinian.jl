using GLMakie
xx = BigFloat[1,1]
A = cf_matrix(3//5; z=BigFloat(0))|>Array
A = [
  1 2; 0 1
  ]
matx(s,t) = ([s,s*t]')*(A) * Float64[1,0]
maty(s,t) = ([s,s*t]')*(A) * Float64[0,1]

tors = x-> x - floor(x)

ss = range(0,1, length=100);
ts = ss./2;
Xs = [matx(s,t) for s in ss, t in ts];
Ys = [maty(s,t) for s in ss, t in ts];
Zs = [t for s in ss, t in ts];
fig = Figure()
Axis3(fig[1, 1], aspect = :equal)
GLMakie.surface!((Xs),(Ys),Zs)
fig

# # GLMakie.activate!(; screen_config...)
# # function mandelbrot(x, y)
# #     z = c = x + y*im
# #     for i in 1:30.0; abs(z) > 2 && return i; z = z^2 + c; end; 0
# # end
# # ,sin(π*ζ)*exp(2*π*im*ϕ)
# function sphere3(ζ,θ)
#   (cos(π*ζ)* exp(2*π*im*θ))
# end
# xs = range(0, 1, length=100)
# ys = range(-1,1, length=100)
# zs = range(-1,1, length=100)
# ws = [sphere3(x,y) for x in xs, y in ys]
# z1s = first.(ws)
# fig = Figure();
# ax1 = Axis3(fig[1,1], inte)

# # CairoMakie.surface!(ax1, real.(ws),imag.(ws),zeros(size(ws)),color=abs.(ws));
# fig


