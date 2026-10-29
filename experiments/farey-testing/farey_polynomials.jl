include("sl2algebra.jl")
# Kx, x = polynomial_ring(QQ, "x")
Kx, xz = polynomial_ring(QQ, :z)
using GLMakie
detA(eQ(1//2))
f, g = FreeGen("f g")
function Φq(Q; tf = 2, tg = 2, tfg = xz+2)
  tr(eQ(Q))(tf,tg,tfg)
end

Φq(5//12)
Φq(1//3)(1)
Q = 4//11

Fp = Φq(Q)
C = 4
Cf = AcbField()
xs= range(-C, C, length=100)
ys= range(-C, C, length=100)
Fs = [Fp(Cf(x + im*y)) for x in xs, y in ys];
Fas = abs.(Fs).|>BigFloat.|>float;
Frs = real.(Fs).|>BigFloat.|>float;
Fis = imag.(Fs).|>BigFloat.|>float;
Farg = angle.(Fs).|>BigFloat.|>float;
extrema(Fas)
fig = Figure()
ax = Axis(fig[1, 1])
contour(xs, ys, (Fis))
# levels = Makie.get_tickvalues(Makie.LinearTicks(30), extrema(Fs)...);
# levels = Makie.get_tickvalues(Makie.LogTicks(7), extrema(Fs)...);
surface(xs, ys, Frs./(1 .+Fas), colormap=:viridis,color= Fis./(1 .+Fas))
surface(xs, ys, Fas./1e8, colormap=:viridis)
contour!(xs, ys, (Fis)/1e8, levels=-1:.5:1,color=:red)
r =C
θr = 2pi/7
lines!([-r*cos(θr),r*cos(θr)],[-r*sin(θr),r*sin(θr)])
x = -2:0.5:2
y = -2:0.5:2
# The curvilinear grid:
xs = [x + 0.01y^3 for x in x, y in y]
ys = [y + 10cos(x/40) for x in x, y in y]
Fs
# Now, for simplicity, we calculate the `zs` values to be
# the radius from the center of the grid (0, 10).
Fs = sqrt.(xs .^ 2 .+ (ys .- 10) .^ 2)

# We can use Makie's tick finders to get some nice looking contour levels:

# and now, we plot!
# fig, ax, srf = surface(xs, ys, fill((0f0), size(Fs)); color=Fs, shading = NoShading, axis = (; type = Axis, aspect = DataAspect()))

# ctr = contour!(ax, xs, ys, Fs; color = :orange, levels = levels, labels = true, labelfont = :bold, labelsize = 12)

fig

