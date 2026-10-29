using Kleinian, Riley, Poincare, Random
import Farey
f,g = Farey.Monoids.FreeGen("f g")
X,Y = rileygroup(2.0im)
word = Farey.farey_word(2//5,f,g)
isometry = wordisometry(word,Dict(f=>X,g=>Y))
@show word matrix(isometry)
stream = randomlimitset((X,Y);samples=256,rng=Xoshiro(2026))
hull = limitsethull(stream)
@show hull.dimension length(vertices(hull)) length(pleatedges(hull))
# Word evaluation alone does not establish that this word is parabolic or a cusp word.
