using Nemo:next_calkin_wilf, ZZ

next_calkin_wilf(ZZ(321)//113)


CWnemo = accumulate(1:50, init=ZZ(0)//1) do x,y
    next_calkin_wilf(x)
end

CW = CWnemo                     .|> 
        Rational                 |>
        x-> filter(l->l<1, x)   .|>
        l->(farey_neighbours(l)..., l|>ContinuedFraction) .|>
        x->x;

CW[1:4]

cf1 = @cf 1

CW .|> x->println("""
\t$(cf1/x[3])\t\t $(x[1])
\t$(x[2])\t $(cf1/x[4])
""");

QQ = CW[12]


# farey_neighbours(CW[20])

function stern_broca(depth::Int)
    tree = [] 

end

for k in 1:10
Q = ContinuedFraction(0, [2, k, 2])
Q |> println
end

