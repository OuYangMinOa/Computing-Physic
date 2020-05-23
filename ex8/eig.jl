using LinearAlgebra
using PyPlot
pygui(true)

function V(x)
    if ((x<=-5) & (x>=5))
        return 0
    elseif ((x>=-5) & (x<=5))
        return -0.5
    else
        return 0
    end # end if 
end
function main()
    N = 500
    up = 10
    down = -10
    Δ = (abs(up-down)/(N))
    x = collect(range(down+Δ,up-Δ,length = N-1))
    A = fill(0.0,(N-1,N-1))
    #V(x) = (atan(30*(x+3))-atan(30*(x-3)))/20
    for i in 1:N-1
        A[i,i] = 1+(Δ^2)*V(x[i])
        if (i!=1)
            A[i,i-1] = -0.5
        end
        if (i!=N-1)
            A[i,i+1] = -0.5
        end # end if
    end # end for
    ge_3 = eigvals(A)
    idx_3 = sortperm(ge_3)
    ge_3 = ge_3[idx_3]./(Δ^2)
    #println(ge_3[1:10])
    vec = (eigvecs(A)[:,idx_3])'
    # for i in 1:10
    #     cla()
    #     plot(x,abs.(vec[i,:]))
    #     pause(0.5)
    # end
    subplots()
    # for i in 1:5
    # cla()
    plot(x,V.(x))
    grid(true)
    nn = 1
    print(ge_3[nn])
    plot(x,(vec[nn,:]),"--",label="E$nn")
    # pause(0.5)
    # end
    # println(sum(vec[1,:].^2))
    # plot(x,(vec[2,:]),"--")
    # plot(x,vec[1,:],"o")
end


@time main()