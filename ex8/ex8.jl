using PyPlot
using LinearAlgebra
pygui(true)
##############
V(x) = (atan(50*(x+2))-atan(50*(x-2)))/3
ħ = 1
q = 5
ϕ(x) = ((4/(pi))^(1/4))*exp(-((x+5)^2)*2)*exp(im*q*x)
le = -40
Ri = 40
N  = 2400
Δ  = abs(le-Ri)/(N)
Cn = Array{ComplexF64,1}(undef,N-1)
x_asix = collect(range(le+Δ,Ri-Δ,length=N-1))

###############
function simpson(Func)
    ans = 0
    for i in 1:N-3
        a = x_asix[i]
        c = x_asix[i+1]
        b = x_asix[i+2]
        o_1 = i
        c_1 = i+1
        b_1 = i+2
        ans += (b-a)/6*( Func(a,o_1)+Func(b,b_1)+4*Func(c,c_1) )
    end # for end 
    return ans
end # end function

function findCn()
    for i in 1:N-1
        Cn[i] =  Vn[i,:]'*ϕ.(x_asix)
    end # end for
end

function time_dependent_persi(t)
    ans = Array{ComplexF64,1}(undef,N-1) # (N-1) list
    for j in 1:N-1 # initialize ans list to 0
        ans[j] = 0.0 + 0.0im
    end
    for i in 1:N-1
        for j in 1:N-1
            ans[j] += Cn[i]*exp(-im*En[i]*t)*Vn[i,j]
        end #end for 
    end # end for
    return ans
end

function H()
    out = fill(0.0im,(N-1,N-1)) # construct (N-1)x(N-1) Matrix
    for i in 1:N-1
        out[i,i] = 1 + (Δ^2)*V(x_asix[i])
        if (i!=1)
            out[i,i-1] = -0.5
        end # end if
        if (i!=N-1)
            out[i,i+1] = -0.5
        end # end if
    end # end for 
    Energy = eigvals(out)/(Δ^2)
    index = sortperm(Energy)
    Vector = eigvecs(out)[:,index]
    return Energy[index], Vector'
end

function main()
    global En, Vn = @time H()
    findCn()
    ψ = V.(x_asix)#./10
    psosos = ϕ.(x_asix)
    grid(true)
    # plot(x_asix,Vn[1,:])
    plot(x_asix,ψ,label="barrier")
    plot(x_asix,abs.(psosos),"x")
    plot(x_asix,abs.(time_dependent_persi(0)),"+")
    legend()
# while true
    # ttt = collect(0:1:500)
    # time_line = fill(0.0im,(length(ttt),N-1))
    # for t in 1:length(ttt)
    #     time_line[t,:] = time_dependent_persi(ttt[t])
    # end
    pause(2)
    println("start")
    # while true
    
####################################################
    # for t in 0:0.05:3
    # cla()
    # title("t = $t")
    # grid(true)
    # ylim(-1.75,1.5)
    # # #xlim(-5,5)
    # this = time_dependent_persi(t)
    # high_point = maximum(abs.(this))
    # plot([le,Ri],[high_point,high_point],"--",label="highest")
    # plot(x_asix,abs.(this),"-",label="papability")
    # plot(x_asix,real.(this).-1.0,"--b",label="real")
    # plot(x_asix,imag.(this).-1.0,"--r",label="imag")
    # legend()
    # ####################################
    # plot(x_asix,ψ,label="barrier")
    # fill_between(x_asix,ψ,0,facecolor="gray")
    # plot(x_asix,psosos,label="ψ")
    # legend()
    # # ##########################################
    # println(t)
    # # plot(x_asix,Vn[2,:])
    # pause(0.05)
    # end # end for
#######################################################
    # end
    println("end")
# end
end

main()