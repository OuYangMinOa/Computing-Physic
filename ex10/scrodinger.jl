using PyPlot
pygui(true)
println("start")
using LinearAlgebra
#######################################################################
global p₀ = 5π
m = 1
β = 1
ħ = 1
τ = m*ħ/(2*β^2)
α = ħ/(2β)
γ(t) = 1+1im*t/τ
norm = (1/(2*m*α^2))^(1/4)
V(x,p=p₀) = step_function(x,p)   #(atan(30*(x-5))-atan(30*(x-9)))/2
ϕ(x) = ((4/(pi))^(1/4))*exp(-((x+5)^2)*2)*exp(im*p₀*x)
ψ(x,t) = norm*1/sqrt(γ(t))*exp(im*p₀*(x - p₀*t/(2m))/ħ)*exp(-((x-p₀*t/m)^2)/(4*α^2*γ(t)))
#########################################################################
global v₀ = 2
function step_function(x)
    x₀ = 3
    λ = 1/sqrt(2*v₀)
    println(λ)
    if (x>=x₀ && x<=x₀+λ)
        return v₀
    else
        return 0
    end # end if 

end 

global cont = true

function handle_close(evt)
    global cont = false
end

function A_mat(Δt,Δx,j,down)
    a = im*Δt/(4*(Δx^2))
    A = fill(complex(0.0,0.0),(j,j))

    #Array{ComplexF64,2}(undef,j,j)

    for i in 1:j
        A[i,i] = 1 + 2a + 2a*(Δx^2)*V(down+(i-1)*Δx)
        if (i>1)
            A[i,i-1] = -a
        end # end if
        if (i<j)
            A[i,i+1] = -a
        end # end if
    end # end for
    return A
end

function B_mat(Δt,Δx,j,down)
    a = im*Δt/(4*(Δx^2))
    B = fill(complex(0.0,0.0),(j,j))
    #Array{ComplexF64,2}(undef,j,j)
    
    for i in 1:j
        B[i,i] = 1 - 2a - 2a*(Δx^2)*V(down+(i-1)*Δx)
        if (i>1)
            B[i,i-1] = a
        end  # end if
        if (i<j)
            B[i,i+1] = a
        end # end if
    end # end for
    return B
end

function loss(y,y♡)
    return log10(sum((y.-y♡).^2)/length(y))
end

function main() 
    up = -40
    down = 40
    N = 500
    #################################################
    Δx = abs(up-down)/N
    x_array = collect(range(up+Δx,down-Δx,length=N-1))
    ################################################
    Δt = 0.05
    T_N = 300
    #################################################
    ψ_zero = Array{ComplexF64,2}(undef,T_N,N-1) # size -> total t X N
    ψ_zero[1,:] = ψ.(x_array,0) #  儲存 t=0 的波函數
    
    A = A_mat(Δt,Δx,N-1,up) # A matrix
    B = B_mat(Δt,Δx,N-1,up) # B matrix

    A_inv_B = (A^-1)*B # A inverse B
    for t in 2:T_N
        ψ_zero[t,:] = A_inv_B * ψ_zero[t-1,:] # ψ(x,t+Δt) = A^-1 ⋅ B ⋅ ψ(x,t) 
    end

    loo = Array{Float64,1}(undef,T_N)
    fig = figure()
    fig.canvas.mpl_connect("close_event",handle_close)
    plot()
    for i =1:5
        cla()
        title("$(5-i)")
        pause(1)
    end #　end for 
    for i in 1:T_N
        cla()
        grid(true)
        xlabel("x")
        ylabel("y")
        title("p = $p₀")#0Δt÷Δx=$(Δt/Δx)")
        ylim(-1.5,1)
        plot(x_array,abs.(ψ_zero[i,:]),label="C.N.")
        plot(x_array,real.(ψ_zero[i,:]).-1,label="real")
        plot(x_array,imag.(ψ_zero[i,:]).-1,label="imag")
        plot(x_array,V.(x_array),"-",label="potential")
        fill_between(x_array,V.(x_array),0,facecolor="gray")
        legend()
        pause(0.01)
        if (!cont)
            break
        end
        loo[i] = loss(abs.(ψ.(x_array,(i-1)*Δt)),abs.(ψ_zero[i,:]))
    end
    subplots()
    title("error")
    xlabel("time")
    ylabel("log(loss)")
    plot(collect(range(1,step = Δt,length=T_N)),loo,"-",label="p₀=$p₀")
    legend()
end

function main2()
    ################################################
    Δt = 0.01
    T_N = 300
    erro = Array{Float64,2}(undef,length(0:0.05:3),T_N)
    for (num,i) = enumerate(0:0.05:3)
        up = -40
        down = 40
        N = trunc(Int64,2*10^i)
        #################################################
        Δx = abs(up-down)/N
        x_array = collect(range(up+Δx,down-Δx,length=N-1))

        #################################################
        ψ_zero = Array{ComplexF64,2}(undef,T_N,N-1) # size -> total t X N
        ψ_zero[1,:] = ψ.(x_array,0) #  儲存 t=0 的波函數
        erro[num,1] = loss(ψ.(x_array,0),ψ_zero[1,:])

        A = A_mat(Δt,Δx,N-1,up) # A matrix
        B = B_mat(Δt,Δx,N-1,up) # B matrix

        A_inv_B = (A^-1)*B # A inverse B
        for t in 2:T_N
            ψ_zero[t,:] = A_inv_B * ψ_zero[t-1,:] # ψ(x,t+Δt) = A^-1 ⋅ B ⋅ ψ(x,t) 
            erro[num,t] = loss(abs.(ψ.(x_array,(t-1)*Δt)),abs.(ψ_zero[t,:]))
        end # end for 
    end # end for 


    fig = figure()
    fig.canvas.mpl_connect("close_event",handle_close)
    plot()
    for i =1:5
        cla()
        title("$(5-i)")
        pause(1)
    end #　end for 

    for (num,i) = enumerate(0:0.05:3)
        cla()
        grid(true)
        xlabel("log(t)")
        ylabel("log(mse)")
        N = trunc(Int64,2*10^i)
        title("N = $(N)")
        plot(log10.(collect(0:T_N-1)*Δt),erro[num,:],"o")
        pause(0.5)
    end # end for 
end # end main2

function mian3()
    arr = collect(0:10)
    totlav = 16
    t_arr = Array{Float64,1}(undef,length(arr))
    err = Array{Float64,2}(undef,length(arr),totlav)
    for (num,i) = enumerate(arr)
        up = -40
        down = 40
        N = 500
        #################################################
        Δx = abs(up-down)/N
        x_array = collect(range(up+Δx,down-Δx,length=N-1))
        ################################################
        Δt = 1.0/(2^i)
        T_N = trunc(Int64,totlav/Δt)+2
        #################################################
        ψ_zero = Array{ComplexF64,2}(undef,T_N,N-1) # size -> total t X N
        ψ_zero[1,:] = ψ.(x_array,0) #  儲存 t=0 的波函數

        t_arr[num] = log10(Δt)
        
        err[num,1] = loss(ψ.(x_array,0),ψ_zero[1,:])
        println(num)
        t = 1
        now_time = 0
        

        A = A_mat(Δt,Δx,N-1,up) # A matrix
        B = B_mat(Δt,Δx,N-1,up) # B matrix

        A_inv_B = (A^-1)*B # A inverse B
        while now_time<totlav
            now_time += Δt
            t += 1
            ψ_zero[t,:] = A_inv_B * ψ_zero[t-1,:] # ψ(x,t+Δt) = A^-1 ⋅ B ⋅ ψ(x,t) 
            if ((now_time)%1==0)
                err[num,trunc(Int64,now_time)] = loss(abs.(ψ.(x_array,(t-1)*Δt)),abs.(ψ_zero[t,:]))
            end # end if
        end

    end

    fig = figure()
    fig.canvas.mpl_connect("close_event",handle_close)
    plot()
    for i =1:5
        cla()
        title("$(5-i)")
        pause(1)
    end #　end for 
    for i in 1:totlav
        cla()
        grid(true)
        xlabel("log(Δt)")
        ylabel("log(mse)")
        title("t = $(i)s")
        plot(t_arr,err[:,i],"o")
        savefig("$i.png")
        pause(0.2)
    end # end for 
end # end main3
# main()

function main4()
    println("Start main 4")
    to_time = 3
    arr = collect(0:0.1:6)
    t_arr = Array{Float64,1}(undef,length(arr))
    err = Array{Float64,1}(undef,length(arr))
    for (num,i) = enumerate(arr)
        up = -20
        down = 20
        N = 500
        #################################################
        Δx = abs(up-down)/N
        x_array = collect(range(up+Δx,down-Δx,length=N-1))
        ################################################
        Δt = 1.0/(2^i)
        #time = to_time
        time = trunc(Int64,to_time/Δt)
        T_N = 300
        #################################################
        ψ_zero = Array{ComplexF64,2}(undef,time+2,N-1) # size -> total t X N
        ψ_zero[1,:] = ψ.(x_array,0) #  儲存 t=0 的波函數
        t_arr[num] = log10(Δt)
        A = A_mat(Δt,Δx,N-1,up) # A matrix
        B = B_mat(Δt,Δx,N-1,up) # B matrix
        A_inv_B = (A^-1)*B # A inverse B
        if (time !=0)
            ψ_zero[2,:] = A_inv_B^time * ψ_zero[1,:] # ψ(x,t+Δt) = A^-1 ⋅ B ⋅ ψ(x,t) 
        else
            ψ_zero[2,:] = ψ_zero[1,:]
        end # end if
        err[num] = loss(abs.(ψ.(x_array,Δt*time)),abs.(ψ_zero[2,:]))
    end # end for
    subplots()
    grid(true)
    xlabel("log(Δt)")
    ylabel("log(mse)")
    title("$to_time second")
    # title("$to_time propagation")
    plot(t_arr,err[:],"o")
end # end function

function main5()
    up = -40
    down = 40
    N = 2000
    Δx = abs(up-down)/N
    #################################################
    Δt = 0.05
    T_N = 300
    #################################################
    arr = collect(0:0.1:10)
    err = Array{ComplexF64,2}(undef,length(arr),T_N) 
    ψ_zero = Array{ComplexF64,3}(undef,T_N,length(arr),N-1)

    x_array = collect(range(up+Δx,down-Δx,length=N-1))
    #################################################
    println("Δt = $Δt , Δx = $Δx")

    for (num,i) in enumerate(arr)
        global p₀ = i

        ################################################
        
        # size -> total t X N
        ψ_zero[1,num,:] = ψ.(x_array,0) #  儲存 t=0 的波函數
        
        err[num,1] = loss(abs.(ψ.(x_array,0)),abs.(ψ_zero[1,num,:]))

        A = A_mat(Δt,Δx,N-1,up) # A matrix
        B = B_mat(Δt,Δx,N-1,up) # B matrix

        A_inv_B = (A^-1)*B # A inverse B

        for t in 2:T_N
            ψ_zero[t,num,:] = A_inv_B * ψ_zero[t-1,num,:] # ψ(x,t+Δt) = A^-1 ⋅ B ⋅ ψ(x,t) 
            err[num,t] = loss(abs.(ψ.(x_array,Δt*(t-1))),abs.(ψ_zero[t,num,:]))
        end

    end # end for 

    plot()
    for i =1:5
        cla()
        title("$(5-i)")
        pause(1)
    end #　end for 

    for yt in 1:T_N
        cla()
        grid(true)
        xlabel("log(p)")
        ylabel("log(mse)")
        title("t = $((yt-1)*Δt)")
        plot(log10.(arr),err[:,yt],"o")
        pause(0.1)
    end # end for

    # cla()
    # for yt in 1:T_N
    #     cla()
    #     grid(true)
    #     xlabel("log(p)")
    #     ylabel("log(mse)")
    #     title("t = $((yt-1)*Δt)")
    #     # title("$to_time propagation")
    #     for (num,i) in enumerate(arr) 
    #         plot(x_array,abs.(ψ_zero[yt,num,:]),"-",label="p = $i")
    #     end # end for 
    #     legend()
    #     pause(0.1)
    # end # end for

end # end function

function main6()
    up = -40
    down = 40
    N = 2000
    #################################################
    Δx = abs(up-down)/N
    x_array = collect(range(up+Δx,down-Δx,length=N-1))
    ################################################
    Δt = 0.05
    T_N = 300
    #################################################
    arr = collect(0:0.1:10)
    err = Array{ComplexF64,2}(undef,length(arr),T_N) 
    ψ_zero = Array{ComplexF64,3}(undef,T_N,length(arr),N-1)

    for (num,i) in enumerate(arr)
        global β = i
        global τ = m*ħ/(2*β^2)
        global α = ħ/(2β)

        ψ_zero[1,num,:] = ψ.(x_array,0) #  儲存 t=0 的波函數

        err[num,1] = loss(abs.(ψ.(x_array,0)),abs.(ψ_zero[1,num,:]))

        A = A_mat(Δt,Δx,N-1,up) # A matrix
        B = B_mat(Δt,Δx,N-1,up) # B matrix

        A_inv_B = (A^-1)*B # A inverse B
        for t in 2:T_N
            ψ_zero[t,num,:] = A_inv_B * ψ_zero[t-1,num,:] # ψ(x,t+Δt) = A^-1 ⋅ B ⋅ ψ(x,t) 
            err[num,t] = loss(abs.(ψ.(x_array,Δt*(t-1))),abs.(ψ_zero[t,num,:]))
        end

    end #end for 

    plot()
    for i =1:5
        cla()
        title("$(5-i)")
        pause(1)
    end #　end for 

    for yt in 1:T_N
        cla()
        grid(true)
        xlabel("log(β)")
        ylabel("log(mse)")
        title("t = $((yt-1)*Δt)")
        plot(log10.(arr),err[:,yt],"o")
        pause(0.1)
    end # end for

    # cla()
    # for yt in 1:T_N
    #     cla()
    #     grid(true)
    #     xlabel("log(β)")
    #     ylabel("log(mse)")
    #     title("t = $((yt-1)*Δt)")
    #     # title("$to_time propagation")
    #     for (num,i) in enumerate(arr) 
    #         plot(x_array,abs.(ψ_zero[yt,num,:]),"-",label="β = $i")
    #     end # end for 
    #     legend()
    #     pause(0.1)
    # end # end for




end # end functionb

#############################################

function H(N,x_asix,Δx)
    out = fill(0.0im,(N-1,N-1)) # construct (N-1)x(N-1) Matrix
    for i in 1:N-1
        out[i,i] = 1 + (Δx^2)*V(x_asix[i])
        if (i!=1)
            out[i,i-1] = -0.5
        end # end if
        if (i!=N-1)
            out[i,i+1] = -0.5
        end # end if
    end # end for 
    Energy = eigvals(out)/(Δx^2)
    index = sortperm(Energy)
    Vector = eigvecs(out)[:,index]
    return Energy[index], Vector'
end

function findCn(N,x_asix,Vn)
    Cn = Array{ComplexF64,1}(undef,N-1)
    for i in 1:N-1
        Cn[i] =  Vn[i,:]'*ψ.(x_asix,0)
    end # end for
    return Cn
end

function time_dependent_persi(t,N,Cn,En,Vn)
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

function cal_tr_re(x_array,arr,r,l)
    r_index = findall(x-> x>r,x_array)[1]
    l_index = findall(x-> x<l,x_array)[end]
    return [arr[r_index:end]'*arr[r_index:end] , arr[1:l_index]'*arr[1:l_index]]
    #   transmiss,  reflection
end


function main7()
    up = -40
    down = 40
    N = 500
    #################################################
    Δx = abs(up-down)/N
    x_array = collect(range(up+Δx,down-Δx,length=N-1))
    ################################################
    Δt = 0.01
    T_N = 300
    #################################################
    
    #############################################
    arr = collect(1:0.1:2)
    hw_8_e = Array{Float64,2}(undef,length(arr),2)
    hw_10_e = Array{Float64,2}(undef,length(arr),2)
    ψ_zero = Array{ComplexF64,3}(undef,length(arr),N-1,2) # size -> total t X N
    #################3
    

    for (num,i) in enumerate(arr)
        
        global p₀ = i

        x₀ = 3
        λ = 1/sqrt(2*v₀)

        ans_t = (5+λ)*m/p₀
        # if ans_t>8
        #     ans_t = 8
        # end # end for 

        prt = trunc(Int64,ans_t/Δt)

        e_v, e_s = H(N,x_array,Δx)
        c = findCn(N,x_array,e_s)


        A = A_mat(Δt,Δx,N-1,up) # A matrix
        B = B_mat(Δt,Δx,N-1,up) # B matrix

        A_inv_B = (A^-1)*B # A inverse B

        ψ_zero[num,:,1] = A_inv_B^prt *ψ.(x_array,0)
        ψ_zero[num,:,2] = time_dependent_persi(ans_t,N,c,e_v,e_s)

        ψ_zero[num,:,1] = ψ_zero[num,:,1]/sqrt(ψ_zero[num,:,1]'*ψ_zero[num,:,1])
        ψ_zero[num,:,2] = ψ_zero[num,:,2]/sqrt(ψ_zero[num,:,2]'*ψ_zero[num,:,2])


        hw_8_e[num,:]  = cal_tr_re(x_array,ψ_zero[num,:,2],5,10)
        hw_10_e[num,:] = cal_tr_re(x_array,ψ_zero[num,:,1],5,10)


    end # end for
    subplots()
    subplot(211)
    grid(true)
    ylabel("%")
    title("HW8")
    plot(arr.^2 ./v₀,hw_8_e[:,1],"-x",label="transmission")
    plot(arr.^2 ./v₀,hw_8_e[:,2],"-x",label="reflection")

    subplot(212)
    grid(true)
    xlabel("p^2")
    ylabel("%")
    title("HW10")
    plot(arr.^2 ./v₀,hw_10_e[:,1],"-x",label="transmission")
    plot(arr.^2 ./v₀,hw_10_e[:,2],"-x",label="reflection")

    savefig("D:\\計算物理\\ex10\\Δx=$Δx,Δt=$Δt.png")
    # plot()
    # for i =1:5
    #     cla()
    #     title("$(5-i)")
    #     pause(1)
    # end #　end for 


    #     subplots()
    #     grid(true)
    #     xlabel("t")
    #     ylabel("log(mse)")
    #     title("MSE V.S. T")
    #     plot(collect(range(1,step = Δt,length = T_N)),hw_8_e,"-",label = "Hw8")
    #     plot(collect(range(1,step = Δt,length = T_N)),hw_10_e,"-",label = "Hw10")
    #     legend()
    # ########################################
    # subplots()
    # plot()
    # for i =1:5
    #     cla()
    #     title("$(5-i)")
    #     pause(1)
    # end #　end for 

    # for (num,i) in enumerate(arr)
    #     cla()
    #     ylim(0,1)
    #     grid(true)
    #     xlabel("x")
    #     ylabel("y")
    #     title("p = $(i)")
    #     plot(x_array,abs.(ψ_zero[num,:,1]),"-",label="HW 10")
    #     plot(x_array,abs.(ψ_zero[num,:,2]),"-",label="HW 8")
    #     plot(x_array,V.(x_array),"-",label="potential")
    #     fill_between(x_array,V.(x_array),0,facecolor="gray")
    #     legend()
    #     savefig("D:\\計算物理\\ex10\\Δx=$Δx,Δt=$Δt,p=$i.png")
    #     pause(0.5)
    # end # end for
end # end function
# @time main()
main7()
function fraw_v()
    x_array = collect(-40:0.1:40)
    grid(true)
    xlabel("x")
    ylabel("y")
    title("potential")
    plot(x_array,V.(x_array),label="potential")
    fill_between(x_array,V.(x_array),0,facecolor="gray")
    legend()
end # end function
# main6()
# main()
# main5()
# fraw_v()
# @time main4()
# mian3()
Ener
