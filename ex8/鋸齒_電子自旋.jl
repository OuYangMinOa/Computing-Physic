using PyPlot
using LinearAlgebra
using Printf
using Dates
using FFTW
pygui(true)
println(repeat("#",50))
############################################
# ZigZag
# 基本函式
# 製造 down + 1~ up 的 亂數矩陣 (長度 num)
function make_rand(num,up,down = 0)
    flag = true
    while(flag)
        global arr
        arr = rand(num)
        arr = sort(map(ceil,(arr .* (up-down) )))
        arr = map(x-> x + down,arr)
        flag_2 = true
        for i in arr
            now_num = i
            if (count(x->x==i,arr) > 1)
                flag_2 = false  
            end # end if
        end # end for
        if (flag_2)
            flag = false
        end # end if
    end # end while
    return arr
end # end function
function click_state(event)
    if (event.button==2)
        println("plot: ky = $(event[:xdata])")
        draw_ZigZag(event.xdata)
    end # end if
end # end function

###########################################
global plot_3d = false
global plot_state = false
global plot_emergy = true
global change_first = true
global plot_scar = false
################## ZigZag ##################
global scar_arr = [ -0.1,0.1    ]
#global scar_arr = [ pi-0.1,pi+0.1]
global ts_1_2 = complex(0,0.08)  # 旋轉  
global t_1_2 = 1 # 躍遷
global a = 1
global graphene_length = 50
global n_3 = graphene_length * 2 #　2:the spin up spin down
global small_3 = n_3 ÷ 2 -1
global rasba_const = complex(0,0.5)
global Magnetic_const = 0.2
global v0 = 0
##########################################
function generate(n)
    gg = []
    if ((n%4 == 2)|(n%4 == 3))
        for i in 2:4:n
            push!(gg,i)
            if (i < n)
            push!(gg,i+1)
        end
        end
    elseif ((n%4 == 1)|(n%4 == 0))
        push!(gg,1)
        for i in 4:4:n
            push!(gg,i)
            if (i<n)
            push!(gg,i+1)
        end
        end
    end
    return gg
end
function condition_2_1_2(i,j,k) # ZigZag型躍遷
        #################### 電子自旋  #################
        if ((abs(i-j)%2)==0) # 判斷自旋方向是否同
            i = (i+1)÷2  # 改變 i,j
            j = (j+1)÷2
        else
            return 0
        end # end if
        ###############################################
        t_1_2 = 1
        # if ((i,j) in [(1,2),(2,1),(n_3,n_3-1),(n_3-1,n_3)] )
        #     t_1_2 = 0.5t_1_2
        # end
        out = 0
        if ((i%4==1))
            if (i-j ==-1)
            out = -t_1_2 - t_1_2*exp(complex(0,1)*k*a)
            elseif (i-j == 1)
            out = -t_1_2
            end
        elseif ((i%4==0))
            if (i-j ==1)
            out = -t_1_2 - t_1_2*exp(complex(0,1)*k*a)
            elseif (i-j == -1)
            out = -t_1_2
            end
        elseif ((i%4==2) )
            if (i-j ==1)
            out = -t_1_2 - t_1_2*exp(-complex(0,1)*k*a)
            elseif (i-j == -1)
            out = -t_1_2
            end
        elseif ((i%4==3))
            if (i-j ==-1)
            out = -t_1_2 - t_1_2*exp(-complex(0,1)*k*a)
            elseif (i-j == 1)
            out = -t_1_2
            end
        end
        return out # + potential(i,j)
end
function condition_1_2(i,j,k) # ZigZag型旋轉
        #################### 電子自旋  #################
        if ((abs(i-j)%2)==0) # 判斷自旋方向是否同
            i = (i+1)÷2  # 改變 i,j
            j = (j+1)÷2
        else
            return 0
        end # end if
        ###############################################
        out = 0
        if ((i%4==1))
            if (i-j ==2)
                out = ts_1_2-ts_1_2*exp(complex(0,1)*k*a)
            elseif (i-j == 0)
                out = ts_1_2*exp(complex(0,1)*k*a)-ts_1_2*exp(-complex(0,1)*k*a)
            elseif (i-j == -2)
                out = ts_1_2 - ts_1_2*exp(complex(0,1)*k*a)
            end
        elseif ((i%4==2))
            if (i-j == 2)
                out = ts_1_2-ts_1_2*exp(-complex(0,1)*k*a)
            elseif (i-j == 0)
                out = ts_1_2*exp(-complex(0,1)*k*a)-ts_1_2*exp(complex(0,1)*k*a)
            elseif (i-j == -2)
                out = ts_1_2 - ts_1_2*exp(-complex(0,1)*k*a)
            end
        elseif ((i%4==3))
            if (i-j == 2)
                out = -ts_1_2+ts_1_2*exp(-complex(0,1)*k*a)
            elseif (i-j == 0)
                out = ts_1_2*exp(complex(0,1)*k*a)-ts_1_2*exp(-complex(0,1)*k*a)
            elseif (i-j == -2)
                out = -ts_1_2 + ts_1_2*exp(-complex(0,1)*k*a)
            end
        elseif ((i%4==0))
            if (i-j == 2)
                out = -ts_1_2+ts_1_2*exp(complex(0,1)*k*a)
            elseif (i-j == 0)
                out = -ts_1_2*exp(complex(0,1)*k*a)+ts_1_2*exp(-complex(0,1)*k*a)
            elseif (i-j == -2)
                out = -ts_1_2 + ts_1_2*exp(complex(0,1)*k*a)
            end
        end
        if ((i==2)&(j==1))
            #println(out)
        end
        return out
end
function rashba_d1(i,j) # d1
    if ((i-j)%2==0)
        return 0
    else
        return 1
    end
end # end function

function rashba_d2(i,j) # d2
    if ((i-j)%2==0)
        return 0
    elseif ((i%2==1) | (j%2==0)) # up down
        return complex(0.5,-sqrt(3)/2)
    elseif ((i%2==0) | (j%2==1)) # down up
        return complex(-0.5,sqrt(3)/2)
    end
end # end function

function rashba_d3(i,j) # d3
    if ((i-j)%2==0)
        return 0
    elseif ((i%2==1) | (j%2==0)) # up down
        return complex(-0.5,-sqrt(3)/2)
    elseif ((i%2==0) | (j%2==1)) # down up
        return complex(-0.5,sqrt(3)/2)
    end
end # end functionW

function rashba_d4(i,j) # d4
    if ((i-j)%2==0)
        return 0
    elseif ((i%2==1) | (j%2==0)) # up down
        return complex(0,-1)
    elseif ((i%2==0) | (j%2==1)) # down up
        return complex(0, 1)
    end
end  # end function

function Rashba(i,j,k)
    original_i = (i+1)÷2
    original_j = (j+1)÷2
    if ((original_i%4 == 1))
        if (original_i - original_j == -1)
            return rashba_d2(i,j)*exp(complex(0,1)*k*a) + rashba_d3(i,j)
        elseif (original_i - original_j == 1)
            return -rashba_d4(i,j)
        end # end if

    elseif ((original_i%4 == 2))
        if (original_i - original_j == -1)
            return rashba_d4(i,j)
        elseif (original_i - original_j == 1)
            return -rashba_d3(i,j) -  rashba_d2(i,j)*exp(-complex(0,1)*k*a)
        end # end if
    elseif ((original_i%4 == 3))
        if (original_i - original_j == -1)
            return rashba_d2(i,j) + rashba_d3(i,j)*exp(-complex(0,1)*k*a)
        elseif (original_i - original_j == 1)
            return -rashba_d4(i,j)
        end # end if
    elseif ((original_i%4 == 0))
        if (original_i - original_j == 1)
            return -rashba_d3(i,j)*exp(-complex(0,1)*k*a) - rashba_d2(i,j)
        elseif (original_i - original_j == -1)
            return rashba_d4(i,j)
        end # end if
    end # end if
    return 0

end  #end function

function Magnetic(i,j,k)
    if ((abs(i-j)%2)==0) # 自旋方向相同
        if (i==j)  # 改變 i,j
            if (i%2==1)
             return Magnetic_const
            else
                return -Magnetic_const
            end
        end # end if
    else # 自旋不同
        return 0
    end # end if
    return 0
end

function potential(i,j,k)
    if (i==j)
        if (i in (1,2,3,4,n_3,n_3-1,n_3-2,n_3-3))
            # println(i)
            return v0
        else
            return v0
        end
    else
        return 0
    end
end
#############################################
function draw_ZigZag(kk = "0")
        global energy = [] # 能帶圖矩陣
        println("ZigZag n = $n_3")
        len_length = 1000 # 切割

        if (kk != "0")
                global plot_emergy = false
                global plot_state = true
                time_array = kk
                println("plot state ky = $kk")
            else
                global plot_emergy = true
                time_array = collect(range(-3π,3π,length = len_length))#-3π,3π
                println("plot dispersion n = $n_3")
                # time_array
        end # end if

        nearest = time_array[argmin(map(abs, time_array))] # 最近的

        for (num,k) in enumerate(time_array)
            H_3 = fill(complex(1.0,1.0),  (n_3,n_3))
            for i in 1:n_3
                for j in 1:n_3
                    out_3 =  condition_2_1_2(i,j,k) + condition_1_2(i,j,k) +potential(i,j,k) #+ test_v0(i,j,k) + Vm(i,j) + sym_zero_one(i,j,k)
                    #+ potential(i,j)
                    #out_3 += rasba_const*Rashba(i,j,k)
                    out_3 += Magnetic(i,j,k)
                    H_3[i,j] = out_3
                end # end for
            end # end for
            ge_3 = eigvals(H_3)
            ge_3 = map(x-> real(x),ge_3)
            idx_3 = sortperm(ge_3)
            ge_3 = ge_3[idx_3]
            push!(energy,ge_3)
            # plot 邊緣 state
            if (plot_state &(k==nearest))
                vec = rotl90(map(abs,eigvecs(H_3)[:,idx_3]))#[end:-1:1,end:-1:1]
                ge_vec = generate(n_3)#collect(1:2:2*n_3)
                new_vec = vec[small_3,:]
                # 畫
                subplots()
                title("ZigZag -3π ~ 3π n = $small_3 k = $kk")
                plot(new_vec[ge_vec],"--")
                plot(new_vec[ge_vec],"o")
            end # end if
        end # end for

        if (plot_emergy)
            #subplots()
            title("ZigZag -3π ~ 3π ")#v₀ = $v_0 vₐ = $Vₐ vₘ = $vmid")
            fig = figure()
            fig.canvas.mpl_connect("button_press_event", click_state)
            positive_e = collect(1:2:n_3)
            negetive_e = collect(2:2:n_3)
            println("energy = $small_3")
            plot(time_array,hcat(energy...)[positive_e,:]',"b--")
            plot(time_array,hcat(energy...)[negetive_e,:]',"r--")
            grid(true)
            #plot(time_array,hcat(energy...)[small_3,:],"r--")
            ################### plot_scar #####################################
            if (plot_scar)
                    
                for i in scar_arr
                    get_min = findmin(map(abs,(map(x -> x-i,time_array))))
                    this_pos = hcat(energy...)[small_3,:][get_min[2]]
                    println(get_min,this_pos)
                scatter(i,this_pos,alpha=1,marker = "o",s = 60,c="black",label="($i,$this_pos)",zorder=3)
                # if (i==-0.1)
                #     annotate("($i,$this_pos)",xy=(i-3.5,this_pos),textcoords="offset points")
                # else
                #     annotate("($i,$this_pos)",xy=(i+0.5,this_pos),textcoords="offset points")
                # end
                end # for
            end # end if
            ########################################################################
        end # end if
end # end function
draw_ZigZag()
