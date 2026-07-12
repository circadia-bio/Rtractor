using Random
using DelimitedFiles

rng                = MersenneTwister()
Random.seed!()

StatsBlock         = 3

##########################################################
############Time, Window and Sampling#####################
##########################################################

Window_Size        = 300

samples            = 18000

x_rand             = zeros(Int64,(samples))
y_rand             = zeros(Int64,(samples))

for count=1:samples
    x_rand[count]  = round(Int64,(rand(rng)*(Window_Size-StatsBlock)))
    y_rand[count]  = round(Int64,(rand(rng)*(Window_Size-StatsBlock)))
end

##########################################################
#####################Arguments S_Max######################
##########################################################

Frac               = 10
Frac2              = 10
Max_Eps            = 0.499999999
Low_Eps            = 0.000000001
Max_Micro          = Int64(2^(StatsBlock*StatsBlock))
Stats              = zeros(Float64,Max_Micro)
Stats_Max          = zeros(Float64,Max_Micro)

pow_vec           = zeros(Int64,(StatsBlock*StatsBlock))
for i=1:(StatsBlock*StatsBlock)
    pow_vec[i]=Int64(2^(i-1))
end

##########################################################
##################Entropy Function########################
##########################################################

function Max_Entropy(Serie)

    S_Max=0.0; Threshold_Max=0.0 ; Threshold=Low_Eps; Var_Eps = (Max_Eps-Low_Eps)/Frac
    
    for i=1:Frac2
        if (i > 1)
            Threshold=Threshold_Max-Var_Eps
            Var_Eps=2*Var_Eps/Frac
        end
        for j=1:Frac
            Stats[:].=0
            for count=1:samples
                Add=0
                for count_y=1:StatsBlock
                    for count_x=1:StatsBlock
                        if (abs(Serie[x_rand[count]+count_x]-Serie[y_rand[count]+count_y]) <= Threshold)
                            a_binary=1
                        else
                            a_binary=0
                        end
                        Add=Add+a_binary*pow_vec[count_x+((count_y-1)*StatsBlock)]
                    end
                end
                Stats[Add+1]+=1
            end

            S=0
            for k=1:Max_Micro
                if (Stats[k] > 0)
                    S+=(-(Stats[k]/(1.0*samples))*(log((Stats[k]/(1.0*samples)))))
                end
            end
            if (S > S_Max)
                S_Max          = S
                Threshold_Max  = Threshold
                Stats_Max[:]   = (Stats[:]./samples)
            end
            Threshold=Threshold+Var_Eps
        end
    end

    return Stats_Max,S_Max,Threshold_Max
        
end

##########################################################
#####################Main Function########################
##########################################################

function main()

    all_cols           = 9
    Amp                = 10
    
    V_Mean_Window      = zeros(Float64,Window_Size)
    MicroStates        = zeros(Float64,Max_Micro)

    
    
    VEC_IND  =["S01","S02","S04","S05","S07","S08","S09","S10","S11","S14","S15","S16","S17","S18","S19","S20","S22","S23","S24","S25","S26","S27","S28","S29","S30"]

    for TERM=1:25

        String_Archive1   = VEC_IND[TERM]
        String_Archive2   = String_Archive1*"_POA"
        Serie_Extra       = readdlm(String_Archive1*"/"*String_Archive2*".txt") #POA #POF #ROA #ROF
    
        r_Ex              = 6
        
        a_Ex              = size(Serie_Extra,1)
        a_s               = a_Ex-r_Ex
        Serie_r           = zeros(a_s,10)
        Serie_r[:,:]      .= Serie_Extra[(r_Ex+1):a_Ex,1:10]
    
        t_end             = a_s
    
        for n_cols=1:all_cols
        
            Component                = n_cols+1
            max_loops2               = floor(Int64,(a_s/(Amp*Window_Size)))-1

            OutPut_Vec               = zeros(Float64,max_loops2,3)
        
            for loop_int=1:max_loops2
            
                for a_count=1:Window_Size
                    V_Mean_Window[a_count]  = Serie_r[floor(Int64,(((loop_int-1)*Amp*Window_Size)+a_count*Amp)),Component]
                end
            
                Maximum_Value,Local_Number  = findmax(V_Mean_Window)
                Minimum_Value,Local_Number  = findmin(V_Mean_Window)
                if ((Maximum_Value-Minimum_Value) != 0.0)
                    V_Mean_Window.=((V_Mean_Window.-Minimum_Value)./(Maximum_Value-Minimum_Value))
                end

                OutPut_Vec[loop_int,1]       = (loop_int*Amp*Window_Size)*0.001
                
                MicroStates[:],OutPut_Vec[loop_int,2],OutPut_Vec[loop_int,3] = Max_Entropy(V_Mean_Window)
            
            end

            writedlm("Quantifiers_$(n_cols)_col_"*String_Archive2*".dat",OutPut_Vec)
            println("S-Eps Calc:"," ",n_cols,"/9")

        end
        
    end

end

main()
