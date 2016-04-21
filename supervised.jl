using DecisionTree

"""
Crossvalidation

WIP
"""
function cross_validate(features::Vector{Matrix{Float64}},
    labels::Vector{Vector{Int}}; method::Symbol=:leave_one_out)
    for i=1:length(features)
        println("evaluating $i")
        F = nothing
        L = nothing
        if(i == 1)
            F = features[2]
            L = labels[2][1:size(F,2)]
            for j=3:length(features)
                F = hcat(F, features[j])
                L = vcat(L, labels[j][1:size(features[j],2)])
            end
        else
            F = features[1]
            L = labels[1][1:size(F,2)]
            for j=2:length(features)
                if(i!=j)
                    F = hcat(F, features[j])
                    L = vcat(L, labels[j][1:size(features[j],2)])
                end
            end
        end
        #println(size(F))
        #println(size(L))
        #println(size(features[i]))
        #println(size(labels[i]))

        validate(F, L, features[i], labels[i][1:size(features[i],2)], trees=100,
        feats=2,
        doPlot=i)
    end
end

"""
validate feature representation
"""
function validate(train::Matrix{Float64}, train_label::Vector{Int},
                  test::Matrix{Float64}, test_label::Vector{Int}; trees::Int=4,
                  feats::Int=20,doPlot=nothing)
    model = build_forest(train_label, train', feats, trees)#, 200, 20, 0.5)
    out = apply_forest(model, test')
    println("classification accuracy = ", mean(out.==test_label))
    if(doPlot != nothing)
        figure(doPlot)
        PyPlot.clf();
        plot(test_label)
        plot(out+0.1)
    end
    out
end


"""
Create a confusion matrix

Arguments:
 - A/B    : Input sequences
 - values : States to build confusion matrix for

All values outside of 'values' are discarded
"""
function make_conf(A::Vector{Int},B::Vector{Int},values::Vector{Int})
    N    = length(values)
    conf = zeros(N,N)
    rm   = Dict{Any, Int}() #Remapper
    for i=1:N
        rm[values[i]] = i
    end

    M = length(A)
    for j=1:M
        a = A[j]
        b = B[j]
        if(a in values && b in values)
            conf[rm[a], rm[b]] += 1
        end
    end
    conf
end
