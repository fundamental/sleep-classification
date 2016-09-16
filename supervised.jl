using DecisionTree

"""
Crossvalidation

WIP
"""
function cross_validate(features::Vector{Matrix{Float64}},
    labels::Vector{Vector{Int}};
    method::Symbol=:leave_one_out, 
    selection::Vector{Int}=collect(1:length(features)),
    feats::Int=20,
    trees::Int=40)
    out = nothing
    for i=selection#[50:end]
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

        out = validate(F, L, features[i], labels[i][1:size(features[i],2)],
        trees=trees, feats=feats, doPlot=i)
    end
    out
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
Validate across datasets
"""
function validate_across(train::Vector{Matrix{Float64}},
                         train_label::Vector{Vector{Int}},
                         test::Vector{Matrix{Float64}},
                         test_label::Vector{Vector{Int}};
                         trees::Int=40,
                         feats::Int=3,
                         doPlot=nothing)
    F = train[1]
    L = train_label[1]
    N = length(train)
    M = length(test)
    for i=2:N
        F = hcat(F, train[i])
        L = vcat(L, train_label[i])
    end
    model = build_forest(L, F', feats, trees)#, 200, 20, 0.5)
    for i=1:M
        out = apply_forest(model, test[i]')
        #println("classification accuracy = ", mean(out.==test_label[i]))
        println(i, " ", mean(out.==test_label[i]))
        if(doPlot != nothing)
            figure(doPlot)
            PyPlot.clf();
            plot(test_label[i])
            plot(out+0.1)
        end
        out
    end
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
