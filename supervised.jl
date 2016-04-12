
"""
Crossvalidation

WIP
"""
function cross_validate(features::Vector{Matrix{Float64}},
    labels::Vector{Matrix{Float64}}; method::Symbol=:leave_one_out)

end

"""
validate feature representation
"""
function validate(train::Matrix{Float64}, train_label::Vector{Int},
                  test::Matrix{Float64}, test_label::Vector{Int})
    model = build_forest(train_label, train', 4,20)#, 200, 20, 0.5)
    out = apply_forest(model, test')
    println("classification accuracy = ", mean(out.==test_label))
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
