
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
