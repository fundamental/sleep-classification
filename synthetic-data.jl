#1000x4000

"""
Generate possible modes for the system to run with
"""
function gen_modes(vals::Vector{Float64}, modes::Int, rows::Int)
    N = rows
    M = modes
    out = zeros(N,M)

    for i=1:N, j=1:M
        out[i,j] = rand(vals)
    end

    return out
end

"""
Generate a random number from 'options' which does not include element
'excluded'
"""
function rand_without{T}(options::Vector{T}, excluded::T)
    out = excluded
    while(out == excluded)
        out = rand(options)
    end
    out
end

function make_irregular_grid(;modes::Int=6, width::Int=4000, height::Int=1000,
                              cols::Int=64, rows::Int=8,
                              levels::Vector{Float64}=[-1.0,0.0,1.0])
    random_modes = gen_modes(levels, modes, rows) 
    #random_modes = gen_modes(rand(6), modes, rows) 

    cdiv = zeros(Int, cols-1)
    rdiv = zeros(Int, rows-1)

    for i=1:cols-1
        cdiv[i] = rand(collect(1:width))
    end

    for i=1:rows-1
        rdiv[i] = rand(collect(1:height))
    end

    cdiv = sort(cdiv)
    rdiv = sort(rdiv)

    col_map = zeros(Int, width)
    row_map = zeros(Int, height)

    current = 1
    for i=1:width
        if(sum(cdiv.==i) != 0)
            current = rand_without([1:modes;], current)
        end
        col_map[i] = current
    end

    current = 1
    for i=1:height
        if(sum(rdiv.==i) != 0)
            current += 1
        end
        row_map[i] = current
    end

    out = zeros(height, width)
    for i=1:height, j=1:width
        #println("grabbing [", row_map[i], ",", col_map[j], "]")
        out[i,j] = random_modes[row_map[i], col_map[j]]
    end
    return out
end

