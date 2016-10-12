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

function pink_noise(x::Vector{Float64})
    b0::Float64 = 0
    b1::Float64 = 0
    b2::Float64 = 0
    b3::Float64 = 0
    b4::Float64 = 0
    b5::Float64 = 0
    b6::Float64 = 0

    pink::Vector{Float64} = zeros(length(x))
    for i=1:length(x)
        white = x[i]
        #from music DSP
        b0 = 0.99886 * b0 + white * 0.0555179;
        b1 = 0.99332 * b1 + white * 0.0750759;
        b2 = 0.96900 * b2 + white * 0.1538520;
        b3 = 0.86650 * b3 + white * 0.3104856;
        b4 = 0.55000 * b4 + white * 0.5329522;
        b5 = -0.7616 * b5 - white * 0.0168980;
        pink[i] = b0 + b1 + b2 + b3 + b4 + b5 + b6 + white * 0.5362;
        b6 = white * 0.115926; 
    end
    pink
end
