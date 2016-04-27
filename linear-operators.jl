function kClassify(data::Matrix{Float64}, centers::Int; plotID::Int=-1,
    maxiter::Int=40)
    k = kmeans(data, centers)

    res = Float64[]
    for i=1:size(data)[2]
        choice = indmin(mapslices(norm, k.centers.-data[:,i],1))
        push!(res, choice)
    end

    if(plotID >= 0)
        figure(plotID)
        plot(res+0.01randn(length(res)))
        figure(plotID+1)
        plot(res[1:end-1].!=res[2:end])
    end
    res
    (k.centers, res)
end

function kClassify_nway_memo(data, centers,
memo_data::Tuple{Matrix{Float64}, Matrix{Float64},
    Vector{Int}, Vector{Float64}}
    ; maxiter::Int=20, plotID::Int=-1,
    )
    #print("F")

    col_id::Vector{Int} = memo_data[3]
    cur_id::Int         = maximum(col_id)

    fast::Matrix{Float64}    = memo_data[2]
    weights::Vector{Float64} = memo_data[4]
    k = kmeans(fast, centers, weights=weights, maxiter=maxiter)
    real_ans = zeros(size(data,2))
    for i=1:size(data,2)
        real_ans[i] = k.assignments[col_id[i]]
    end
    return real_ans
end
    
function kClassify_nway(data, centers; plotID::Int=-1, maxiter::Int=20)
    col_id::Vector{Int} = ones(Int, size(data,2))
    cur_id::Int         = 1
    for i=2:size(data,2)
        if(data[:,i] != data[:,i-1])
            cur_id += 1
        end
        col_id[i] = cur_id
    end
    #println("same col = ", same_col, " of ", size(data,2))

    #Fast path
    ev_id = 1
    if(cur_id < 0.4size(data,2))
        #print("f")
        fast = zeros(size(data,1), cur_id)
        weights = zeros(cur_id)
        for i=1:cur_id
            weights[i] = length(find(col_id.==i))
        end
        for i=1:size(data,2)
            if(col_id[i] == ev_id)
                fast[:,ev_id] = data[:,i]
                ev_id += 1
            end
        end
        k = kmeans(fast, centers, weights=weights, maxiter=maxiter)
        real_ans = zeros(size(data,2))
        for i=1:size(data,2)
            real_ans[i] = k.assignments[col_id[i]]
        end
        return (real_ans, (data, fast, col_id, weights))
    else
        #print("s(", cur_id./(1.0*size(data,2)),")")
        k = kmeans(data, centers, maxiter=maxiter)
        return (k.assignments, nothing)
    end



    k = kmeans(data, centers, maxiter=maxiter)
    return k.assignments
    freq = zeros(centers)
    res = Float64[]
    for i=1:size(data)[2]
        choice = k.assignments[i];#indmin(sum((k.centers.-data[:,i]).^2,1))
        push!(res, choice)
        freq[choice] += 1
    end
    per = reverse(sortperm(freq))
    for i=1:size(data)[2]
        res[i] = per[round(Int,res[i])]
    end

    #figure(plotID)
    #plot(res+0.01randn(length(res)))
    #figure(plotID+1)
    #plot(res[1:end-1].!=res[2:end])
    res
end

"""
Perform a median over a 1D vector

Arguments:

- vec - 1D vector
- med - window of the vector
"""
function linearMedian{T}(vec::Vector{T}, med::Int)
    #return img
    @assert isodd(med)
    scale = round(Int,(med-1)/2)
    result = copy(vec)
    N = length(vec)
    getIt(x,N) = Int(x<1?1:x>N?N:x)

    opts = Array(Float64, med)
    for x=1:length(vec)
        ii = 1
        if(x < med || x > (N-med))
            for i=-scale:scale
                opts[ii] = vec[getIt(x+i,N)]
                ii += 1
            end
            result[x] = median(opts)
        else
            prev = result[x-1]
            if(sign(vec[x-(scale+1)] - prev) == sign(vec[x+scale] - prev))
                result[x] = prev
            else
                for i=-scale:scale
                    opts[ii] = vec[x+i]
                    ii += 1
                end
                result[x] = median(opts)
            end
        end
    end
    result
end

