module Utils

using GraphUtils

export
    interpolate,
    linear_interp,
    prune_sorted_array

interpolate(a,b,t) = (1 - t)*a + t*b

function linear_interp(v,t)
    idx = find_index_in_sorted_array(v,t)-1
    idx = max(1, min(length(v)-1, idx))
    if isapprox(v[idx+1] - v[idx],0.0)
        return idx,0.0
    end
    Δt = (t - v[idx]) / (v[idx+1] - v[idx])
    return idx, Δt
end

"""
    `prune_array(arr,v,t1,t2)`

    prunes `arr` on either side based on where t1 and t2 fall in v.
"""
function prune_sorted_array(v,idx1,idx2,Δt1,Δt2)
    @assert [idx1,Δt1] <= [idx2,Δt2] "[idx1,Δt1] > [idx2,Δt2]"
    vp = [interpolate(v[idx1],v[idx1+1],Δt1)]
    for idx in min(idx1+1,length(v)):min(idx2,length(v))
        push!(vp,v[idx])
    end
    push!(vp,interpolate(v[idx2],v[idx2+1],Δt2))
    return vp
end

end
