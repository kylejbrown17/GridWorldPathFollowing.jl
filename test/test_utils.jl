let
    v = collect(1:10)
    t = 5.5
    idx,Δt = linear_interp(v,t)
    @test idx == 5

end
let
    v = collect(1:10)
    t1 = 2.5
    t2 = 5.5
    idx1,Δt1 = linear_interp(v,t1)
    idx2,Δt2 = linear_interp(v,t2)
    vp = prune_sorted_array(v,idx1,idx2,Δt1,Δt2)
    @test length(vp) == 5
end
