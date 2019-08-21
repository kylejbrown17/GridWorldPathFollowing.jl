let
    model = UnicycleKinematicModel(1.0,2.0)
    state = UnicycleKinematicState(VecE2(0.0,0.0),VecE2(0.0,1.0),0.0)
    action = UnicycleKinematicCommand(-1.0,1.0)
    dt = 1.0
    nsteps=10
    @show state.pos
    @show state.heading
    for i in 1:12
        state = next_state(model,state,action,dt;nsteps=nsteps)
    end
    @show state.pos
    @show state.heading
end
