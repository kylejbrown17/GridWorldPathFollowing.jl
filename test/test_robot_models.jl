let
    start_pt = VecE2(0.0,0.0)
    start_time = 0.0
    cell_width = 1.0
    transition_time = 1.0

    action_sequence = [WEST,WEST,NORTH,NORTH,WAIT,EAST,NORTH,EAST]
    grid_path = construct_grid_world_path(start_pt,start_time,
        action_sequence,cell_width,transition_time)
    traj = construct_trajectory(grid_path)
    verify(traj)

    controller = SwitchingController()
    state = [0.0,0.0,0.0]
    t = 0.1
    u = get_action(controller,traj,state,t)

    model = UnicycleModel()
    dt = 0.1
    simulate(model,controller,traj,state,0.0,2.0,dt)
end
let
    t1 = 0.0
    t2 = 1.0
    interval = TimeInterval(t1,t2)
    t = 0.5
    dt = get_Δt(interval,t)
    p1 = VecE2(0.0,0.0)
    p2 = VecE2(0.0,1.0)
    heading = (p2-p1)/norm(p2-p1)
    arclength = norm(p2-p1)
    seg1 = StraightTrajectory(p1,heading,arclength,TimeInterval(t1,t2))

    pos_vec = p1+dt*(p2-p1)
    velocity = (p2-p1)/(t2-t1)
    yaw_rate = 0.0

    t3 = 2.0
    t2_mid = (t3+t2)/2
    center = VecE2(0.0,1.0)
    radius = 0.0
    θ1 = 0.0
    Δθ = π/2
    θ2 = θ1 + Δθ
    seg2 = ArcTrajectory(center,p2,get_heading(seg1,get_end_time(seg1)),Δθ,TimeInterval(t2,t3))

    traj_ = Trajectory()
    push!(traj_,seg1)
    push!(traj_,seg2)
    verify(traj_)

    # Construct a DenseTrajectory that encodes the exact same information as the
    # piecewise constant Trajectory
    traj = Trajectory()
    for seg in traj_.segments
        t_vec = [get_start_time(seg), get_end_time(seg)]
        accel = [0.0]
        vel = [norm(get_vel(seg,get_start_time(seg))), norm(get_vel(seg,get_end_time(seg)))]
        pos = [0.0, get_length(seg)]
        push!(traj,DenseTrajectory(seg,t_vec,accel,vel,pos))
    end

    controller = SwitchingController()
    state = [0.001,0.001,0.001]
    t = 0.1
    u = get_action(controller,traj,state,t)

    model = UnicycleModel()
    dt = 0.1
    simulate(model,controller,traj,state,0.01,2.0,dt)
end
