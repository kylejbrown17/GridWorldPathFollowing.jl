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
    t_vec = Vector{Float64}()
    vel = Vector{Float64}()
    pos = Vector{Float64}()
    s = 0.0
    for seg in traj_.segments
        push!(t_vec, get_start_time(seg))
        push!(t_vec,   get_end_time(seg))
        push!(vel, norm(get_vel(seg,get_start_time(seg))))
        push!(vel, norm(get_vel(seg,get_end_time(seg))))
        push!(pos, s)
        s += get_length(seg)
        push!(pos, s)
    end
    accel = zeros(length(t_vec)-1)
    traj = DenseTrajectory(traj_,t_vec,accel,vel,pos)

    controller = SwitchingController()
    state = [0.0,0.0,0.0]
    t = 0.1
    u = get_action(controller,traj,state,t)

    model = UnicycleModel()
    dt = 0.1
    simulate(model,controller,traj,state,0.0,2.0,dt)
end
