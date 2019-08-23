let
    t1 = 0.0
    t2 = 1.0
    t = 0.5
    @test_throws AssertionError verify(TimeInterval(t2,t1))
    interval = TimeInterval(t1,t2)
    @test get_start_time(interval) == t1
    @test get_end_time(interval) == t2
    @test get_Δt(interval,t) == (t - t1) / (t2 - t1)
end
let
    TrajectoryPoint()
    @test_throws AssertionError verify(TrajectoryPoint(t=NaN))

    x = [0.0,1.0]
    y = [1.0,0.0]
    θ = [0.0,π]
    vx = [0.0,1.0]
    vy = [1.0,0.0]
    ω = [1.0,2.0]
    t = [0.0,1.0]
    pt1 = TrajectoryPoint(
        pos=VecE2(x[1],y[1]),
        heading=VecE2(cos(θ[1]),sin(θ[1])),
        vel = VecE2(vx[1],vy[1]),
        yaw_rate = ω[1],
        t = t[1]
        )
    pt2 = TrajectoryPoint(
        pos=VecE2(x[2],y[2]),
        heading=VecE2(cos(θ[2]),sin(θ[2])),
        vel = VecE2(vx[2],vy[2]),
        yaw_rate = ω[2],
        t = t[2]
        )
    t_idx = 0.5
    pt = interpolate(pt1,pt2,t_idx)
    @test array_isapprox(pt.pos, (1-t_idx)*pt1.pos + t_idx*pt2.pos)
    @test isapprox(atan(pt.heading), (1-t_idx)*θ[1] + t_idx*θ[2])
    @test array_isapprox(pt.vel, (1-t_idx)*pt1.vel + t_idx*pt2.vel)
    @test isapprox(pt.yaw_rate, (1-t_idx)*pt1.yaw_rate + t_idx*pt2.yaw_rate)
    @test isapprox(pt.t, (1-t_idx)*pt1.t + t_idx*pt2.t)
end
# StraightTrajectory
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
    traj = StraightTrajectory(p1,heading,arclength,TimeInterval(t1,t2))

    position = p1+dt*(p2-p1)
    velocity = (p2-p1)/(t2-t1)
    yaw_rate = 0.0

    @test get_start_time(traj) == get_start_time(interval)
    @test get_end_time(traj) == get_end_time(interval)
    @test get_Δt(traj,t) == dt

    @test get_time_from_arc_length(traj,arclength*(t-t1)/(t2-t1)) == t

    @test get_length(traj) == arclength
    @test get_dist(traj, (t1+t2)/2) == arclength/2

    @test array_isapprox(get_start_pt(traj),p1)
    @test array_isapprox(get_end_pt(traj),p2)
    @test array_isapprox(get_position(traj,t), position)
    @test array_isapprox(get_heading(traj,t), heading)
    @test array_isapprox(get_vel(traj,t), velocity)
    @test get_yaw_rate(traj,t) == yaw_rate
    @test get_time_from_pt(traj,get_position(traj,t)) == t

    traj_pt1 = get_trajectory_point_by_time(traj,t)
    @test array_isapprox(traj_pt1.pos, position)
    @test array_isapprox(traj_pt1.heading, heading)
    @test array_isapprox(traj_pt1.vel, velocity)
    @test traj_pt1.yaw_rate == yaw_rate

    traj_pt2 = get_trajectory_point_by_pt(traj, get_position(traj,t))
    @test array_isapprox(traj_pt2.pos, position)
    @test array_isapprox(traj_pt2.heading, heading)
    @test array_isapprox(traj_pt2.vel, velocity)
    @test traj_pt2.yaw_rate == yaw_rate
end
# ArcTrajectory
let
    t1 = 0.0
    t2 = 1.0
    interval = TimeInterval(t1,t2)
    t = 0.5
    dt = get_Δt(interval,t)

    center = VecE2(0.0,0.0)
    radius = 1.0
    θ1 = 0.0
    Δθ = π/2
    θ2 = θ1 + Δθ
    heading1 = VecE2(-sin(θ1),cos(θ1))
    if sign(Δθ) != 0
        heading1 = VecE2(-sin(θ1),cos(θ1)) * sign(Δθ)
    end
    θ = θ1 + dt*(θ2-θ1)
    p1 = center + radius*VecE2(cos(θ1),sin(θ1))
    p2 = center + radius*VecE2(cos(θ2),sin(θ2))
    arclength = radius * abs(θ2-θ1)

    traj = ArcTrajectory(center,p1,heading1,Δθ,interval)
    # traj1 = ArcTrajectory(center,radius,θ1,Δθ,interval)
    # @test traj1.θ1 == traj.θ1

    position = center+radius*VecE2(cos(θ),sin(θ))
    heading = VecE2(-sin(θ),cos(θ))
    if sign(Δθ) != 0
        heading = VecE2(-sin(θ),cos(θ)) * sign(Δθ)
    end
    velocity = heading * arclength / (t2-t1)
    yaw_rate = (θ2-θ1) / (t2-t1)

    @test get_start_time(traj) == get_start_time(interval)
    @test get_end_time(traj) == get_end_time(interval)
    @test get_Δt(traj,t) == dt

    @test get_length(traj) == arclength
    @test get_dist(traj, (t1+t2)/2) == arclength/2

    @test array_isapprox(get_start_pt(traj),p1)
    @test array_isapprox(get_end_pt(traj),p2)
    @test array_isapprox(get_position(traj,t), position)
    @test array_isapprox(get_heading(traj,t), heading)
    @test array_isapprox(get_vel(traj,t), velocity)
    @test get_yaw_rate(traj,t) == yaw_rate
    # @test get_time_from_pt(traj,get_position(traj,t)) == t

    traj_pt1 = get_trajectory_point_by_time(traj,t)
    @test array_isapprox(traj_pt1.pos, position)
    @test array_isapprox(traj_pt1.heading, heading)
    @test array_isapprox(traj_pt1.vel, velocity)
    @test traj_pt1.yaw_rate == yaw_rate

    traj_pt2 = get_trajectory_point_by_pt(traj, get_position(traj,t))
    @test array_isapprox(traj_pt2.pos, position)
    @test array_isapprox(traj_pt2.heading, heading)
    @test array_isapprox(traj_pt2.vel, velocity)
    @test traj_pt2.yaw_rate == yaw_rate
end
# WaitTrajectory
let
    t1 = 0.0
    t2 = 1.0
    interval = TimeInterval(t1,t2)
    t = 0.5
    dt = get_Δt(interval,t)
    p1 = VecE2(0.0,0.0)
    p2 = p1
    heading = VecE2(0.0,1.0)
    arclength = norm(p2-p1)
    traj = WaitTrajectory(p1,heading,TimeInterval(t1,t2))

    position = p1+dt*(p2-p1)
    velocity = (p2-p1)/(t2-t1)
    yaw_rate = 0.0

    @test get_start_time(traj) == get_start_time(interval)
    @test get_end_time(traj) == get_end_time(interval)
    @test get_Δt(traj,t) == dt

    @test get_length(traj) == arclength
    @test get_dist(traj, (t1+t2)/2) == arclength/2

    @test array_isapprox(get_start_pt(traj),p1)
    @test array_isapprox(get_end_pt(traj),p2)
    @test array_isapprox(get_position(traj,t), position)
    @test array_isapprox(get_heading(traj,t), heading)
    @test array_isapprox(get_vel(traj,t), velocity)
    @test get_yaw_rate(traj,t) == yaw_rate
    @test get_time_from_pt(traj,get_position(traj,t)) == get_start_time(traj)

    traj_pt1 = get_trajectory_point_by_time(traj,t)
    @test array_isapprox(traj_pt1.pos, position)
    @test array_isapprox(traj_pt1.heading, heading)
    @test array_isapprox(traj_pt1.vel, velocity)
    @test traj_pt1.yaw_rate == yaw_rate

    traj_pt2 = get_trajectory_point_by_pt(traj, get_position(traj,t))
    @test array_isapprox(traj_pt2.pos, position)
    @test array_isapprox(traj_pt2.heading, heading)
    @test array_isapprox(traj_pt2.vel, velocity)
    @test traj_pt2.yaw_rate == yaw_rate
end
# PivotTrajectory
let
    t1 = 0.0
    t2 = 1.0
    interval = TimeInterval(t1,t2)
    t = 0.5
    dt = get_Δt(interval,t)

    center = VecE2(0.0,0.0)
    radius = 0.0
    θ1 = 0.0
    Δθ = 0.0
    θ2 = θ1 + Δθ
    heading1 = VecE2(-sin(θ1),cos(θ1))
    if sign(Δθ) != 0
        heading1 = VecE2(-sin(θ1),cos(θ1)) * sign(Δθ)
    end
    θ = θ1 + dt*(θ2-θ1)
    p1 = center + radius*VecE2(cos(θ1),sin(θ1))
    p2 = center + radius*VecE2(cos(θ2),sin(θ2))
    arclength = radius * abs(θ2-θ1)

    traj = PivotTrajectory(p1,heading1,Δθ,interval)

    position = center+radius*VecE2(cos(θ),sin(θ))
    heading = VecE2(-sin(θ),cos(θ))
    if sign(Δθ) != 0
        heading = VecE2(-sin(θ),cos(θ)) * sign(Δθ)
    end
    velocity = heading * arclength / (t2-t1)
    yaw_rate = (θ2-θ1) / (t2-t1)

    @test get_start_time(traj) == get_start_time(interval)
    @test get_end_time(traj) == get_end_time(interval)
    @test get_Δt(traj,t) == dt

    @test get_length(traj) == arclength
    @test get_dist(traj, (t1+t2)/2) == arclength/2

    @test array_isapprox(get_start_pt(traj),p1)
    @test array_isapprox(get_end_pt(traj),p2)
    @test array_isapprox(get_position(traj,t), position)
    @test array_isapprox(get_heading(traj,t), heading)
    @test array_isapprox(get_vel(traj,t), velocity)
    @test get_yaw_rate(traj,t) == yaw_rate
    # @test get_time_from_pt(traj,get_position(traj,t)) == t

    traj_pt1 = get_trajectory_point_by_time(traj,t)
    @test array_isapprox(traj_pt1.pos, position)
    @test array_isapprox(traj_pt1.heading, heading)
    @test array_isapprox(traj_pt1.vel, velocity)
    @test traj_pt1.yaw_rate == yaw_rate

    traj_pt2 = get_trajectory_point_by_pt(traj, get_position(traj,t))
    @test array_isapprox(traj_pt2.pos, position)
    @test array_isapprox(traj_pt2.heading, heading)
    @test array_isapprox(traj_pt2.vel, velocity)
    @test traj_pt2.yaw_rate == yaw_rate
end
# Trajectory
let
    t1 = 0.0
    t2 = 1.0
    t = 0.5
    p1 = VecE2(0.0,0.0)
    p2 = VecE2(0.0,1.0)
    heading = (p2-p1)/norm(p2-p1)
    arclength = norm(p2-p1)
    seg1 = StraightTrajectory(p1,heading,arclength,TimeInterval(t1,t2))

    t3 = 2.0
    center = VecE2(-1.0,1.0)
    radius = 1.0
    θ1 = 0.0
    Δθ = π/2
    θ2 = θ1 + Δθ
    seg2 = ArcTrajectory(center,p2,get_heading(seg1,get_end_time(seg1)),Δθ,TimeInterval(t2,t3))

    traj = Trajectory()
    push!(traj,seg1)
    push!(traj,seg2)
    verify(traj)

    get_dist(traj,t1)
    get_dist(traj,t3)
    get_position(traj,t1)
    get_position(traj,t3)
    get_heading(traj,t1)
    get_heading(traj,t3)
    get_vel(traj,t1)
    get_vel(traj,t3)
    get_yaw_rate(traj,t1)
    get_yaw_rate(traj,t3)
end
# hierarchical Trajectory
let
end
# construct_trajectory
let
    start_pt = VecE2(0.0,0.0)
    start_time = 0.0
    cell_width = 1.0
    transition_time = 1.0

    action_sequence = [WEST,NORTH]
    grid_path = construct_grid_world_path(start_pt,start_time,
        action_sequence,cell_width,transition_time)
    traj = construct_trajectory(grid_path)
    verify(traj)

    action_sequence = [WEST,WEST,NORTH,NORTH,WAIT,EAST,NORTH,EAST]
    grid_path = construct_grid_world_path(start_pt,start_time,
        action_sequence,cell_width,transition_time)
    traj = construct_trajectory(grid_path)
    verify(traj)

    action_sequence = [WAIT,WEST,WEST,NORTH,SOUTH,NORTH,WAIT,EAST,NORTH,EAST]
    grid_path = construct_grid_world_path(start_pt,start_time,
        action_sequence,cell_width,transition_time)
    traj = construct_trajectory(grid_path)
    verify(traj)

    action_sequence = [WEST,NORTH,SOUTH,NORTH,WAIT,EAST,NORTH,EAST,WAIT]
    grid_path = construct_grid_world_path(start_pt,start_time,
        action_sequence,cell_width,transition_time)
    traj = construct_trajectory(grid_path)
    verify(traj)
end
# DenseTrajectory
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

    @test get_time_from_arc_length(traj_,arclength*(t-t1)/(t2-t1)) == t
    # @test get_time_from_arc_length(traj,arclength/2) == t
    @test get_index_time(traj,t) == t

    @test get_start_time(traj) == get_start_time(interval)
    @test get_dist(traj, (t1+t2)/2) == arclength/2

    @test array_isapprox(get_start_pt(traj),p1)
    @test array_isapprox(get_position(traj,t), pos_vec)
    @test array_isapprox(get_heading(traj,t), heading)
    @test array_isapprox(get_vel(traj,t), velocity)
    @test isapprox(get_yaw_rate(traj,t), yaw_rate)

    traj_pt1 = get_trajectory_point_by_time(traj,t)
    @test array_isapprox(traj_pt1.pos, pos_vec)
    @test array_isapprox(traj_pt1.heading, heading)
    @test array_isapprox(traj_pt1.vel, velocity)
    @test isapprox(traj_pt1.yaw_rate, yaw_rate)

    # Now verify that the yaw_rate matches in the zero-velocity segment
    @test isapprox(get_yaw_rate(traj,t2_mid), get_yaw_rate(traj_,t2_mid))
end
# Speed Profile Optimization
let
    # start_pt = VecE2(0.0,0.0)
    # start_time = 0.0
    # action_sequence = [SOUTH,EAST,WAIT,NORTH,EAST,EAST,EAST,SOUTH,SOUTH,WEST,NORTH,WEST,SOUTH,WEST,WEST]
    # # action_sequence = [EAST,NORTH,EAST]
    # cell_width = 1.0
    # transition_time = 2.0
    # grid_path = construct_grid_world_path(start_pt,start_time,
    #     action_sequence,cell_width,transition_time)
    #
    # traj = construct_trajectory(grid_path)
    # verify(traj)
    #
    # t_vec, accel, vel, pos = optimize_velocity_profile(traj)
    #
    # dense_traj = DenseTrajectory(traj,t_vec,accel,vel,pos)
    # verify(dense_traj)
    # t = 0.5
    # get_length(dense_traj)
    # get_dist(dense_traj,t)
    # get_position(dense_traj,t)
    # get_heading(dense_traj,t)
    # get_vel(dense_traj,t)
    # get_yaw_rate(dense_traj,t)
    #
    # get_trajectory_point_by_time(dense_traj,t)
end
let
    t1 = 0.0
    t2 = 1.0
    t = 0.5
    p1 = VecE2(0.0,0.0)
    p2 = VecE2(0.0,1.0)
    heading = (p2-p1)/norm(p2-p1)
    arclength = norm(p2-p1)
    seg1 = StraightTrajectory(p1,heading,arclength,TimeInterval(t1,t2))

    t3 = 2.0
    center = VecE2(-1.0,1.0)
    radius = 1.0
    θ1 = 0.0
    Δθ = π/2
    θ2 = θ1 + Δθ
    seg2 = ArcTrajectory(center,p2,get_heading(seg1,get_end_time(seg1)),Δθ,TimeInterval(t2,t3))

    traj_ = Trajectory()
    push!(traj_,seg1)
    push!(traj_,seg2)
    verify(traj_)

    Δt = (get_end_time(traj_)-get_start_time(traj_))/100
    SimpleTrajectory(traj_,100)
    traj = SimpleTrajectory(traj_,Δt)
    verify(traj)
    get_trajectory_point_by_time(traj,t)
    get_dist(traj,t)
    get_position(traj,t)
    get_vel(traj,t)
    get_heading(traj,t)
    get_yaw_rate(traj,t)
end
