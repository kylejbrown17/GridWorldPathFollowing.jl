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
end
# ConstSpeedStraightTrajectory
let
    t1 = 0.0
    t2 = 1.0
    interval = TimeInterval(t1,t2)
    t = 0.5
    dt = get_Δt(interval,t)
    p1 = VecE2(0.0,0.0)
    p2 = VecE2(0.0,1.0)
    traj = ConstSpeedStraightTrajectory(p1,p2,TimeInterval(t1,t2))
    arclength = norm(p2-p1)

    position = p1+dt*(p2-p1)
    heading = (p2-p1)/norm(p2-p1)
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
# ConstSpeedArcTrajectory
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
    traj = ConstSpeedArcTrajectory(center,radius,θ1,Δθ,θ2,interval)

    θ = θ1 + dt*(θ2-θ1)
    p1 = center + radius*VecE2(cos(θ1),sin(θ1))
    p2 = center + radius*VecE2(cos(θ2),sin(θ2))
    arclength = radius * abs(θ2-θ1)

    position = center+radius*VecE2(cos(θ),sin(θ))
    heading = VecE2(-sin(θ),cos(θ))
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
# Trajectory
let
    t1 = 0.0
    t2 = 1.0
    t = 0.5
    p1 = VecE2(0.0,0.0)
    p2 = VecE2(0.0,1.0)
    # traj = ConstSpeedStraightTrajectory(p1,p2,TimeInterval(t1,t2))

    t3 = 2.0
    center = VecE2(0.0,0.0)
    radius = 1.0
    θ1 = 0.0
    Δθ = π/2
    θ2 = θ1 + Δθ
    # traj = ConstSpeedArcTrajectory(center,radius,θ1,θ2,TimeInterval(t2,t3))

    traj = Trajectory(
        Vector{TrajectoryPrimitive}(
        [   ConstSpeedStraightTrajectory(p1,p2,TimeInterval(t1,t2)),
            ConstSpeedArcTrajectory(center,radius,θ1,Δθ,θ2,TimeInterval(t2,t3)) ])
        )
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
let
    start_pt = VecE2(0.0,0.0)
    start_time = 0.0
    cell_width = 1.0
    transition_time = 1.0

    action_sequence = [LEFT,LEFT,UP,UP,WAIT,RIGHT,UP,RIGHT]
    grid_path = construct_grid_world_path(start_pt,start_time,
        action_sequence,cell_width,transition_time)
    traj = construct_trajectory(grid_path)
    verify(traj)

    action_sequence = [WAIT,LEFT,LEFT,UP,DOWN,UP,WAIT,RIGHT,UP,RIGHT]
    grid_path = construct_grid_world_path(start_pt,start_time,
        action_sequence,cell_width,transition_time)
    traj = construct_trajectory(grid_path)
    verify(traj)

    action_sequence = [LEFT,UP,DOWN,UP,WAIT,RIGHT,UP,RIGHT,WAIT]
    grid_path = construct_grid_world_path(start_pt,start_time,
        action_sequence,cell_width,transition_time)
    traj = construct_trajectory(grid_path)
    verify(traj)
end
