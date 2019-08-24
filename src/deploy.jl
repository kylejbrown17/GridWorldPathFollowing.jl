module Deploy

using Vec
using LinearAlgebra

using ..Utils
using ..GridPaths
using ..Trajectories
using ..RobotModels

export
    warmup

"""
    `warmup()`

    This function is meant to be called before running a speed-critical
    application, so that all required functions will be compiled prior to
    the first call in the core application.
"""
function warmup()
    println("GridWorldPathFollowing warming up ...")
    pt0 = VecE2(0.0,0.0) # start position
    t0 = 0.0 # start time
    actions = [SOUTH,EAST,WAIT,NORTH,EAST,WAIT,WAIT,EAST,EAST,SOUTH]
    w_cell = 1.0 # cell width
    t_cell = 2.0 # transition time
    grid_path1 = construct_grid_world_path(pt0,t0,actions,w_cell,t_cell)
    base_traj1 = construct_trajectory(grid_path1)
    Δt = 4.0
    grid_path2 = construct_grid_world_path(
        get_end_pt(base_traj1),get_end_time(base_traj1)+Δt,
        [EAST,NORTH,WAIT,EAST,SOUTH],w_cell,t_cell)
    base_traj2 = construct_trajectory(grid_path2)
    base_traj = stitch_trajectories(base_traj1,base_traj2;buffer=[0.1,0.1])
    # optimize velocity profile
    # traj, t_vec, accel, vel, pos = optimize_velocity_profile(base_traj);
    traj = optimize_velocity_profile_traj_only(base_traj)
    # simulate closed-loop controller
    sim_model = UnicycleModel()
    controller = SwitchingController()
    construct_unicycle_controller(sim_model,grid_path1,controller)
    tf = get_end_time(traj)
    dt = 0.05
    state_pt = get_trajectory_point_by_time(traj, t0)
    initial_state = [pt0.x, pt0.y, atan(state_pt.heading)] + .1*(1 .- rand(3))
    states, cmds = simulate(sim_model,controller,traj,initial_state,t0,tf,dt)
    println("GridWorldPathFollowing warm-up complete!")
    return
end

end
