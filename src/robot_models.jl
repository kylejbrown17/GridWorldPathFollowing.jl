module RobotModels

using Parameters
using Vec, LinearAlgebra

using ..Utils
using ..Trajectories
using ..GridPaths

export
    forward_euler_integration,
    sat,

    UnicycleModel,
    UnicycleState,
    UnicycleAction,

    dynamics,

    Controller,
    TrackingController,
    StabilizeController,
    PivotController,
    SwitchingController,

    get_action,
    simulate

"""
    `forward_euler_integration`
"""
function forward_euler_integration(model,state,u,dt;nsteps=10)
    s = deepcopy(state)
    for i in 1:nsteps
        sdot = dynamics(model,s,u)
        s += sdot*dt/nsteps
    end
    return s
end

""" saturation function """
sat(x,δ) = abs(x) <= δ ? x : sign(x)*δ

abstract type AbstractRobotModel end

struct UnicycleModel <: AbstractRobotModel end
const UnicycleState = Vector{Float64}
const UnicycleAction = Vector{Float64}
function dynamics(model::UnicycleModel,state,cmd)
    x = state[1]
    y = state[2]
    θ = state[3]
    w = cmd[1]
    v = cmd[2]
    sdot = [v*cos(θ),v*sin(θ),w]
end
################################################################################
################################## Controllers #################################
################################################################################

abstract type Controller end

################################################################################
############################## TrackingController ##############################
################################################################################
"""
    `TrackingController`

    To be employed when the robot needs to track a time-varying reference signal

    From Lee et al., "Tracking Control of
    Unicycle-Modeled Mobile Robots Using a Saturation Feedback Controller"
    [link]{https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=911382}
"""
@with_kw struct TrackingController <: Controller
    dt  ::Float64 = 0.0
    a   ::Float64 = 4.0
    k0  ::Float64 = 4.0
    μ   ::Float64 = 1.0
    γ   ::Float64 = 0.1
    ϵ   ::Float64 = 0.1
    k1  ::Float64 = 4.0
    b   ::Float64 = 4.0
end


################################################################################
############################### PivotController ################################
################################################################################
"""
    `PivotController`

    To be employed when the robot needs to turn in place.
"""
@with_kw struct PivotController <: Controller
    kw::Float64 = 0.5
    kp::Float64 = 2.0
end
################################################################################
############################# StabilizeController ##############################
################################################################################
"""
    `StabilizeController`

    To be employed when the Robot needs to stabilize about a particular point
"""
@with_kw struct StabilizeController <: Controller
    kw::Float64 = 2.0
    kp::Float64 = 4.0
    kϕ::Float64 = 15.0
end

################################################################################
############################# SwitchingController ##############################
################################################################################
"""
    SwitchingController

    A composite controller that switches between operating modes based on the
    active trajectory type.
"""
@with_kw struct SwitchingController <: Controller
    tracker     ::TrackingController    = TrackingController()
    pivoter     ::PivotController       = PivotController()
    stabilizer  ::StabilizeController   = StabilizeController()
end

"""
    This function follows Lee et al., "Tracking Control of
    Unicycle-Modeled Mobile Robots Using a Saturation Feedback Controller"

    [link]{https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=911382}

    Inputs:
    * `target` - the current reference point [xr, yr, θr] of the tracked trajectory
    * `ff`     - the current feedforward command [vr, wr] of the tracked trajectory
    * `state`  - the current state [x, y, θ] of the robot
"""
function get_action(controller::TrackingController,target::UnicycleState,ff::UnicycleAction,state::UnicycleState,t::Float64)
    # println("get_action(controller::TrackingController,target::UnicycleState,ff::UnicycleAction,state::UnicycleState,t::Float64)")
    a = controller.a   # 0 < a < vmax - sup_{t >= 0} abs(vr)
    k0 = controller.k0 # k0 > 0
    μ = controller.μ   # μ = 0 or 1
    γ = controller.γ   # 0 < γ < 1
    ϵ = controller.ϵ   # 0 < ϵ < 1/(1+γ)
    k1 = controller.k1 # k1 > 0
    b = controller.b   # b > 0
    # state (global frame)
    x,y,θ = state[1],state[2],state[3]
    # reference (global frame)
    xr,yr,θr = target[1],target[2],target[3]
    wr,vr = ff[1],ff[2]
    # errors in robot frame
    e = [ cos(θ) sin(θ) 0;
         -sin(θ) cos(θ) 0;
              0      0  1]*[xr-x,yr-y,get_angular_offset(θ,θr)]
    xe,ye,θe = e[1],e[2],e[3]
    # dirty trick to avoid problems with waiting
    # if (abs(vr) <= 0.000001) && (abs(wr) <= 0.000001)
    #     return [0.0,0.0]
    # end
    # change of coordinates
    x0 = θe
    x1 = ye
    x2 = -xe
    # u0 = wr - w
    # u1 = v - vr*cos(x0)
    V1 = x1^2 + x2^2
    u1 = -sat(k0*x2,a)
    h = 1 + γ*cos(μ*t)
    h_dot = -μ*γ*sin(μ*t)
    x0_ = x0 + ϵ*h*x1 / (1 + V1^(1/2))
    α = 1 - ϵ*h*x2 / (1 + V1^(1/2))
    β = ϵ*(
        (h_dot*x1+h*wr*x2+h*vr*sin(x0))/(1 + V1^(1/2))
        - (h*x1/( ((1 + V1^(1/2))^2)*V1^(1/2) ))
        * (x1*vr*sin(x0) - sat(k0*x2,a)*x2)
    )
    u0 = -β/α - sat(k1*x0_,b)
    # transform back
    w = -u0 + wr
    v = u1 + vr*cos(x0)
    return [w, v]
end
"""
    `get_action(controller::PivotController,target,ff,state,t)`

    A simple control law to track about a time-varying heading while reducing
    the position error when possible. The idea is to track the heading signal
    and drive forward or backward to reduce the position error.

    Inputs:
    * `target` - the current reference point [xr, yr, θr] of the tracked trajectory
    * `ff`     - the current feedforward command [vr, wr] of the tracked trajectory
    * `state`  - the current state [x, y, θ] of the robot
"""
function get_action(controller::PivotController,target::Vector{Float64},ff::Vector{Float64},state::Vector{Float64},t::Float64)
    # println("get_action(controller::PivotController,target::Vector{Float64},ff::Vector{Float64},state::Vector{Float64},t::Float64)")
    kw = controller.kw
    kp = controller.kp
    # state (global frame)
    x,y,θ = state[1],state[2],state[3]
    # reference (global frame)
    xr,yr,θr = target[1],target[2],target[3]
    wr,vr = ff[1],ff[2]
    # control law
    dp = dot([xr-x,yr-y],[cos(θ),sin(θ)]) # distance to target position projected onto heading vector
    dθ = get_angular_offset(θ,θr) # heading error
    w = wr + kw*dθ
    v = vr + kp*dp
    return [w, v]
end
"""
    `get_action(controller::stabilizeController,target,ff,state,t)`

    A simple control law to stabilize about a static pose. 

    Inputs:
    * `target` - the current reference point [xr, yr, θr] of the tracked trajectory
    * `ff`     - the current feedforward command [vr, wr] of the tracked trajectory
    * `state`  - the current state [x, y, θ] of the robot
"""
function get_action(controller::StabilizeController,target::Vector{Float64},ff::Vector{Float64},state::Vector{Float64},t::Float64)
    # println("get_action(controller::StabilizeController,target::Vector{Float64},ff::Vector{Float64},state::Vector{Float64},t::Float64)")
    kw = controller.kw
    kp = controller.kp
    kϕ = controller.kϕ
    # state (global frame)
    x,y,θ = state[1],state[2],state[3]
    # reference (global frame)
    xr,yr,θr = target[1],target[2],target[3]
    wr,vr = ff[1],ff[2]
    # control law
    dx = xr-x
    dy = yr-y
    dp = dot([dx,dy],[cos(θ),sin(θ)]) # distance to target position projected onto heading vector
    dθ = ([cos(θ),sin(θ),0] × [cos(θr),sin(θr),0])[end] # heading error
    ϕ = atan(dy,dx)
    dϕ = ([cos(θ),sin(θ),0] × [cos(ϕ),sin(ϕ),0])[end] # angular displacement from target position
    w = wr + kϕ*dϕ*dp + kw*dθ
    v = vr + kp*dp
    return [w, v]
end

function get_action(controller::C,ref::TrajectoryPoint,state::UnicycleState,t::Float64) where {C<: Controller}
    # println("get_action(controller::C,ref::TrajectoryPoint,state::UnicycleState,t::Float64) where {C<: Controller}")
    target = [ref.pos.x, ref.pos.y, atan(ref.heading)]
    ff = [ref.yaw_rate, norm(ref.vel)]
    get_action(controller,target,ff,state,t)
end
function get_action(controller::C,traj::T,state::UnicycleState,t::Float64) where {C<: Controller, T <: AbstractTrajectory}
    # println("get_action(controller::C,traj::T,state::UnicycleState,t::Float64) where {C<: Controller, T <: AbstractTrajectory}")
    ref = get_trajectory_point_by_time(traj,t)
    get_action(controller,ref,state,t)
end

function get_action(controller::SwitchingController,traj::Trajectory,state::UnicycleState,t::Float64)
    # println("get_action(controller::SwitchingController,traj::Trajectory,state::UnicycleState,t::Float64)")
    get_action(controller,get_active_segment(traj,t),state,t)
end
const TrackingTraj = Union{StraightTrajectory,ArcTrajectory,DenseTrajectory{StraightTrajectory},DenseTrajectory{ArcTrajectory}}
function get_action(controller::SwitchingController,traj::T,state::UnicycleState,t::Float64) where {T<:TrackingTraj}
    # println("get_action(controller::SwitchingController,traj::DenseTrajectory,state::UnicycleState,t::Float64)")
    get_action(controller.tracker,traj,state,t)
end
const WaitTraj = Union{WaitTrajectory,DenseTrajectory{WaitTrajectory}}
function get_action(controller::SwitchingController,traj::T,state::UnicycleState,t::Float64) where {T<:WaitTraj}
    get_action(controller.stabilizer,traj,state,t)
end
const PivotTraj = Union{PivotTrajectory,DenseTrajectory{PivotTrajectory}}
function get_action(controller::SwitchingController,traj::T,state::UnicycleState,t::Float64) where {T<:PivotTraj}
    get_action(controller.pivoter,traj,state,t)
end

################################################################################
################################## Simulation ##################################
################################################################################

function simulate(model::UnicycleModel,controller,traj,state,t0,tf,dt)
    states = [state]
    cmds = Vector{Vector{Float64}}()
    t = t0
    while t <= tf
        # get target state in global frame
        u           = get_action(controller,traj,state,t)
        # target_pt   = get_trajectory_point_by_time(traj, t)
        # u           = get_action(controller,target_pt,state,t)
        state     = forward_euler_integration(model,state,u,dt)
        push!(states, state)
        push!(cmds,u)
        t += dt
    end
    return states, cmds
end

function optimize_controller_params(sim_model,traj,initial_state;
        dt=0.1,
        k0_range=1.0:0.25:4.0,
        k1_range=1.0:0.25:4.0
    )
    t0 = get_start_time(traj)
    tf = get_end_time(traj)

    opt_vec = [NaN,NaN]
    opt_cost = Inf
    for k0 in k0_range
        for k1 in k1_range
            controller = TrackingController(k0=k0,k1=k1)
            states, cmds = simulate(sim_model,controller,traj,initial_state,t0,tf,dt)
            err = 0.0
            for (t,state) in zip(collect(t0:dt:tf),states)
                target_pt = get_trajectory_point_by_time(traj, t)
                target = [target_pt.pos.x, target_pt.pos.y, atan(target_pt.heading)]
                err += norm(target[1:2] - state[1:2])^2
            end
            if err < opt_cost
                opt_vec = [k0,k1]
            end
        end
    end
    opt_vec
end

end # end RobotModels
