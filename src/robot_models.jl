module RobotModels

using Parameters
using Vec, LinearAlgebra

using ..Trajectories
using ..GridPaths

abstract type AbstractRobotModel end
abstract type AbstractState end
abstract type AbstractAction end

export
    forward_euler_integration,
    TrackingController,
    sat,
    get_action,
    UnicycleModel,
    dynamics,
    simulate
    # UnicycleKinematicState,
    # UnicycleKinematicCommand,
    # next_state

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
"""
    `TrackingController`

    From Lee et al., "Tracking Control of
    Unicycle-Modeled Mobile Robots Using a Saturation Feedback Controller"
    [link]{https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=911382}
"""
@with_kw struct TrackingController
    dt  ::Float64 = 0.0
    a   ::Float64 = 4.0
    k0  ::Float64 = 1.0
    μ   ::Float64 = 1.0
    γ   ::Float64 = 0.1
    ϵ   ::Float64 = 0.1
    k1  ::Float64 = 1.0
    b   ::Float64 = 4.0
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
function get_action(controller::TrackingController,target,ff,state,t)
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
    if (abs(vr) <= 0.000001) && (abs(wr) <= 0.000001) # dirty trick to avoid problems with waiting
        return [0.0,0.0]
    end
    # errors in robot frame
    e = [ cos(θ) sin(θ) 0;
         -sin(θ) cos(θ) 0;
              0      0  1]*[xr-x,yr-y,get_angular_offset(θ,θr)]
    xe,ye,θe = e[1],e[2],e[3]
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
function get_action(controller::TrackingController,ref::TrajectoryPoint,state,t)
    target = [ref.pos.x, ref.pos.y, atan(ref.heading)]
    ff = [ref.yaw_rate, norm(ref.vel)]
    get_action(controller,target,ff,state,t)
end

struct UnicycleModel <: AbstractRobotModel end
function dynamics(model::UnicycleModel,state,cmd)
    x = state[1]
    y = state[2]
    θ = state[3]
    w = cmd[1]
    v = cmd[2]
    sdot = [v*cos(θ),v*sin(θ),w]
end

function simulate(model::UnicycleModel,controller::TrackingController,traj,state,t0,tf,dt)
    states = [state]
    cmds = Vector{Vector{Float64}}()
    t = t0
    while t <= tf
        # get target state in global frame
        target_pt = get_trajectory_point_by_time(traj, t)
        u         = get_action(controller,target_pt,state,t)
        state     = forward_euler_integration(model,state,u,dt)
        push!(states, state)
        push!(cmds,u)
        t += dt
    end
    return states, cmds
end

# struct UnicycleKinematicModel <: AbstractRobotModel
#     radius      ::Float64
#     track_width ::Float64
# end
# struct UnicycleKinematicState <: AbstractState
#     pos::VecE2
#     heading::VecE2
#     t::Float64
# end
# struct UnicycleKinematicCommand <: AbstractAction
#     left_wheel_rate::Float64
#     right_wheel_rate::Float64
# end
# function next_state(
#     m::UnicycleKinematicModel,
#     s::UnicycleKinematicState,
#     a::UnicycleKinematicCommand,
#     dt::Float64;
#     nsteps=1)
#
#     v_left = a.left_wheel_rate * m.radius
#     v_right = a.right_wheel_rate * m.radius
#     v = (v_left + v_right)/2
#     # r = (v_left + v_right)*(m.track_width/2) / (v_right - v_left)
#     k = (v_right - v_left) / ((v_left + v_right)*m.track_width/2)
#     yaw_rate = (v_right - v_left) / m.track_width
#     d = v*dt
#     θ = atan(s.heading)
#     Δxy = VecE2(0.0,0.0)
#     Δθ = d*k
#     for i in 1:nsteps
#         Δθ = yaw_rate * dt/nsteps
#         Δxy =
#         Δxy += (d/nsteps)*VecE2(cos(θ+Δθ/nsteps),sin(θ+Δθ/nsteps))
#         θ = θ + Δθ/nsteps
#     end
#     return UnicycleKinematicState(s.pos+Δxy,VecE2(cos(θ),sin(θ)),s.t+dt)
# end

# @with_kw struct TrackingController
#     dt  ::Float64 = 0.0
#     a   ::Float64 = 1.0
#     k0  ::Float64 = 1.0
#     μ   ::Float64 = 1.0
#     γ   ::Float64 = 0.1
#     ϵ   ::Float64 = 0.1
#     k1  ::Float64 = 1.0
#     b   ::Float64 = 1.0
# end

end # end RobotModels
