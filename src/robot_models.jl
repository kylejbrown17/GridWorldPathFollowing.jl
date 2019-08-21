module RobotModels

using Vec, LinearAlgebra

using ..Trajectories

abstract type AbstractRobotModel end
abstract type AbstractState end
abstract type AbstractAction end

export
    UnicycleKinematicModel,
    UnicycleKinematicState,
    UnicycleKinematicCommand,
    next_state

struct UnicycleKinematicModel <: AbstractRobotModel
    radius      ::Float64
    track_width ::Float64
end
struct UnicycleKinematicState <: AbstractState
    pos::VecE2
    heading::VecE2
    t::Float64
end
struct UnicycleKinematicCommand <: AbstractAction
    left_wheel_rate::Float64
    right_wheel_rate::Float64
end
function next_state(
    m::UnicycleKinematicModel,
    s::UnicycleKinematicState,
    a::UnicycleKinematicCommand,
    dt::Float64;
    nsteps=1)

    v_left = a.left_wheel_rate * m.radius
    v_right = a.right_wheel_rate * m.radius
    v = (v_left + v_right)/2
    # r = (v_left + v_right)*(m.track_width/2) / (v_right - v_left)
    k = (v_right - v_left) / ((v_left + v_right)*m.track_width/2)
    yaw_rate = (v_right - v_left) / m.track_width
    d = v*dt
    θ = atan(s.heading)
    Δxy = VecE2(0.0,0.0)
    Δθ = d*k
    for i in 1:nsteps
        Δθ = yaw_rate * dt/nsteps
        Δxy += (d/nsteps)*VecE2(cos(θ+Δθ/nsteps),sin(θ+Δθ/nsteps))
        θ = θ + Δθ/nsteps
    end
    return UnicycleKinematicState(s.pos+Δxy,VecE2(cos(θ),sin(θ)),s.t+dt)
end

end # end RobotModels
