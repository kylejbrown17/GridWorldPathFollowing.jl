module RobotModels

using Vec, LinearAlgebra

using ..Trajectories

abstract type AbstractRobotModel end
abstract type AbstractState end
abstract type AbstractAction end

struct UnicycleRobotModel <: AbstractRobotModel
    wheel_radius::VecE2
    track_width::Float64
end
struct UnicycleState <: AbstractState
    pos::VecE2
    heading::VecE2
    vel::VecE2
    yaw_rate::Float64
    t::Float64
end
struct UnicycleCommand <: AbstractAction
    left_wheel_rate::Float64
    right_wheel_rate::Float64
end
function next_state(m::UnicycleRobotModel, s::UnicycleState, a::UnicycleCommand, dt::Float64)
    v_left = a.left_wheel_rate * m.radius
    v_right = a.right_wheel_rate * m.radius
    v = (v_left + v_right)/2
    r = (v_right - v_left)/m.track_width
    yaw_rate
end
