module Trajectories

using Parameters
using Vec, LinearAlgebra
using ..GridPaths

export
    TimeInterval,
    verify,
    get_start_time,
    get_end_time,
    get_Δt,
    TrajectoryPoint,
    TrajectoryPrimitive,
    EmptyPrimitive,
    get_Δt,
    get_start_pt,
    get_end_pt,
    get_position,
    get_heading,
    get_vel,
    get_yaw_rate,
    get_time_from_pt,
    get_trajectory_point_by_time,
    get_trajectory_point_by_pt,

    ConstSpeedStraightTrajectory,
    ConstSpeedArcTrajectory,

    Trajectory,
    get_active_segment_idx,
    construct_trajectory

"""
    `TimeInterval`

    Defines the start and end time (the speed profile, essentially) for a
    particular TrajectoryPrimitive
"""
struct TimeInterval
    t1::Float64
    t2::Float64
end
get_start_time(interval::TimeInterval) = interval.t1
get_end_time(interval::TimeInterval) = interval.t2
get_Δt(i::TimeInterval,t::Float64) = (t - i.t1) / (i.t2 - i.t1)
function verify(i::TimeInterval)
    @assert (
        get_start_time(i) < get_end_time(i)
        ) "end time is before or coincident with start time"
end
"""
    `TrajectoryPoint`

    The feed forward command (2D rigid body state) associated with a particular
    point on a trajectory.
"""
@with_kw struct TrajectoryPoint
    pos::VecE2          = VecE2()
    heading::VecE2      = VecE2()
    vel::VecE2          = VecE2()
    yaw_rate::Float64   = 0.0
    t::Float64          = NaN
end
function verify(pt::TrajectoryPoint)
    @assert (!isnan(pt.t)) "time is NaN"
    @assert (!any(isnan,pt.pos)) "pos is NaN"
    @assert (!any(isnan,pt.heading)) "heading is NaN"
    @assert (!any(isnan,pt.vel)) "vel is NaN"
    @assert (!isnan(pt.yaw_rate)) "yaw_rate is NaN"
end
abstract type TrajectoryPrimitive end
struct EmptyPrimitive <: TrajectoryPrimitive end
"""
    `get_start_time`

    Returns the start time of a given TrajectoryPrimitive.
"""
function get_start_time(traj::T) where {T <: TrajectoryPrimitive}
    return get_start_time(traj.interval)
end
"""
    `get_end_time`

    Returns the end time of a given TrajectoryPrimitive.
"""
function get_end_time(traj::T) where {T <: TrajectoryPrimitive}
    return get_end_time(traj.interval)
end
"""
    `get_Δt`

    Returns the fraction Δt by which t interpolates `get_start_time(traj)` and
    `get_end_time(traj)`
"""
function get_Δt(traj::T,t::Float64) where {T <: TrajectoryPrimitive}
    get_Δt(traj.interval,t)
end
"""
    `get_trajectory_point_by_time`

    Returns the TrajectoryPoint associated with a particular time along a
    TrajectoryPrimitive.
"""
function get_trajectory_point_by_time(traj::T,t::Float64) where {T <: TrajectoryPrimitive}
    pos         = get_position(traj,t)
    heading     = get_heading(traj,t)
    vel         = get_vel(traj,t)
    yaw_rate    = get_yaw_rate(traj,t)
    TrajectoryPoint(pos,heading,vel,yaw_rate,t)
end
"""
    `get_trajectory_point_by_pt`

    Returns the TrajectoryPoint associated with a particular location along a
    TrajectoryPrimitive.
"""
function get_trajectory_point_by_pt(traj::T,pt::VecE2) where {T <: TrajectoryPrimitive}
    get_trajectory_point_by_time(traj,get_time_from_pt(traj,pt))
end
"""
    `get_start_pt`

    Returns the start point of a given TrajectoryPrimitive.
"""
function get_start_pt end
"""
    `get_end_pt`

    Returns the end point of a given TrajectoryPrimitive.
"""
function get_end_pt end
"""
    `get_position`

    Return the position associated with a particular time along a
    TrajectoryPrimitive
"""
function get_position end
"""
    `get_heading`

    Return the position associated with a particular time along a
    TrajectoryPrimitive
"""
function get_heading end
"""
    `get_vel`

    Return the speed associated with a particular time along a
    TrajectoryPrimitive
"""
function get_vel end
"""
    `get_yaw_rate`

    Return the yaw rate associated with a particular time along a
    TrajectoryPrimitive
"""
function get_yaw_rate end

"""
    `ConstSpeedStraightTrajectory`
"""
struct ConstSpeedStraightTrajectory <: TrajectoryPrimitive
    pt1::VecE2 # start point
    pt2::VecE2 # end point
    interval::TimeInterval
end
get_start_pt(traj::ConstSpeedStraightTrajectory) = traj.pt1
get_end_pt(traj::ConstSpeedStraightTrajectory) = traj.pt2
function verify(traj::ConstSpeedStraightTrajectory)
    @assert (!any(isnan, get_start_pt(traj))) "start point is NaN"
    @assert (!any(isnan, get_end_pt(traj))) "end point is NaN"
    verify(traj.interval)
end
function get_position(traj::ConstSpeedStraightTrajectory,t::Float64)
    p1 = get_start_pt(traj)
    p2 = get_end_pt(traj)
    return p1 + get_Δt(traj,t) * (p2 - p1)
end
function get_heading(traj::ConstSpeedStraightTrajectory,t::Float64)
    p1 = get_start_pt(traj)
    p2 = get_end_pt(traj)
    return (p2-p1) / norm(p2-p1)
end
function get_vel(traj::ConstSpeedStraightTrajectory, t::Float64)
    vec = get_end_pt(traj) - get_start_pt(traj)
    dt = get_end_time(traj) - get_start_time(traj)
    return vec / dt
end
function get_yaw_rate(traj::ConstSpeedStraightTrajectory,t::Float64)
    return 0.0
end
function get_time_from_pt(traj::ConstSpeedStraightTrajectory,pt::VecE2)
    p1 = get_start_pt(traj)
    p2 = get_end_pt(traj)
    t  = proj(pt-p1,p2-p1,Float64) + get_start_time(traj)
end

"""
    `ConstSpeedArcTrajectory`
"""
struct ConstSpeedArcTrajectory <: TrajectoryPrimitive
    center::VecE2
    radius::Float64 #
    θ1::Float64 # start angle
    θ2::Float64 # end angle
    interval::TimeInterval
end
function verify(traj::ConstSpeedArcTrajectory)
    @assert !any(isnan, traj.center) "center is NaN"
    @assert (0.0 <= traj.θ1 <= 2π) "θ1 is outside of [0,2π]"
    @assert (0.0 <= traj.θ2 <= 2π) "θ2 is outside of [0,2π]"
    @assert (traj.radius >= 0.0) "radius is negative"
    verify(traj.interval)
end
get_start_pt(traj::ConstSpeedArcTrajectory) = traj.center + traj.radius*VecE2(cos(traj.θ1),sin(traj.θ1))
get_end_pt(traj::ConstSpeedArcTrajectory) = traj.center + traj.radius*VecE2(cos(traj.θ2),sin(traj.θ2))
function get_position(traj::ConstSpeedArcTrajectory,t::Float64)
    θ = get_Δt(traj,t) * (traj.θ2 - traj.θ1)
    return traj.center + traj.radius*VecE2(cos(θ),sin(θ))
end
function get_heading(traj::ConstSpeedArcTrajectory,t::Float64)
    θ = get_Δt(traj,t) * (traj.θ2 - traj.θ1)
    return VecE2(-sin(θ),cos(θ))
end
function get_vel(traj::ConstSpeedArcTrajectory, t::Float64)
    arclength = traj.radius * abs(traj.θ2-traj.θ1)
    dt = get_end_time(traj) - get_start_time(traj)
    return get_heading(traj,t) * arclength / dt
end
function get_yaw_rate(traj::ConstSpeedArcTrajectory,t::Float64)
    dt = get_end_time(traj) - get_start_time(traj)
    return (traj.θ2-traj.θ1) / dt
end
function get_time_from_pt(traj::ConstSpeedArcTrajectory,pt::VecE2)
    θ1 = traj.θ1
    θ2 = traj.θ2
    θ = atan(pt-traj.center)
    t = (θ-θ1)/(θ2-θ1)
end

"""
    `Trajectory`
"""
@with_kw struct Trajectory
    segments::Vector{TrajectoryPrimitive} = Vector{TrajectoryPrimitive}()
end
function verify(traj::Trajectory)
    for (i,p) in enumerate(traj.segments)
        if i < length(traj.segments)
            @assert (
                get_start_time(traj.segments[i+1]) == get_end_time(p)
            ) "segment start and end times are not properly aligned"
        end
        verify(p)
    end
end
function construct_trajectory(path::GridWorldPath)
    pt = path.start_pt
    t = path.start_time
    traj = Trajectory()

    i = 1
    while i < length(path.waypoints)
        w1 = path.waypoints[i]
        w2 = path.waypoints[i+1]

        if w1.transition == w2.transition
            seg = ConstSpeedStraightTrajectory(pt,w1.pt,TimeInterval(t,w1.t))
            pt = get_end_pt(seg)
            t = get_end_time(seg)
            i += 1
        else
            # Turn
            center = (pt + w2.pt)/2
            radius = path.cellwidth / (2*sqrt(2))
            # WAIT
            # REVERSE
        end
    end
    # for i in 1:length(path.waypoints)-1
    #     w1 = path.waypoints[i]
    #     w2 = path.waypoints[i+1]
    #     # initial straight segment
    #
    #     # curved segment or two straights
    #
    #     if w1.transition == w2.transition
    #         seg = ConstSpeedStraightTrajectory(pt,w1.pt,TimeInterval(t,w1.t))
    #     end
    # end
end
function get_start_time(traj::Trajectory)
    if length(traj.segments) > 0
        return get_start_time(traj.segments[1])
    else
        return NaN
    end
end
function get_end_time(traj::Trajectory)
    if length(traj.segments) > 0
        return get_end_time(traj.segments[end])
    else
        return NaN
    end
end
"""
    `get_active_segment_idx`

    Returns the integer index identifying which mode of the trajectory is active
"""
function get_active_segment_idx(traj::Trajectory,t::Float64)
    for (i,p) in enumerate(traj.segments)
        if get_start_time(p) > t
            break
        end
        if get_end_time(p) <= t
            continue
        end
        return i
    end
    return -1
end
get_position(traj::Trajectory,t::Float64,idx::Int)  = get_position(traj.segments[idx],t)
get_heading(traj::Trajectory,t::Float64,idx::Int)   = get_heading(traj.segments[idx],t)
get_vel(traj::Trajectory,t::Float64,idx::Int)       = get_vel(traj.segments[idx],t)
get_yaw_rate(traj::Trajectory,t::Float64,idx::Int)  = get_yaw_rate(traj.segments[idx],t)
function get_position(traj::Trajectory,t::Float64)
    idx = get_active_segment_idx(traj::Trajectory,t::Float64)
    if idx != -1
        return get_position(traj,t,idx)
    end
    return VecE2(NaN,NaN)
end
function get_heading(traj::Trajectory,t::Float64)
    idx = get_active_segment_idx(traj::Trajectory,t::Float64)
    if idx != -1
        return get_heading(traj,t,idx)
    end
    return VecE2(NaN,NaN)
end
function get_vel(traj::Trajectory,t::Float64)
    idx = get_active_segment_idx(traj::Trajectory,t::Float64)
    if idx != -1
        return get_vel(traj,t,idx)
    end
    return VecE2(NaN,NaN)
end
function get_yaw_rate(traj::Trajectory,t::Float64)
    idx = get_active_segment_idx(traj::Trajectory,t::Float64)
    if idx != -1
        return get_yaw_rate(traj,t,idx)
    end
    return NaN
end

end
