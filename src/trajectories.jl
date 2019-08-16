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
    get_length,
    get_dist,
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
    `get_length`

    Return the length of a TrajectoryPrimitive
"""
function get_length end
"""
    `get_dist`

    Return the distance along a TrajectoryPrimitive from the beginning up to a
    given time.
"""
function get_dist end
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
function get_length(traj::ConstSpeedStraightTrajectory)
    norm(get_end_pt(traj) - get_start_pt(traj))
end
function get_dist(traj::ConstSpeedStraightTrajectory, t::Float64)
    dt = get_end_time(traj) - get_start_time(traj)
    return get_length(traj)*(t-get_start_time(traj)) / dt
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
    θ1::Float64 # start angle of vector from center to pt1
    Δθ::Float64
    θ2::Float64 # end angle of vector from center to pt1
    interval::TimeInterval
end
function ConstSpeedArcTrajectory(center::VecE2,radius::Float64,θ1::Float64,
    Δθ::Float64,interval::TimeInterval)
    return ConstSpeedArcTrajectory(center,radius,θ1,Δθ,mod(θ1+Δθ,2π),interval)
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
function get_length(traj::ConstSpeedArcTrajectory)
    arclength = traj.radius * abs(traj.Δθ)
end
function get_dist(traj::ConstSpeedArcTrajectory, t::Float64)
    dt = get_end_time(traj) - get_start_time(traj)
    return get_length(traj) * (t-get_start_time(traj)) / dt
end
function get_position(traj::ConstSpeedArcTrajectory,t::Float64)
    θ = traj.θ1 + get_Δt(traj,t) * traj.Δθ
    return traj.center + traj.radius*VecE2(cos(θ),sin(θ))
end
function get_heading(traj::ConstSpeedArcTrajectory,t::Float64)
    θ = traj.θ1 + get_Δt(traj,t) * traj.Δθ
    return VecE2(-sin(θ),cos(θ))
end
function get_vel(traj::ConstSpeedArcTrajectory, t::Float64)
    dt = get_end_time(traj) - get_start_time(traj)
    return get_heading(traj,t) * get_length(traj) / dt
end
function get_yaw_rate(traj::ConstSpeedArcTrajectory,t::Float64)
    dt = get_end_time(traj) - get_start_time(traj)
    return traj.Δθ / dt
end
function get_time_from_pt(traj::ConstSpeedArcTrajectory,pt::VecE2)
    θ = atan(pt-traj.center)
    t = abs(get_angular_offset(traj.θ1,θ))/abs(traj.Δθ)
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
function get_start_pt(traj::Trajectory)
    if length(traj.segments) > 0
        return get_start_pt(traj.segments[1])
    else
        return NaN
    end
end
function get_end_pt(traj::Trajectory)
    if length(traj.segments) > 0
        return get_end_pt(traj.segments[end])
    else
        return NaN
    end
end
get_length(traj::Trajectory) = sum([get_length(seg) for seg in traj.segments])
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
function get_dist(traj::Trajectory,t::Float64)
    idx = get_active_segment_idx(traj::Trajectory,t::Float64)
    if idx != -1
        d = 0.0
        for i in 1:idx-1
            d += get_length(traj.segments[i])
        end
        d += get_dist(traj.segments[idx],t)
        return d
    elseif t <= get_start_time(traj)
        return 0.0
    elseif t >= get_end_time(traj)
        return get_length(traj)
    end
    return NaN
end
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
"""
    `construct_trajectory`

    Returns a Trajectory that matches the sequence of vertices visited by `path`
    Inputs:
    - beginning robot state (`TrajectoryPoint`)
    - instructions (`GridWorldPath`)
"""
function construct_trajectory(path::GridWorldPath)#,traj_pt::TrajectoryPoint;ϵ = 0.01)
    N = length(path.waypoints)
    @assert (N > 0) "path is empty"

    traj = Trajectory()
    # trajectory head: position, orientation, time
    pos     = path.start_pt
    ctr_pos = pos
    h       = VecE2(NaN,NaN)
    t       = path.start_time
    ctr_t   = t

    i = 1
    while i <= N
        w      = path.waypoints[i]
        a      = w.transition
        h_next = get_heading_vector(a)
        if i == 1
            # Stay in place if a == WAIT
            while (a == WAIT) && (i < N)
                # increment head (time only)
                t = w.t
                # skip to next action
                i += 1
                w = path.waypoints[i]
                a = w.transition
            end
            ctr_pos = w.pt
            ctr_t   = w.t
            # straight to boundary between cells
            seg = ConstSpeedStraightTrajectory(pos,(pos+ctr_pos)/2,TimeInterval(t,(t+ctr_t)/2))
            push!(traj.segments, seg)
        elseif a == WAIT
            seg = ConstSpeedStraightTrajectory(pos,ctr_pos,TimeInterval(t,ctr_t))
            while i <= N
                w = path.waypoints[i]
                a = w.transition
                if a != WAIT
                    # straight to center of cell
                    push!(traj.segments, seg)
                    pos = get_end_pt(traj)
                    t   = get_end_time(traj)
                    h_next = get_heading_vector(a)
                    # turn in place if necessary
                    if abs(dot(h,h_next) != 1) # abs(Δθ) > 0
                        Δθ = get_angular_offset(atan(h),atan(h_next))
                        seg = ConstSpeedArcTrajectory(
                            pos,0.0,atan(h),Δθ,atan(h_next),TimeInterval(t,ctr_t))
                        pos = get_end_pt(traj)
                        t   = get_end_time(traj)
                    end
                    ctr_pos = w.pt
                    ctr_t = w.t
                    seg = ConstSpeedStraightTrajectory(pos,(pos+ctr_pos)/2,TimeInterval(t,(t+ctr_t)/2))
                    push!(traj.segments, seg)
                    break
                else
                    ctr_pos = w.pt
                    ctr_t = w.t
                    i += 1
                end
            end
        elseif dot(h,h_next) == 1 # abs(Δθ) == 0
            # straight to next boundary
            seg = ConstSpeedStraightTrajectory(pos,(ctr_pos + w.pt)/2,TimeInterval(t,(ctr_t + w.t)/2))
            push!(traj.segments, seg)
        elseif dot(h,h_next) == 0 # abs(Δθ) == π/2
            # Turn to next boundary
            radius = path.cellwidth / 2
            center = pos + radius * h_next
            θ1 = mod(atan(h_next) + π, 2π)
            Δθ = get_angular_offset(θ1,atan(h))
            θ2 = θ1 + Δθ
            seg = ConstSpeedArcTrajectory(center,radius,θ1,Δθ,θ2,TimeInterval(t,(ctr_t + w.t)/2))
            push!(traj.segments, seg)
        elseif dot(h,h_next) == -1 # abs(Δθ) == π
            # reverse: straight to interior point, then backward to boundary
            seg = ConstSpeedStraightTrajectory(pos,ctr_pos,TimeInterval(t,ctr_t))
            push!(traj.segments, seg)
            seg = ConstSpeedStraightTrajectory(ctr_pos,pos,TimeInterval(ctr_t,(ctr_t + w.t)/2))
            push!(traj.segments, seg)
        end
        # push!(traj.segments, seg)
        pos = get_end_pt(traj)
        h   = get_heading_vector(a)
        t   = get_end_time(traj)
        ctr_pos = w.pt
        ctr_t   = w.t
        i += 1
    end
    # final half piece
    seg = ConstSpeedStraightTrajectory(pos,ctr_pos,TimeInterval(t,ctr_t))
    push!(traj.segments, seg)
    return traj
end

end
