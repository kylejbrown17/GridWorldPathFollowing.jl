module Trajectories

using Parameters
using Vec, LinearAlgebra
using Convex, SCS, ECOS
using NearestNeighbors
using GraphUtils

using ..Utils
using ..GridPaths

export
    TimeInterval,
    verify,
    get_start_time,
    get_end_time,
    get_Δt,
    TrajectoryPoint,

    AbstractTrajectory,
    TrajectoryPrimitive,
    AbstractTrackingTrajectory,
    AbstractStabilizingTrajectory,
    AbstractPivotingTrajectory,
    CompositeTrajectory,
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
    verify_abstract_traj,

    get_time_from_arc_length,
    get_time_from_pt,
    get_trajectory_point_by_time,
    get_trajectory_point_by_pt,

    StraightTrajectory,
    ArcTrajectory,
    WaitTrajectory,
    PivotTrajectory,
    ReverseTrajectory,

    Trajectory,
    get_active_segment_idx,
    get_active_segment,
    cap_trajectory,
    push_with_reverse_flag!,
    construct_trajectory,
    optimize_velocity_profile,
    optimize_velocity_profile_traj_only,

    DenseTrajectory,
    get_index_time,
    SimpleTrajectory,

    stitch_trajectories

"""
    `TimeInterval`

    Defines the start and end time (the speed profile, essentially) for a
    particular AbstractTrajectory
"""
struct TimeInterval
    t1::Float64
    t2::Float64
end
get_start_time(interval::TimeInterval) = interval.t1
get_end_time(interval::TimeInterval) = interval.t2
function get_Δt(i::TimeInterval,t::Float64)
    dt = (i.t2 - i.t1)
    if dt == 0.0
        return 0.0
    end
    return (t - i.t1) / dt
 end
function verify(i::TimeInterval)
    @assert (
        get_start_time(i) <= get_end_time(i)
        ) "end time is before start time"
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
function Utils.interpolate(pt1::TrajectoryPoint,pt2::TrajectoryPoint,t::Float64)
    θ1 = atan(pt1.heading)
    θ2 = atan(pt2.heading)
    Δθ = get_angular_offset(θ1,θ2)
    θ = interpolate(θ1,θ1+Δθ,t)
    heading = VecE2(cos(θ),sin(θ))
    TrajectoryPoint(
        pos         = interpolate(pt1.pos, pt2.pos, t),
        heading     = heading,
        vel         = interpolate(pt1.vel, pt2.vel, t),
        yaw_rate    = interpolate(pt1.yaw_rate, pt2.yaw_rate, t),
        t           = interpolate(pt1.t, pt2.t, t)
    )
end
function verify(pt::TrajectoryPoint)
    @assert (!isnan(pt.t)) "time is NaN"
    @assert (!any(isnan,pt.pos)) "pos is NaN"
    @assert (!any(isnan,pt.heading)) "heading is NaN"
    @assert (!any(isnan,pt.vel)) "vel is NaN"
    @assert (!isnan(pt.yaw_rate)) "yaw_rate is NaN"
end

abstract type AbstractTrajectory end
abstract type TrajectoryPrimitive <: AbstractTrajectory end
abstract type AbstractTrackingTrajectory <: TrajectoryPrimitive end
abstract type AbstractStabilizingTrajectory <: TrajectoryPrimitive end
abstract type AbstractPivotingTrajectory <: TrajectoryPrimitive end
abstract type CompositeTrajectory <: AbstractTrajectory end
struct EmptyPrimitive <: AbstractTrajectory end
"""
    `get_start_time`

    Returns the start time of a given AbstractTrajectory.
"""
function get_start_time(traj::T) where {T <: TrajectoryPrimitive}
    return get_start_time(traj.interval)
end
"""
    `get_end_time`

    Returns the end time of a given AbstractTrajectory.
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
    AbstractTrajectory.
"""
function get_trajectory_point_by_time(traj::T,t::Float64) where {T <: AbstractTrajectory}
    pos         = get_position(traj,t)
    heading     = get_heading(traj,t)
    vel         = get_vel(traj,t)
    yaw_rate    = get_yaw_rate(traj,t)
    TrajectoryPoint(pos,heading,vel,yaw_rate,t)
end
"""
    `get_time_from_arc_length(traj::T,s::Float64)`

    Returns the time associated with a particular arc length
"""
function get_time_from_arc_length(traj::T,s::Float64) where {T <: AbstractTrajectory}
    Δt = get_end_time(traj) - get_start_time(traj)
    if get_length(traj) == 0.0
        return get_start_time(traj)
    end
    return get_start_time(traj) + Δt*(s / get_length(traj))
end
"""
    `get_trajectory_point_by_pt`

    Returns the TrajectoryPoint associated with a particular location along a
    AbstractTrajectory.
"""
function get_trajectory_point_by_pt(traj::T,pt::VecE2) where {T <: AbstractTrajectory}
    get_trajectory_point_by_time(traj,get_time_from_pt(traj,pt))
end
"""
    `get_closest_pt(traj::T, pt::VecE2) where {T <: AbstractTrajectory}`

    Finds the closest point on `traj` to `pt` (based on Euclidean distance)
"""
function get_closest_pt(traj::T, pt::VecE2) where {T <: AbstractTrajectory}
    t = get_time_from_pt(traj,pt)
    t = min(max(t,get_start_time(traj)),get_end_time(traj))
    get_trajectory_point_by_time(traj, t)
end
"""
    `get_start_pt`

    Returns the start point of a given AbstractTrajectory.
"""
function get_start_pt end
"""
    `get_end_pt`

    Returns the end point of a given AbstractTrajectory.
"""
function get_end_pt end
"""
    `get_length`

    Return the length of a AbstractTrajectory
"""
function get_length end
"""
    `get_dist`

    Return the distance along a AbstractTrajectory from the beginning up to a
    given time.
"""
function get_dist end
"""
    `get_position`

    Return the position associated with a particular time along a
    AbstractTrajectory
"""
function get_position end
"""
    `get_heading`

    Return the position associated with a particular time along a
    AbstractTrajectory
"""
function get_heading end
"""
    `get_vel`

    Return the speed associated with a particular time along a
    AbstractTrajectory
"""
function get_vel end
"""
    `get_yaw_rate`

    Return the yaw rate associated with a particular time along a
    AbstractTrajectory
"""
function get_yaw_rate end
"""
    `get_primitive`
"""
function get_primitive end
function verify_abstract_traj(traj::T) where {T<:AbstractTrajectory}
    t_mid = (get_start_time(traj)+get_end_time(traj))/2.0
    @assert !isnan(get_yaw_rate(traj,t_mid)) "get_yaw_rate is NaN"
    @assert !isnan(get_dist(traj,t_mid)) "get_dist is NaN"
    @assert !any(isnan,get_position(traj,t_mid)) "get_position is NaN"
    @assert !any(isnan,get_heading(traj,t_mid)) "get_heading is NaN"
    @assert !any(isnan,get_vel(traj,t_mid)) "get_vel is NaN"
end

"""
    `get_active_segment(traj::T,t::Float64) where {T <: AbstractTrajectory}`

    Returns the trajectory active at time `t`. For a base motion primitive (e.g.
    `StraightTrajectory`, `ArcTrajectory`, etc.) This is just `traj`.
    The function must be implemented for Composite Trajectory types such that
    a call to `get_active_segment` will recurse all the way down to the "leaf"
    trajectory primitives, thus permitting hierarchical definition of
    trajectories.
"""
get_active_segment(traj::T,t::Float64) where {T <: AbstractTrajectory} = traj

"""
    `StraightTrajectory`
"""
struct StraightTrajectory <: AbstractTrackingTrajectory
    pt1::VecE2 # start point
    # pt2::VecE2 # end point
    heading::VecE2
    length::Float64
    interval::TimeInterval
end
function StraightTrajectory(pt1::VecE2,pt2::VecE2,interval::TimeInterval)
    L = norm(pt2-pt1)
    StraightTrajectory(pt1,(pt2-pt1)/L,L,interval)
end
get_start_pt(traj::StraightTrajectory) = traj.pt1
get_end_pt(traj::StraightTrajectory) = traj.pt1 + traj.length * traj.heading
function verify(traj::StraightTrajectory)
    @assert (!any(isnan, get_start_pt(traj))) "start point is NaN"
    @assert (!any(isnan, get_end_pt(traj))) "end point is NaN"
    verify_abstract_traj(traj)
    verify(traj.interval)
end
function get_length(traj::StraightTrajectory)
    norm(get_end_pt(traj) - get_start_pt(traj))
end
function get_dist(traj::StraightTrajectory, t::Float64)
    dt = get_end_time(traj) - get_start_time(traj)
    return get_length(traj) * get_Δt(traj,t)
end
function get_position(traj::StraightTrajectory,t::Float64)
    p1 = get_start_pt(traj)
    p2 = get_end_pt(traj)
    return p1 + get_Δt(traj,t) * (p2 - p1)
end
function get_heading(traj::StraightTrajectory,t::Float64)
    # p1 = get_start_pt(traj)
    # p2 = get_end_pt(traj)
    # return (p2-p1) / norm(p2-p1)
    return traj.heading
end
function get_vel(traj::StraightTrajectory, t::Float64)
    # vec = get_end_pt(traj) - get_start_pt(traj)
    vec = get_heading(traj,t) * get_length(traj)
    dt = get_end_time(traj) - get_start_time(traj)
    if get_length(traj) == 0
        return VecE2(0.0,0.0)
    end
    return vec / dt
end
function get_yaw_rate(traj::StraightTrajectory,t::Float64)
    return 0.0
end
function get_time_from_pt(traj::StraightTrajectory,pt::VecE2)
    p1 = get_start_pt(traj)
    p2 = get_end_pt(traj)
    s = proj(pt-p1,get_heading(traj,get_end_time(traj)),Float64) * get_length(traj)
    get_time_from_arc_length(traj,s)
end

"""
    `ArcTrajectory`
"""
struct ArcTrajectory <: AbstractTrackingTrajectory
    center::VecE2
    radius::Float64 #
    heading1::VecE2
    θ1::Float64 # start angle of vector from center to pt1
    Δθ::Float64
    θ2::Float64 # end angle of vector from center to pt1
    interval::TimeInterval
end
function ArcTrajectory(center::VecE2,pt1::VecE2,heading1::VecE2,
    Δθ::Float64,interval::TimeInterval)
    radius = norm(pt1-center)
    sgn = sign(cross([(pt1-center)..., 0],[heading1..., 0 ])[end])
    if sgn == 0
        θ1 = atan(heading1) - π/2
    else
        θ1 = atan(heading1) - sgn*π/2
    end
    return ArcTrajectory(center,radius,heading1,mod(θ1,2π),Δθ,mod(θ1,2π)+Δθ,interval)
end
function ArcTrajectory(center::VecE2,radius::Float64,θ1::Float64,
    Δθ::Float64,interval::TimeInterval)
    if sign(Δθ) == 0
        heading1 = VecE2(NaN,NaN) # undefined
    else
        heading1 = VecE2(-sin(θ1),cos(θ1)) * sign(Δθ)
    end
    return ArcTrajectory(center,radius,heading1,mod(θ1,2π),Δθ,mod(θ1,2π)+Δθ,interval)
end
function verify(traj::ArcTrajectory)
    @assert !any(isnan, traj.center) "center is NaN"
    @assert (0.0 <= traj.θ1 <= 2π) "θ1 is outside of [0,2π]"
    # @assert (0.0 <= traj.θ2 <= 2π) "θ2 is outside of [0,2π]"
    @assert (traj.radius >= 0.0) "radius is negative"
    verify_abstract_traj(traj)
    verify(traj.interval)
end
get_start_pt(traj::ArcTrajectory) = traj.center + traj.radius*VecE2(cos(traj.θ1),sin(traj.θ1))
get_end_pt(traj::ArcTrajectory) = traj.center + traj.radius*VecE2(cos(traj.θ2),sin(traj.θ2))
function get_length(traj::ArcTrajectory)
    arclength = traj.radius * abs(traj.Δθ)
end
function get_dist(traj::ArcTrajectory, t::Float64)
    dt = get_end_time(traj) - get_start_time(traj)
    return get_length(traj) * get_Δt(traj,t) #(t-get_start_time(traj)) / dt
end
function get_position(traj::ArcTrajectory,t::Float64)
    θ = traj.θ1 + get_Δt(traj,t) * traj.Δθ
    return traj.center + traj.radius*VecE2(cos(θ),sin(θ))
end
function get_heading(traj::ArcTrajectory,t::Float64)
    Δθ = get_Δt(traj,t) * traj.Δθ
    heading = rot(traj.heading1,Δθ)
end
function get_vel(traj::ArcTrajectory, t::Float64)
    dt = get_end_time(traj) - get_start_time(traj)
    return get_heading(traj,t) * get_length(traj) / dt
end
function get_yaw_rate(traj::ArcTrajectory,t::Float64)
    dt = get_end_time(traj) - get_start_time(traj)
    return traj.Δθ / dt
end
function get_time_from_pt(traj::ArcTrajectory,pt::VecE2)
    if traj.Δθ == 0
        return get_start_time(traj)
    end
    θ = atan(pt-traj.center)
    Δθ = get_angular_offset(traj.θ1,θ)
    while sign(Δθ) != sign(traj.Δθ)
        Δθ += 2π*sign(traj.Δθ)
    end
    t = Δθ/traj.Δθ
end

"""
    ` WaitTrajectory <: AbstractTrajectory`

    When the robot should not be moving
"""
struct WaitTrajectory <: AbstractStabilizingTrajectory
    pt1::VecE2 # start point
    heading::VecE2
    interval::TimeInterval
end
function StraightTrajectory(traj::WaitTrajectory)
    StraightTrajectory(traj.pt1,traj.heading,0.0,traj.interval)
end
get_start_pt(traj::WaitTrajectory) = get_start_pt(StraightTrajectory(traj))
get_end_pt(traj::WaitTrajectory) = get_end_pt(StraightTrajectory(traj))
verify(traj::WaitTrajectory) = verify(StraightTrajectory(traj))
get_length(traj::WaitTrajectory) = 0.0
get_dist(traj::WaitTrajectory, t::Float64) = 0.0
get_position(traj::WaitTrajectory,t::Float64) = traj.pt1
get_heading(traj::WaitTrajectory,t::Float64) = traj.heading
get_vel(traj::WaitTrajectory, t::Float64) = VecE2(0.0,0.0)
get_yaw_rate(traj::WaitTrajectory,t::Float64) = get_yaw_rate(StraightTrajectory(traj),t)
get_time_from_pt(traj::WaitTrajectory,pt::VecE2) = get_time_from_pt(StraightTrajectory(traj),pt)

"""
    `PivotTrajectory`

    When the robot simply turns in place.
"""
struct PivotTrajectory <: AbstractPivotingTrajectory
    pt1::VecE2
    heading1::VecE2
    Δθ::Float64
    interval::TimeInterval
end
function ArcTrajectory(traj::PivotTrajectory)
    ArcTrajectory(traj.pt1,traj.pt1,traj.heading1,traj.Δθ,traj.interval)
end
get_start_pt(traj::PivotTrajectory) = get_start_pt(ArcTrajectory(traj))
get_end_pt(traj::PivotTrajectory) = get_end_pt(ArcTrajectory(traj))
verify(traj::PivotTrajectory) = verify(ArcTrajectory(traj))
get_length(traj::PivotTrajectory) = get_length(ArcTrajectory(traj))
get_dist(traj::PivotTrajectory, t::Float64) = get_dist(ArcTrajectory(traj), t)
get_position(traj::PivotTrajectory,t::Float64) = traj.pt1
get_heading(traj::PivotTrajectory,t::Float64) = get_heading(ArcTrajectory(traj),t)
get_vel(traj::PivotTrajectory, t::Float64) = VecE2(0.0,0.0)
get_yaw_rate(traj::PivotTrajectory,t::Float64) = get_yaw_rate(ArcTrajectory(traj),t)
get_time_from_pt(traj::PivotTrajectory,pt::VecE2) = get_time_from_pt(ArcTrajectory(traj),pt)

for T = (:StraightTrajectory, :ArcTrajectory, :WaitTrajectory, :PivotTrajectory)
    @eval get_primitive(traj::$T) = traj
end

export
    get_reversed_trajectory_point_by_time

"""
    `ReverseTrajectory`

    When the robot is supposed to be driving in reverse
"""
struct ReverseTrajectory{T} <: TrajectoryPrimitive
    traj::T
end
for op = (:get_start_time, :get_end_time, :get_start_pt, :get_end_pt, :verify, :get_length)
    @eval $op(traj::R) where {R<:ReverseTrajectory} = $op(traj.traj)
end
get_dist(traj::R,t::Float64) where {R<:ReverseTrajectory}      = get_dist(traj.traj, t)
get_position(traj::R,t::Float64) where {R<:ReverseTrajectory}   = get_position(traj.traj, t)
get_heading(traj::R,t::Float64) where {R<:ReverseTrajectory}    = get_heading(traj.traj, t)
get_vel(traj::R,t::Float64) where {R<:ReverseTrajectory}       = get_vel(traj.traj, t)
get_yaw_rate(traj::R,t::Float64) where {R<:ReverseTrajectory}   = get_yaw_rate(traj.traj, t)
get_time_from_pt(traj::R, pt::VecE2) where {R<:ReverseTrajectory} = get_time_from_pt(traj.traj, pt)
get_primitive(traj::R) where {R<:ReverseTrajectory} = get_primitive(traj.traj)
function get_reversed_trajectory_point_by_time(traj::R,t::Float64) where {R<:ReverseTrajectory}
    ref = get_trajectory_point_by_time(traj,t)
    TrajectoryPoint(
        pos     = ref.pos,
        heading = ref.heading,
        vel     = ref.vel,
        yaw_rate= ref.yaw_rate,
        t       = ref.t
    )
end

"""
    `DenseTrajectory`

    Contains an underlying `AbstractTrajectory` with extra information about the
    acceleration profile (and the resultant velocity and position profiles).
"""
struct DenseTrajectory{T<:AbstractTrajectory} <: CompositeTrajectory
    traj::T
    t_vec::Vector{Float64}
    accel::Vector{Float64}
    vel::Vector{Float64}
    pos::Vector{Float64}
end
for op = (:get_start_time, :get_end_time, :get_start_pt, :get_end_pt, :get_length)
    @eval $op(traj::R) where {R<:DenseTrajectory} = $op(traj.traj)
end
function verify(traj::DenseTrajectory)
    verify(traj.traj)
    @assert isapprox(traj.t_vec[1], get_start_time(traj.traj);atol=1e-10) "start times do not match"
    @assert isapprox(traj.t_vec[end], get_end_time(traj.traj);atol=1e-10) "end times do not match"
    @assert isapprox(traj.pos[end], get_length(traj.traj);atol=1e-10) "lengths do not match"
    @assert length(traj.vel) == length(traj.pos) == length(traj.t_vec) "t_vec, vel, pos have different lengths"
    @assert length(traj.accel) == max(0,length(traj.t_vec)-1) "vel and pos have different lengths"
    verify_abstract_traj(traj)
end
function get_dist(traj::DenseTrajectory, t::Float64)
    idx, Δe = linear_interp(traj.t_vec, t)
    # This could be replaced with precise integration of accel and vel signals...
    if traj.pos[idx] == traj.pos[idx+1]
        return traj.pos[idx]
    end
    s = interpolate(traj.pos[idx], traj.pos[idx+1], Δe)
    return s
end
"""
    `get_index_time(traj::DenseTrajectory, t::Float64)`
"""
function get_index_time(traj::DenseTrajectory, t::Float64)
    idx, Δe = linear_interp(traj.t_vec, t)
    ds = (traj.pos[idx+1] - traj.pos[idx])
    s = traj.pos[idx]+Δe*ds
    dt = (traj.t_vec[idx+1] - traj.t_vec[idx])
    if isapprox(ds,0.0;atol=1e-5)
        return traj.t_vec[idx] + Δe*dt
    end
    return get_time_from_arc_length(traj.traj, s)
end
function get_position(traj::DenseTrajectory, t::Float64)
    t_idx = get_index_time(traj,t)
    get_position(traj.traj, t_idx)
end
function get_heading(traj::DenseTrajectory, t::Float64)
    t_idx = get_index_time(traj,t)
    get_heading(traj.traj, t_idx)
end
function get_vel(traj::DenseTrajectory, t::Float64)
    idx, Δe = linear_interp(traj.t_vec, t)
    v = interpolate(traj.vel[idx], traj.vel[idx+1],Δe)
    return v * get_heading(traj,t)
end
function get_accel(traj::DenseTrajectory, t::Float64)
    idx, Δe = linear_interp(traj.t_vec, t)
    a = interpolate(traj.accel[idx], traj.accel[idx+1],Δe)
    return a
end
function get_yaw_rate(traj::DenseTrajectory,t::Float64)
    v_true = norm(get_vel(traj, t))
    t_idx = get_index_time(traj,t)
    v_nominal = norm(get_vel(traj.traj, t_idx))
    yw_nominal = get_yaw_rate(traj.traj, t_idx)
    if v_nominal != 0
        return yw_nominal * (v_true / v_nominal)
    else
        return yw_nominal
    end
end
get_primitive(traj::R) where {R<:DenseTrajectory} = get_primitive(traj.traj)

"""
    `Trajectory`
"""
@with_kw struct Trajectory <: CompositeTrajectory
    t_vec::Vector{Float64}               = Vector{Float64}() # time vector
    s_vec::Vector{Float64}               = Vector{Float64}() # arc length vector
    segments::Vector{AbstractTrajectory} = Vector{AbstractTrajectory}()
end
function Trajectory(segments::Vector{AbstractTrajectory})
    t_vec = map(seg->get_start_time(seg), segments)
    push!(t_vec, get_end_time(segments[end]))
    s_vec = cumsum([0.0, map(seg->get_length(seg), segments)...])
    Trajectory(t_vec,s_vec,segments)
end
Trajectory(v::Vector{T}) where {T <: AbstractTrajectory} = Trajectory(Vector{AbstractTrajectory}(v))
function Base.push!(traj::Trajectory,seg::T) where {T <: AbstractTrajectory}
    if length(traj.segments) == 0
        push!(traj.t_vec, get_start_time(seg))
        push!(traj.s_vec, 0.0)
    end
    push!(traj.t_vec, get_end_time(seg))
    push!(traj.s_vec, traj.s_vec[end] + get_length(seg))
    push!(traj.segments, seg)
    traj
end
function push_with_reverse_flag!(traj::Trajectory,seg::T,flag::Bool=false) where {T <: AbstractTrajectory}
    if flag
        push!(traj,ReverseTrajectory(seg))
    else
        push!(traj,seg)
    end
    traj
end
function verify(traj::Trajectory)
    for (i,p) in enumerate(traj.segments)
        if i < length(traj.segments)
            @assert (
                get_start_time(traj.segments[i+1]) == get_end_time(p) == traj.t_vec[i+1]
            ) "segment start and end times are not properly aligned"
            @assert (
                norm(get_start_pt(traj.segments[i+1]) - get_end_pt(p)) < 0.0000001
            ) "segment start and end positions are not properly aligned"
            @assert (
                length(traj.s_vec) == length(traj.t_vec)
            ) "s_vec and t_vec are not the same length"
        end
        verify(p)
    end
    verify_abstract_traj(traj)
end
get_length(traj::Trajectory) = sum([get_length(seg) for seg in traj.segments])
"""
    `get_active_segment_idx`

    Returns the integer index identifying which mode of the trajectory is active
"""
function get_active_segment_idx(traj::Trajectory,t::Float64)
    idx,Δe = linear_interp(traj.t_vec,t)
    return idx
end
function get_active_segment(traj::Trajectory,t::Float64)
    idx = get_active_segment_idx(traj,t)
    get_active_segment(traj.segments[idx],t)
end
function get_time_from_arc_length(traj::Trajectory,s::Float64)
    idx,Δe = linear_interp(traj.s_vec,s)
    seg = traj.segments[idx]
    get_time_from_arc_length(seg, s-traj.s_vec[idx])
end
for op = (:get_dist, :get_position, :get_heading, :get_vel, :get_yaw_rate)
    @eval $op(traj::T,t::Float64) where {T<:Trajectory} = $op(get_active_segment(traj,t),t)
end
for op = (:get_start_time, :get_start_pt)
    @eval begin
        function $op(traj::T) where {T<:Trajectory}
            if length(traj.segments) > 0
                return $op(traj.segments[1])
            else
                return NaN
            end
        end
    end
end
for op = (:get_end_time, :get_end_pt)
    @eval begin
        function $op(traj::T) where {T<:Trajectory}
            if length(traj.segments) > 0
                return $op(traj.segments[end])
            else
                return NaN
            end
        end
    end
end

"""
    `cap_trajectory(traj::Trajectory)`
"""
function cap_trajectory(traj::Trajectory)
    start = WaitTrajectory(
        get_start_pt(traj),
        get_heading(traj,get_start_time(traj)),
        TimeInterval(get_start_time(traj),get_start_time(traj))
        )
    finish = WaitTrajectory(
        get_end_pt(traj),
        get_heading(traj,get_end_time(traj)),
        TimeInterval(get_end_time(traj),get_end_time(traj))
        )
    new_traj = Trajectory()
    push!(new_traj,start)
    for seg in traj.segments
        push!(new_traj,seg)
    end
    push!(new_traj,finish)
    new_traj
end

"""
    `construct_trajectory(path::GridWorldPath;cap=true)`

    Returns a Trajectory that matches the sequence of vertices visited by `path`
    Inputs:
    - path::GridWorldPath - contains information about the initial starting
    position and time of the robot, as well as the series of transitions that
    the robot must make.
    - cap::Bool - when `true`, calls the `cap_trajectory` function on the traj
    before returning it.
"""
function construct_trajectory(path::GridWorldPath;cap::Bool=true)
    N = length(path.waypoints)
    @assert (N > 0) "path is empty"

    traj = Trajectory()
    # trajectory head: position, orientation, time
    pos     = path.start_pt
    ctr_pos = pos
    h       = VecE2(NaN,NaN)
    t       = path.start_time
    ctr_t   = t
    reverse_flag = false

    i = 1
    while i <= N
        w = path.waypoints[i]
        a = w.transition
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
            seg = StraightTrajectory(pos,get_heading_vector(a),path.cellwidth/2,TimeInterval(t,(t+ctr_t)/2))
            push_with_reverse_flag!(traj, seg, reverse_flag)
        elseif a == WAIT
            ctr_pos_0 = ctr_pos
            ctr_t_0 = ctr_t
            while i <= N
                w = path.waypoints[i]
                a = w.transition
                h_next = get_heading_vector(a)
                if a != WAIT
                    if abs(dot(h,h_next)) != 1 # abs(Δθ) > 0
                        radius = path.cellwidth / 2
                        center = pos + radius * h_next
                        Δθ = get_angular_offset(atan(h),atan(h_next))
                        # halfway t->ctr_t_0
                        seg = ArcTrajectory(center,get_end_pt(traj),get_heading(traj,get_end_time(traj)),Δθ/2,TimeInterval(t,ctr_t_0))
                        push_with_reverse_flag!(traj, seg, reverse_flag)
                        # wait ctr_t_0 -> ctr_t
                        # seg = ArcTrajectory(center,get_end_pt(traj),get_heading(traj,get_end_time(traj)),0.0,TimeInterval(ctr_t_0,ctr_t))
                        seg = WaitTrajectory(get_end_pt(traj),get_heading(traj,get_end_time(traj)),TimeInterval(ctr_t_0,ctr_t))
                        push_with_reverse_flag!(traj, seg, reverse_flag)
                        # halfway again ctr_t -> (ctr_t+w.t)/2
                        seg = ArcTrajectory(center,get_end_pt(traj),get_heading(traj,get_end_time(traj)),Δθ/2,TimeInterval(ctr_t,(ctr_t + w.t)/2))
                        push_with_reverse_flag!(traj, seg, reverse_flag)
                    else # straight and wait
                        # straight to center of cell
                        seg = StraightTrajectory(pos,h,path.cellwidth/2,TimeInterval(t,ctr_t_0))
                        push_with_reverse_flag!(traj, seg, reverse_flag)
                        # wait
                        # seg = StraightTrajectory(ctr_pos_0,h_next,0.0,TimeInterval(ctr_t_0,ctr_t))
                        seg = WaitTrajectory(get_end_pt(traj),get_heading(traj,get_end_time(traj)),TimeInterval(ctr_t_0,ctr_t))
                        push_with_reverse_flag!(traj, seg, reverse_flag)
                        # straight on
                        if dot(h,h_next) < 0
                            reverse_flag = !reverse_flag
                        end
                        seg = StraightTrajectory(ctr_pos,h_next,path.cellwidth/2,TimeInterval(ctr_t,(ctr_t+w.t)/2))
                        push_with_reverse_flag!(traj, seg, reverse_flag)
                    end
                    break
                else
                    ctr_pos = w.pt
                    ctr_t = w.t
                    # Handle the case where the robot is just waiting until the end
                    if i == N
                        # straight to center of cell
                        seg = StraightTrajectory(pos,h,path.cellwidth/2,TimeInterval(t,ctr_t_0))
                        push_with_reverse_flag!(traj, seg, reverse_flag)
                        # wait
                        # seg = WaitTrajectory(get_end_pt(traj),get_heading(traj,get_end_time(traj)),TimeInterval(ctr_t_0,ctr_t))
                        # push_with_reverse_flag!(traj, seg, reverse_flag)
                        if cap
                            return cap_trajectory(traj)
                        end
                        return traj
                    end
                    i += 1
                end
            end
        elseif dot(h,h_next) == 1 # abs(Δθ) == 0
            # straight to next boundary
            seg = StraightTrajectory(pos,h,path.cellwidth/2,TimeInterval(t,ctr_t))
            push_with_reverse_flag!(traj, seg, reverse_flag)
            seg = StraightTrajectory(ctr_pos,h_next,path.cellwidth/2,TimeInterval(ctr_t,(ctr_t + w.t)/2))
            push_with_reverse_flag!(traj, seg, reverse_flag)
        elseif dot(h,h_next) == 0 # abs(Δθ) == π/2
            # Turn to next boundary
            radius = path.cellwidth / 2
            center = pos + radius * h_next
            Δθ = get_angular_offset(atan(h),atan(h_next))
            seg = ArcTrajectory(center,pos,get_heading(traj,get_end_time(traj)),Δθ/2,TimeInterval(t,ctr_t))
            push_with_reverse_flag!(traj, seg, reverse_flag)
            seg = ArcTrajectory(center,pos,get_heading(traj,get_end_time(traj)),Δθ/2,TimeInterval(ctr_t,(ctr_t + w.t)/2))
            push_with_reverse_flag!(traj, seg, reverse_flag)
        elseif dot(h,h_next) == -1 # abs(Δθ) == π
            # reverse: straight to interior point, then backward to boundary
            # TODO this is still a problem!!! Need a flag telling the robot to
            # reverse its coordinate frame
            seg = StraightTrajectory(pos,ctr_pos,TimeInterval(t,ctr_t))
            push_with_reverse_flag!(traj, seg, reverse_flag)
            seg = WaitTrajectory(ctr_pos,h_next,TimeInterval(ctr_t,ctr_t))
            push_with_reverse_flag!(traj, seg, reverse_flag)
            # switch flag
            reverse_flag = !reverse_flag
            seg = StraightTrajectory(ctr_pos,pos,TimeInterval(ctr_t,(ctr_t + w.t)/2))
            push_with_reverse_flag!(traj, seg, reverse_flag)
        end
        pos = get_end_pt(traj)
        h   = get_heading_vector(a)
        t   = get_end_time(traj)
        ctr_pos = w.pt
        ctr_t   = w.t
        i += 1
    end
    # final half piece
    seg = StraightTrajectory(pos,get_heading(traj,get_end_time(traj)),
        path.cellwidth/2,TimeInterval(t,ctr_t))
    push_with_reverse_flag!(traj, seg, reverse_flag)
    if cap
        return cap_trajectory(traj)
    end
    return traj
end

struct SimpleTrajectory <: AbstractTrajectory
    pts ::Vector{TrajectoryPoint}
    t   ::Vector{Float64}
    s   ::Vector{Float64}
end
function SimpleTrajectory(traj::T,t_range) where {T <: AbstractTrajectory}
    pts     = Vector{TrajectoryPoint}()
    t_vec   = Vector{Float64}()
    s       = Vector{Float64}()
    for (i,t) in enumerate(t_range)
        push!(pts, get_trajectory_point_by_time(traj,t))
        push!(s, get_dist(traj,t))
        push!(t_vec, t)
    end
    SimpleTrajectory(pts,t_vec,s)
end
function SimpleTrajectory(traj::T,Δt::Float64) where{T <: AbstractTrajectory}
    t0 = get_start_time(traj)
    tf = get_end_time(traj)
    t_range = collect(t0:Δt:tf)
    SimpleTrajectory(traj,t_range)
end
function SimpleTrajectory(traj::T,N::Int) where{T <: AbstractTrajectory}
    t0 = get_start_time(traj)
    tf = get_end_time(traj)
    t_range = [t0 + (tf-t0)*(i/(N-1)) for i in 0:N-1]
    SimpleTrajectory(traj,t_range)
end
function verify(traj::SimpleTrajectory)
    @assert length(traj.t) == length(traj.s) == length(traj.pts) "lengths do not match"
    verify_abstract_traj(traj)
end
get_start_time(traj::SimpleTrajectory)  = traj.t[1]
get_end_time(traj::SimpleTrajectory)    = traj.t[end]
get_start_pt(traj::SimpleTrajectory)    = traj.pts[1]
get_end_pt(traj::SimpleTrajectory)      = traj.pts[end]
function get_trajectory_point_by_time(traj::SimpleTrajectory,t::Float64)
    idx,Δe = linear_interp(traj.t,t)
    pt = interpolate(traj.pts[idx],traj.pts[idx+1],Δe)
end
function get_dist(traj::SimpleTrajectory, t::Float64)
    idx, Δe = linear_interp(traj.t, t)
    interpolate(traj.s[idx],traj.s[idx+1], Δe)
end
function get_position(traj::SimpleTrajectory,t::Float64)
    get_trajectory_point_by_time(traj,t).pos
end
function get_vel(traj::SimpleTrajectory,t::Float64)
    get_trajectory_point_by_time(traj,t).vel
end
function get_heading(traj::SimpleTrajectory,t::Float64)
    get_trajectory_point_by_time(traj,t).heading
end
function get_yaw_rate(traj::SimpleTrajectory,t::Float64)
    get_trajectory_point_by_time(traj,t).yaw_rate
end

"""
    `optimize_velocity_profile(traj::Trajectory)`

    Given the constraints on robot position as a function of time, compute
    a dynamically feasible velocity profile that minimizes (the objective can
    be modified) the squared control effort.

    Params:
    * m::Int - the number of discrete control command blocks per trajectory.
    * n::Int - the polynomial order of the command signal within each command
    window

    Returns:
    - t_vec ::Vector{Float64} - a vector of time indices of length N*m + 1
    - accel ::Vector{Float64} - a vector of accel values of length N*m
    - vel   ::Vector{Float64} - a vector of speed values of length N*m + 1
    - pos   ::Vector{Float64} - a vector of position values of length N*m + 1

    TODO: make the optimization objective tunable via keyword parameters
"""
function optimize_velocity_profile(traj::Trajectory;
    m::Int=8,
    n::Int = 0,
    a_max::Float64 = 2.0,
    verbose=false
    )
    # # t = 0.5
    # # n = 3
    # # a = ones(n) # coefficients
    # # Polynomial basis vector
    # tv(t,n) = [t^j for j in 0:n]
    # # a(t) = a[1] + a[2]*t + a[3]*t^2 + ... + a[n]*t^(n-1)
    # #      = a'*tv
    # # Differentiation
    # Dmat(n) = vcat([[i*(i == j-1) for j in 1:n+1]' for i in 1:n+1]...)
    # # d/dt a(t) = a[2] + 2*a[3]*t + ... + (n-1)*a[n]*t^(n-2)
    # #           = (Dmat*a)'*tv
    # # Integration from 0 to t
    # Imat(n,d=1) = vcat([[(1/prod([k for k in j:j+d-1]))*(i == j+1) for j in 1:n+1]' for i in 2:n+2]...);
    # # ∫ a(t) = a[1]*t + (1/2)*a[2]*t^2 + ... + (1/n)*a[n]*t^n
    # #        = (Imat*a)'*t*tv

    N = length(traj.segments)
    stop_idxs = map(seg->norm(get_vel(seg,get_start_time(seg))) == 0, traj.segments) # true if v == 0 along this trajectory
    d = [get_length(seg) for seg in traj.segments]
    cd = cumsum(d)
    t = [[get_start_time(seg) for seg in traj.segments]..., get_end_time(traj)]
    dt = diff(t)
    v0 = 0.0
    vT = 0.0
    s0 = 0.0
    sT = get_length(traj)

    # n = 1 # order of acceleration polynomial: a[i] = a_0 + a_1*t + a_2*t^2 + ... a_n*t^n
    # m = 8 # num accel windows (m per seg, equally spaced in time)
    # each trajectory has a m-phase acceleration profile
    a = [Variable(m) for t in 1:N] # accel cmd coefficients
    v = [Variable(m+1) for t in 1:N]
    s = [Variable(m+1) for t in 1:N]

    constraints = [
        [a[i] <= a_max for i in 1:N]...,
        [a[i] >= -a_max for i in 1:N]...,
        v[1][1] == v0,    # initial conditions
        v[end][end] == vT,  # final conditions
        [v[i] >= 0.0 for i in 1:N]...,
        [v[i] == 0.0 for i in 1:N if stop_idxs[i] == true]...,
        [v[i][end] == v[i+1][1] for i in 1:N-1]...,
        [v[i][j+1] == v[i][j] + a[i][j]*(dt[i]/m) for i in 1:N for j in 1:m]...,# dynamics
        s[1][1] == s0,
        # s[end][end] == sT,
        [s[i][end] == cd[i] for i in 1:N]...,
        [s[i][end] == s[i+1][1] for i in 1:N-1]...,
        [s[i][j+1] == s[i][j] + v[i][j]*dt[i]/m + (1/2)*a[i][j]*(dt[i]/m)^2 for i in 1:N for j in 1:m]...,
    ]

    problem = minimize(sum([sumsquares(a[i]) for i in 1:N]),constraints);
    solve!(problem, ECOSSolver(verbose=verbose))
    # solve!(problem, SCSSolver(verbose=verbose))

    @assert !(problem.status in Set([:Infeasible,:Error,:Unbounded,:Indeterminate])) string("Optimization failed: problem.status == ",problem.status)

    new_traj = Trajectory()
    s0 = 0.0
    for (i,seg) in enumerate(traj.segments)
        t_vec = t[i] .+ (dt[i]/m)*collect(0:m)
        accel = a[i].value[:]
        vel = v[i].value[:]
        pos = s[i].value[:] .- s0
        push!(new_traj, DenseTrajectory(deepcopy(seg),t_vec,accel,vel,pos)) # using deepcopy to prevent the old pointers from lingering when embedding in c
        s0 += pos[end]
    end
    # return new_traj
    t_vec = [0, vcat([t[i] .+ (dt[i]/m)*collect(1:m) for i in 1:N]...)...]
    accel = vcat([a[i].value[:] for i in 1:N]...)
    vel = [vcat([v[i].value[1:end-1] for i in 1:N]...)..., v[end].value[end]]
    pos = [vcat([s[i].value[1:end-1] for i in 1:N]...)..., s[end].value[end]]

    return new_traj, t_vec, accel, vel, pos
end

function optimize_velocity_profile_traj_only(traj::Trajectory;
    m::Int=8,
    n::Int = 0,
    a_max::Float64 = 2.0,
    verbose=false
    )
    new_traj, t_vec, accel, vel, pos = optimize_velocity_profile(traj;
        m=m,n=n,a_max=a_max,verbose=verbose)
    return new_traj
end

"""
    `stitch_trajectories(traj1::Trajectory,traj2::Trajectory;buffer=[0.0,0.0])`

    Combines two trajectories with a `PivotTrajectory` sandwiched between two
    `WaitTrajectory` segments.
"""
function stitch_trajectories(traj1::Trajectory,traj2::Trajectory;buffer=[0.0,0.0])
    @assert norm(get_end_pt(traj1) - get_start_pt(traj2)) < 0.000001 "start and end points do not align"
    @assert buffer[1]+buffer[2] <= get_start_time(traj2) - get_end_time(traj1)
    θ1 = atan(get_heading(traj1,get_end_time(traj1)))
    θ2 = atan(get_heading(traj2,get_start_time(traj2)))
    Δθ = get_angular_offset(θ1,θ2)
    traj = Trajectory()
    for seg in traj1.segments
        push!(traj, seg)
    end
    push!(traj, WaitTrajectory(
        get_end_pt(traj),
        get_heading(traj,get_end_time(traj)),
        TimeInterval(get_end_time(traj),get_end_time(traj)+buffer[1]))
        )
    push!(traj, PivotTrajectory(
        get_end_pt(traj),
        get_heading(traj,get_end_time(traj)),
        Δθ,TimeInterval(get_end_time(traj),get_start_time(traj2)-buffer[2]))
        )
    push!(traj, WaitTrajectory(
        get_end_pt(traj),
        get_heading(traj,get_end_time(traj)),
        TimeInterval(get_end_time(traj),get_start_time(traj2)))
        )
    for seg in traj2.segments
        push!(traj, seg)
    end
    return traj
end

end # module Trajectories
