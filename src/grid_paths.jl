module GridPaths

using Vec
using Parameters

export
    GridTransition,
        WEST,EAST,NORTH,SOUTH,WAIT,
    get_translation,
    get_heading_angle,
    get_heading_vector,
    get_angular_offset,
    GridWaypoint,
    GridWorldPath,
    construct_grid_world_path

@enum GridTransition begin
    EAST    = 0
    NORTH   = 1
    WEST    = 2
    SOUTH   = 3
    WAIT    = 4
end
"""
    `get_translation`

    Returns the direction associated with a particular `GridTransition`
"""
function get_translation(a::GridTransition)
    if a == WEST
        return VecE2(-1.0, 0.0)
    elseif a == EAST
        return VecE2(1.0, 0.0)
    elseif a == NORTH
        return VecE2(0.0,1.0)
    elseif a == SOUTH
        return VecE2(0.0,-1.0)
    elseif a == WAIT
        return VecE2(0.0,0.0)
    end
    return VecE2(0.0,0.0)
end
"""
    `get_heading_angle`

    Returns the heading angle corresponding to a `GridTransition`s
"""
function get_heading_angle(a::GridTransition)
    if a == WAIT
        return 0.0
    elseif a == EAST
        return 0.0
    elseif a == NORTH
        return π/2
    elseif a == WEST
        return π
    elseif a == SOUTH
        return 3π/2
    end
    return NaN
end
"""
    `get_heading_vector`

    Returns the heading vector corresponding to a `GridTransition`s
"""
function get_heading_vector(a::GridTransition)
    if a == WAIT
        return VecE2(0.0,0.0)
    elseif a == EAST
        return VecE2(1.0,0.0)
    elseif a == NORTH
        return VecE2(0.0,1.0)
    elseif a == WEST
        return VecE2(-1.0,0.0)
    elseif a == SOUTH
        return VecE2(0.0,-1.0)
    end
    return VecE2(NaN,NaN)
end
"""
    `get_angular_offset`

    Returns the angular offset between two `GridTransition`s
"""
function get_angular_offset(θ1,θ2)
    Δθ = θ2 - θ1
    if Δθ > π
        Δθ = Δθ - 2π
    elseif Δθ < -π
        Δθ = Δθ + 2π
    end
    return Δθ
end

struct GridWaypoint
    pt::VecE2
    t::Float64
    transition::GridTransition
end
GridWaypoint(x::Float64,y::Float64,t::Float64,a::GridTransition) = GridWaypoint(VecE2(x,y),t,a)

"""
    `GridWorldPath`
"""
@with_kw struct GridWorldPath
    start_pt::VecE2                 = VecE2()
    start_time::Float64             = 0.0
    waypoints::Vector{GridWaypoint} = Vector{GridWaypoint}()
    # environment parameters
    cellwidth::Float64              = 1.0
end
function construct_grid_world_path(
    start_pt::VecE2,
    start_time::Float64,
    cellwidth::Float64,
    transition_time::Float64,
    transitions::Vector{GridTransition})
    waypoints = Vector{GridWaypoint}()
    pt = start_pt
    t = start_time
    for a in transitions
        next_pt = pt + cellwidth * get_translation(a)
        next_t = t + transition_time
        push!(waypoints, GridWaypoint(next_pt, next_t, a))
        pt = next_pt
        t = next_t
    end
    GridWorldPath(
        start_pt=start_pt,
        start_time=start_time,
        waypoints=waypoints,
        cellwidth=cellwidth)
end
function construct_grid_world_path(
    start_pt::VecE2,
    start_time::Float64,
    transitions::Vector{GridTransition},
    cellwidth::Float64,
    transition_time::Float64=1.0)
    construct_grid_world_path(start_pt,start_time,cellwidth,transition_time,transitions)
end
function construct_grid_world_path(
    x0::Float64,
    y0::Float64,
    t0::Float64,
    cellwidth::Float64,
    transition_time::Float64,
    transitions::Vector{GridTransition}
    )
    construct_grid_world_path(
        VecE2(x0,y0),
        t0,
        cellwidth,
        transition_time,
        transitions
        )
end
function construct_grid_world_path(
    x0::Float64,
    y0::Float64,
    t0::Float64,
    cellwidth::Float64,
    transition_time::Float64,
    transitions::Vector{Int}
    )
    construct_grid_world_path(
        x0, y0, t0, cellwidth, transition_time,
        map(i->GridTransition(i),transitions)
        )
end

end # End of module Paths
