module GridPaths

using Vec

export
    GridTransition,
        LEFT,RIGHT,UP,DOWN,WAIT,
    get_translation,
    get_angular_offset,
    GridWaypoint,
    GridWorldPath,
    construct_grid_world_path

@enum GridTransition begin
    RIGHT   = 1
    UP      = 2
    LEFT    = 3
    DOWN    = 4
    WAIT    = 5
end
"""
    `get_translation`

    Returns the direction associated with a particular `GridTransition`
"""
function get_translation(a::GridTransition)
    if a == LEFT
        return VecE2(-1.0, 0.0)
    elseif a == RIGHT
        return VecE2(1.0, 0.0)
    elseif a == UP
        return VecE2(0.0,1.0)
    elseif a == DOWN
        return VecE2(0.0,-1.0)
    elseif a == WAIT
        return VecE2(0.0,0.0)
    end
    return VecE2(0.0,0.0)
end
"""
    `get_angular_offset`

    Returns the angular offset between two `GridTransition`s
"""
function get_angular_offset(a1::GridTransition,a2::GridTransition)
    if (a1 == WAIT || a2 == WAIT)
        return 0.0
    else
        d = Int(a2)-Int(a1)
        if d > 2
            d = d - 4
        elseif d < -2
            d = d + 4
        end
        return d * π/2
    end
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
struct GridWorldPath
    start_pt::VecE2
    start_time::Float64
    waypoints::Vector{GridWaypoint}
    # environment parameters
    cellwidth::Float64
end
function construct_grid_world_path(
    start_pt::VecE2,
    start_time::Float64,
    transitions::Vector{GridTransition},
    cellwidth::Float64,
    transition_time::Float64=1.0)
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
    GridWorldPath(start_pt, start_time, waypoints, cellwidth)
end

end # End of module Paths
