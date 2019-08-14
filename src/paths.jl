module Paths

export
    GridTransition,
    GridWaypoint,
    GridWorldPath

@enum GridTransition begin
    LEFT        = 1
    RIGHT       = 2
    FORWARD     = 3
    DOWN        = 4
    WAIT_UP     = 5
    WAIT_LEFT   = 6
    WAIT_RIGHT  = 7
    WAIT_DOWN   = 8
end

struct GridWaypoint
    x::Float64
    y::Float64
    transition::GridTransition
end

struct GridWorldPath
    waypoints::Vector{GridWaypoint}
end

end # End of module RubberDucks (This comment is nice to distinguish from an end)
