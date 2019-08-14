let
    for i in 1:8
        GridTransition(i)
    end
end

let
    GridWaypoint(0.0,0.0,GridTransition(1))
end

let
    GridWorldPath([
        GridWaypoint(0.0,0.0,GridTransition(1)),
        GridWaypoint(0.0,1.0,GridTransition(1)),
    ])
end
