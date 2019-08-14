let
    for i in 1:5
        GridTransition(i)
    end
end
let
    GridWaypoint(0.0,0.0,0.0,GridTransition(1))
end
let
    construct_grid_world_path(VecE2(0.0,0.0), 0.0, 
        [LEFT,LEFT,UP,UP,RIGHT,WAIT,DOWN], 1.0, 1.0)
end
