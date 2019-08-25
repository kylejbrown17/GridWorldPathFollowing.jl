let
    # for i in 1:5
    #     GridTransition(i)
    # end
end
let
    GridWaypoint(0.0,0.0,0.0,GridTransition(1))
end
let
    construct_grid_world_path(VecE2(0.0,0.0), 0.0,
        [WEST,WEST,NORTH,NORTH,EAST,WAIT,SOUTH], 1.0, 1.0)
    construct_grid_world_path(
        0.0,0.0,0.0,1.0,1.0,
        [2,2,1,1,0,4,3])
end
