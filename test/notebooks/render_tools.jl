using DataFrames, Gadfly

function plot_traj_with_steps(traj,grid_path;Δt=0.01,pad=0.01)
    t_range = get_start_time(traj):Δt:get_end_time(traj)
    cw = grid_path.cellwidth/2 - pad
    xpts = [map(w->w.pt.x, grid_path.waypoints)..., grid_path.start_pt.x]
    ypts = [map(w->w.pt.y, grid_path.waypoints)..., grid_path.start_pt.y]
    plot(
        layer( x=[get_position(traj,t).x for t in t_range], y=[get_position(traj,t).y for t in t_range],
            Geom.path, Theme(default_color="black") ),
        layer( x=[get_start_pt(seg).x for seg in traj.segments], y=[get_start_pt(seg).y for seg in traj.segments],
            size=[8pt], Geom.point, Theme(default_color="red") ),
        layer( x=[get_end_pt(seg).x for seg in traj.segments], y=[get_end_pt(seg).y for seg in traj.segments],
            size=[8pt],Geom.point, Theme(default_color="yellow") ),
        layer( xmin=xpts .- cw, xmax=xpts .+ cw, ymin=ypts .- cw, ymax=ypts .+ cw, Geom.rect ),
        Coord.cartesian(fixed=true)
    )
end
function compare_velocity_profiles(traj1,traj2;Δt=0.05)
    t0 = max(get_start_time(traj1),get_start_time(traj2))
    tf = min(get_end_time(traj1),get_end_time(traj2))
    @assert t0 < tf "trajectories do not overlap"
    t_vec = collect(get_start_time(base_traj):Δt:get_end_time(base_traj))
    df = DataFrame(
        t=t_vec, 
        v1=[norm(get_vel(traj1,t)) for t in t_vec],
        v2=[norm(get_vel(traj2,t)) for t in t_vec],
    )
    plot(
        layer(df,x="t",y="v2",Geom.path,Theme(default_color="red")),
        layer(df,x="t",y="v1",Geom.path,Theme(default_color="black")),
        Guide.manual_color_key("",["v1","v2"],["black","red"])
    )
end
function plot_vel_and_yaw_rate(traj;Δt=0.1)
    t_vec = collect(get_start_time(traj):Δt:get_end_time(traj))
    df = DataFrame(
        t=t_vec, 
        theta=[atan(get_heading(traj,t)) for t in t_vec],
        s=[get_dist(traj,t) for t in t_vec], 
        v=[norm(get_vel(traj,t)) for t in t_vec] )
    vstack( plot(df,x="t", y="theta",Geom.path), 
        plot(df,x="t",y="v",Geom.path,Theme(default_color="red")) )
end
function plot_accel_vel_pos(base_traj,traj,t_vec,accel,vel,pos)
    df = DataFrame(
        t = t_vec,
        a = [accel..., 0.0],
        v = vel,
        v1=map(t->norm(get_vel(base_traj,t)),t_vec),
        s = pos
    )
    df2 = DataFrame(
        t = [[get_start_time(seg) for seg in traj.segments]..., get_end_time(traj)],
        s = [0, cumsum([get_length(seg) for seg in traj.segments])...]
    )
#     vstack(
        plot(
            layer(df,x="t",y="a",Geom.step,Theme(default_color="green")),
            layer(df,x="t",y="v",Geom.path,Theme(default_color="red")),
            layer(df,x="t",y="v1",Geom.path,Theme(default_color="black")),
            Guide.manual_color_key("",["accel","vel_opt","vel"],["green","red","black"])
        )
#         plot(
#             layer(df,x="t",y="s",Geom.path,Theme(default_color="blue")),
#             layer(df2,x="t",y="s",Geom.point,Theme(default_color="orange")),
#             Guide.manual_color_key("",["distance","constraints"],["blue","orange"])
#         )
#     )
end
function summarize_simulation(traj,states,cmds,time_vec)
    if length(time_vec) < length(states)
        time_vec = [time_vec..., 2*time_vec[end]-time_vec[end-1]]
    end
    errors = Vector{Vector{Float64}}()
    for (t,state) in zip(time_vec,states)
        target_pt = get_trajectory_point_by_time(traj, t)
        target = [target_pt.pos.x, target_pt.pos.y, atan(target_pt.heading)]
        push!(errors, [target[1]-state[1],target[2]-state[2],get_angular_offset(state[3],target[3])])
    end
    df = DataFrame(
        t=time_vec,
        w=[map(u->u[1],cmds)...,0.0],
        v=[map(u->u[2],cmds)...,0.0],
        wr=map(t->get_yaw_rate(traj,t),time_vec),
        vr=map(t->norm(get_vel(traj,t)),time_vec),
        x=map(s->s[1],states[1:length(time_vec)]),
        y=map(s->s[2],states[1:length(time_vec)]),
        θ=map(s->wrap_to_pi(s[3]),states[1:length(time_vec)]),
        xr=map(t->get_position(traj,t).x,time_vec),
        yr=map(t->get_position(traj,t).y,time_vec),
        θr=map(t->atan(get_heading(traj,t)),time_vec),
        ex=map(e->e[1],errors),
        ey=map(e->e[2],errors),
        eθ=map(e->e[3],errors)
    )
    hstack(
            vstack(
                plot(
                    layer(df,x="x",y="y",Geom.path,Theme(default_color="blue")),
                    layer(x=[states[1][1]],y=[states[1][2]],Geom.point,Theme(default_color="blue")),
                    layer(df,x="xr",y="yr",Geom.path,Theme(default_color="red")),
                    layer(x=[get_position(traj,time_vec[1]).x],
                        y=[get_position(traj,time_vec[1]).y],Geom.point,Theme(default_color="red")),
                    Coord.Cartesian(fixed=true),
                    Guide.manual_color_key("",["true","ref"],["blue","red"])
                )
            ),
            vstack(
                plot(
                    layer(df,x="t",y="w",Geom.path,Theme(default_color="red")),
                    layer(df,x="t",y="wr",Geom.path,Theme(default_color="black")),
                    Guide.manual_color_key("",["w","wref"],["red","black"]) ),
                plot(
                    layer(df,x="t",y="v",Geom.path,Theme(default_color="blue")),
                    layer(df,x="t",y="vr",Geom.path,Theme(default_color="black")),
                    Guide.manual_color_key("",["v","vref"],["blue","black"]) ),
                plot(
                    layer(df,x="t",y="ex",Geom.path,Theme(default_color="blue")),
                    layer(df,x="t",y="ey",Geom.path,Theme(default_color="red")),
                    layer(df,x="t",y="eθ",Geom.path,Theme(default_color="green")),
                    Guide.manual_color_key("",["ex","ey","eθ"],["blue","red","green"]) )
            )
        )
end