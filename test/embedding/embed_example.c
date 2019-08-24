#include <julia.h>
#include <stdio.h>

JULIA_DEFINE_FAST_TLS()
// Compile this progam into an executable called embed_example by running the following command:
// ~/.local/opt/julia-1.0.2/share/julia/julia-config.jl --cflags --ldflags --ldlibs | xargs gcc -o ~/Desktop/embed_example ~/Desktop/embed_example.c
enum TRANSITION{EAST,NORTH,WEST,SOUTH,WAIT};

int main(int argc, char *argv[])
{
  /* required: setup the Julia context */
  jl_init();
  // (void)jl_gc_enable(0);

  // Initialize refs dict to protect objects from being Garbage Collected
  // jl_value_t* refs = jl_eval_string("refs = IdDict()");
  // jl_function_t* setindex = jl_get_function(jl_base_module, "setindex!");
  // jl_datatype_t* reft = (jl_datatype_t*)jl_eval_string("Base.RefValue{Any}");

  /* run Julia commands */
  jl_eval_string("using Vec");
  jl_eval_string("using GridWorldPathFollowing");
  jl_eval_string("GridWorldPathFollowing.warmup()");
  // Get the optimized trajectory
  jl_module_t *GridWorldPathFollowing = (jl_module_t *)jl_eval_string("GridWorldPathFollowing");
  jl_function_t *construct_trajectory = jl_get_function(GridWorldPathFollowing, "construct_trajectory");
  jl_value_t *grid_path = jl_eval_string("grid_path = construct_grid_world_path(VecE2(0.0,0.0),0.0,[SOUTH,EAST,WAIT,NORTH,EAST,WAIT,WAIT,EAST,EAST,SOUTH],0.5,4.0)");
  if (jl_exception_occurred())
    printf("Exception occured: %s \n", jl_typeof_str(jl_exception_occurred()));
  // jl_value_t *base_traj = jl_eval_string("base_traj = construct_trajectory(construct_grid_world_path(VecE2(0.0,0.0),0.0,[SOUTH,EAST,WAIT,NORTH,EAST,WAIT,WAIT,EAST,EAST,SOUTH],0.5,4.0))");
  jl_value_t *base_traj = jl_call1(construct_trajectory, grid_path);
  // optimize_velocity_profile currently passes old objects to the new_traj. This may cause GC problems
  jl_function_t *optimize_velocity_profile_traj_only = jl_get_function(GridWorldPathFollowing, "optimize_velocity_profile_traj_only");
  jl_value_t *traj = jl_call1(optimize_velocity_profile_traj_only, base_traj);

  // (void)jl_eval_string("println(\"HELLO\")");
  // (void)jl_eval_string("println(typeof(refs))");
  // (void)jl_eval_string("println(get_length(base_traj))");
  // (void)jl_eval_string("println(get_end_time(base_traj))");
  // (void)jl_eval_string("println(get_vel(base_traj,0.5))");

  // (void)jl_gc_enable(1);
  // (void)jl_gc_collect();
  // // initialize the controller
  // jl_value_t *controller = jl_eval_string("SwitchingController()");

  // jl_value_t *model = jl_eval_string("UnicycleModel()");

  // const double start_pt[] = {0.0,0.0};
  // jl_value_t *start_time = jl_box_float64(0.0);
  // enum TRANSITION transitions[] = {SOUTH,EAST,WAIT,NORTH,EAST,WAIT,WAIT,EAST,EAST,SOUTH};
  // jl_value_t *cellwidth = jl_box_float64(0.5);

  // jl_value_t *new_line = jl_eval_string("\"\n\"");
  // jl_function_t *print_function = jl_get_function(jl_base_module, "print");
  //
  // jl_eval_string("using LinearAlgebra");
  // jl_module_t *LinearAlgebra = (jl_module_t *)jl_eval_string("LinearAlgebra");
  // jl_function_t *norm_function = jl_get_function(LinearAlgebra, "norm");
  // // jl_function_t *norm_function = jl_get_function(jl_main_module, "LinearAlgebra.norm");
  // jl_value_t *my_arr = jl_eval_string("[1.0, 0.0]");
  // jl_value_t *norm_val = jl_call1(norm_function, my_arr);
  // double norm_val_unboxed = jl_unbox_float64(norm_val);
  // printf("norm_val_unboxed = %f\n",norm_val_unboxed);
  //
  // (void)jl_call1(print_function, start_time);
  // (void)jl_call2(print_function, start_time, new_line);
  // printf("transitions: %i %i %i\n", transitions[0],transitions[1],transitions[2]);

  (void)jl_eval_string("println(sqrt(2.0))");
  // jl_eval_string("using GridWorldPathFollowing");
  (void)jl_eval_string("a = 84");
  (void)jl_eval_string("println(a)");

  jl_value_t *ret = jl_eval_string("sqrt(2.0)");

  if (jl_typeis(ret, jl_float64_type)) {
          double ret_unboxed = jl_unbox_float64(ret);
              printf("sqrt(2.0) in C: %e \n", ret_unboxed);
  }
  else {
          printf("ERROR: unexpected return type from sqrt(::Float64)\n");
  }

  /* strongly recommended: notify Julia that the
   * program is about to terminate. this allows
   * Julia time to cleanup pending write requests
   * and run all finalizers
   * */
  // JL_GC_POP();
  // JL_GC_POP();
  jl_atexit_hook(0);
  return 0;
}
