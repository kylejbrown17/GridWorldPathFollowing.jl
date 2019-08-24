#include <julia.h>
#include <stdio.h>

JULIA_DEFINE_FAST_TLS()
// Compile this progam into an executable called embed_example by running the following command:
// ~/.local/opt/julia-1.0.2/share/julia/julia-config.jl --cflags --ldflags --ldlibs | xargs gcc -o ~/Desktop/embed_example ~/Desktop/embed_example.c
// enum TRANSITION{EAST,NORTH,WEST,SOUTH,WAIT};

int main(int argc, char *argv[])
{
  /* required: setup the Julia context */
  jl_init();

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

  // initialize the controller and robot model
  jl_value_t *controller = jl_eval_string("SwitchingController()");
  jl_value_t *model = jl_eval_string("UnicycleModel()");
  if (jl_exception_occurred())
    printf("Exception occured: %s \n", jl_typeof_str(jl_exception_occurred()));

  double state[] = {0.0,0.1,0.0};
  jl_value_t* array_type = jl_apply_array_type((jl_value_t*)jl_float64_type, 1);
  jl_array_t *state_vec = jl_ptr_to_array_1d(array_type, state, 10, 0);

  jl_function_t *get_action = jl_get_function(GridWorldPathFollowing, "get_action");
  jl_value_t *t = jl_box_float64(3.0); // time
  jl_value_t *f_args[4] = {NULL};
  f_args[0] = controller;
  f_args[1] = traj;
  f_args[2] = (jl_value_t*)state_vec;
  f_args[3] = t;
  jl_array_t *cmd = (jl_array_t*)jl_call(get_action,f_args,4);
  double *cmdData = (double*)jl_array_data(cmd);

  printf("cmd: w=%f,v=%f\n", cmdData[0], cmdData[1]);

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
