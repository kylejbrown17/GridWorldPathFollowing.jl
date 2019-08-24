#include <julia.h>
#include <stdio.h>

// Compile this progam into an executable called embed_example by running the following command:
// ~/.local/opt/julia-1.0.2/share/julia/julia-config.jl --cflags --ldflags --ldlibs | xargs gcc -o ~/Desktop/embed_example ~/Desktop/embed_example.c
enum transition{EAST,NORTH,WEST,SOUTH,WAIT};

int main(int argc, char *argv[])
{
    /* required: setup the Julia context */
    jl_init();

    /* run Julia commands */
    // jl_eval_string("using Pkg");
    jl_eval_string("using GridWorldPathFollowing");
    jl_eval_string("warmup()");

    jl_eval_string("controller = SwitchingController()");
    jl_eval_string("model  = UnicycleModel()");

    const double[] start_pt = {0.0,0.0};
    const double start_time = 0.0;
    const transition[] transitions = {SOUTH,EAST,WAIT,NORTH,EAST,WAIT,WAIT,EAST,EAST,SOUTH};
    printf("transitions: %i %i %i", transitions[0],transitions[1],transitions[2]);


    jl_eval_string("grid_path=construct_grid_world_path()")

    (void)jl_eval_string("println(sqrt(2.0))");
    // jl_eval_string("using GridWorldPathFollowing");
    jl_eval_string("a = 84");
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
    jl_atexit_hook(0);
    return 0;
}
