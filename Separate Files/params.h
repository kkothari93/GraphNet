#ifndef __params__
#define __params__

#define DIM 2
#define TIME_STEP 1e-2
#define SIM_TIME 200
#define TOL 1e-6
#define STEPS int(SIM_TIME/TIME_STEP)
#define L_MEAN 250.0f
#define L_STD 75.0f
#define MAXBOUND 500.0f
#define PAD MAXBOUND*1.03
#define SACBONDS false
#define IMPLEMENT_PBC true
#define FLDR_STRING "set15"
#define CRACKED true
#define RATE_DAMAGE false

#if CRACKED
#define PROB_REMOVAL 1.0
#else
#define PROB_REMOVAL 0.0
#endif

#define __constants__

#define kB 1.38064852e-5
#define b_poly 0.1
#define T 300.0
#define ae 0.1
#define delxe 0.15
#define af 1.0
#define delxf 0.25

#endif
