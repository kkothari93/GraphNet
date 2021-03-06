#ifndef __params__
#define __params__

#define DIM 2
#define TIME_STEP 1e-3
#define SIM_TIME 40
#define TOL 1e-6
#define STEPS int(SIM_TIME/TIME_STEP)
#define L_MEAN 120.0f
#define L_STD 12.0f
#define MAXBOUND_X 500.0f
#define MAXBOUND_Y 100.0
#define PAD_X MAXBOUND_X*1.03
#define PAD_Y MAXBOUND_Y*1.03
#define SACBONDS false
#define IMPLEMENT_PBC false
#define FLDR_STRING "set96_mednetwork_1"
#define CRACKED false
#define RATE_DAMAGE true

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
#define delxe 0.005
#define af 0.3
#define delxf 0.25

#endif
