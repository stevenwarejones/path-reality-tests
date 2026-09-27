// Optional deterministic accelerator. It has no I/O, randomness, or fitting logic.
#include <cstdint>
extern "C" void propagate(const int64_t* offsets, const int8_t* words, int64_t count,
                          const double* gates, const double* initial, const double* effect,
                          int alternate, double* output) {
  // gates[branch][gate][row][column]; order of coordinates: mass,x,y,z,leak.
  for(int64_t i=0; i<count; ++i) {
    double sum=0.;
    for(int branch=0; branch<2; ++branch) {
      double state[5]; for(int a=0;a<5;++a) state[a]=initial[a];
      int phase=branch;
      for(int64_t t=offsets[i]; t<offsets[i+1]; ++t) {
        const double* g=gates+(phase*3+words[t])*25;
        double next[5];
        for(int a=0;a<5;++a) {
          double v=0.; for(int b=0;b<5;++b) v+=g[a*5+b]*state[b];
          next[a]=v;
        }
        for(int a=0;a<5;++a) state[a]=next[a];
        if(alternate) phase=1-phase;
      }
      for(int a=0;a<5;++a) sum+=effect[a]*state[a]/2.;
    }
    output[i]=sum;
  }
}
