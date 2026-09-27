// Stationary qubit-only accelerator; shared lexicographic prefixes are reused.
// No fitting, sampling or I/O. All gates are the same physical affine maps used
// by the independent complex Kraus implementation.
#include <cstdint>
#include <vector>
extern "C" void qubit_prefix(const int64_t* offsets, const int8_t* words,
    const int64_t* reuse, int64_t count, const double* gates,
    const double* initial, const double* effect, double* output) {
  int64_t longest=0;
  for(int64_t i=0;i<count;++i) {
    int64_t length=offsets[i+1]-offsets[i];
    if(length>longest)longest=length;
  }
  std::vector<double> states(3*(longest+1));
  for(int a=0;a<3;++a)states[a]=initial[a+1];
  for(int64_t i=0;i<count;++i) {
    int64_t length=offsets[i+1]-offsets[i];
    for(int64_t j=reuse[i];j<length;++j) {
      const double* G=gates+25*words[offsets[i]+j];
      const double* r=&states[3*j];
      double* next=&states[3*(j+1)];
      for(int a=0;a<3;++a) {
        const double* row=G+5*(a+1);
        next[a]=row[0]+row[1]*r[0]+row[2]*r[1]+row[3]*r[2];
      }
    }
    const double* r=&states[3*length];
    output[i]=effect[0]+effect[1]*r[0]+effect[2]*r[1]+effect[3]*r[2];
  }
}
