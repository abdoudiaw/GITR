#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include "flags.h"

bool Flags::initialize_flags(libconfig::Config &cfg,std::string s) 
{
  // std::string base = "flags.";
  // int flag = 0;
  // bool success = getVariable_cfg<int>(cfg, base+s); // assuming getVariable_cfg can return a bool for success
  // if (!success) {
  //   std::cout << "Error: could not find flag " << s << " in config file." << std::endl;
  //   exit(1);
  // }

  // return flag;
  // int flag = getVariable_cfg<int> (cfg, base+s);
  // if(flag > 0) return true;
  // else return false;
}
