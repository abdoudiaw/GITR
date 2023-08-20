#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif


gitr_precision CouletteManfredi(gitr_precision distance, gitr_precision alpha); 

gitr_precision BrooksModel(gitr_precision distance, gitr_precision alpha, gitr_precision larmorRadius, gitr_precision Te);