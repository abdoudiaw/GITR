#include <map>
#include <string>

#if USE_DOUBLE
using gitr_precision = double;
#else
using gitr_precision = float;
#endif

struct MaterialProperties {
    std::string name;
    gitr_precision surfaceBindingEnergy;

    MaterialProperties() : name(""), surfaceBindingEnergy(0.0) {} // This is where the default constructor should be

    MaterialProperties(const std::string& name_, gitr_precision surfaceBindingEnergy_)
        : name(name_), surfaceBindingEnergy(surfaceBindingEnergy_) {}
};

std::map<int, MaterialProperties> materialData = {
    {1,  MaterialProperties("H", 0.1)},   
    {2,  MaterialProperties("He", 0.2)},
    {3,  MaterialProperties("Li", 0.3)},
    {74, MaterialProperties("W", 8.4)}
};
