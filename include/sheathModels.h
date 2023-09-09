#include <vector>
#include <cmath>  
#include "Boundary.h"

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif


enum class SheathModel {
    STANGEBY,
    COULETTE_MANFREDI
};

class ISheathModel {
public:
    virtual gitr_precision compute(gitr_precision minDistance, Boundary* boundaryVector, int minIndex) = 0;
};

class StangebyModel : public ISheathModel {
public:
    gitr_precision compute(gitr_precision minDistance, Boundary* boundaryVector, int minIndex) override {
        gitr_precision angle = boundaryVector[minIndex].angle;   
        gitr_precision fd  = boundaryVector[minIndex].fd;
        gitr_precision potential = boundaryVector[minIndex].potential;
        gitr_precision debyeLength = boundaryVector[minIndex].debyeLength;
        gitr_precision larmorRadius = boundaryVector[minIndex].larmorRadius;

        gitr_precision E_component_DS = fd / (2.0 * debyeLength) * std::exp(-minDistance / (2.0 * debyeLength));
        gitr_precision E_component_CD = (1.0 - fd) / larmorRadius * std::exp(-minDistance / larmorRadius);
        gitr_precision Emag = potential * (E_component_DS + E_component_CD);

        return Emag;
}
};

class CouletteManfrediModel : public ISheathModel {
public:
    gitr_precision compute(gitr_precision minDistance, Boundary* boundaryVector, int minIndex) override {
    gitr_precision alpha = boundaryVector[minIndex].angle; 
    
    gitr_precision C1_alpha = -0.00281407f * alpha - 2.31655435f;
    gitr_precision C2_alpha = 0.00640402f * alpha + 0.01023915f;

    gitr_precision te = boundaryVector[minIndex].te;
    gitr_precision debyeLength = boundaryVector[minIndex].debyeLength;

    gitr_precision normalizationFactor = gitr_constants::k * te * 11605.0/ gitr_constants::e;    // k * te /e 

    gitr_precision thresholdForSheathModel =  200; // 200 debye lengths

    if ( minDistance/debyeLength > thresholdForSheathModel ){
        return 0;
    }
    else{
        gitr_precision result = C1_alpha * C2_alpha *  std::exp(-C2_alpha * abs(minDistance  / debyeLength )) * normalizationFactor / debyeLength ;
        return result;
    }
}
};

gitr_precision sheathModel(SheathModel model, gitr_precision minDistance, Boundary* boundaryVector, int minIndex) {
    ISheathModel* sheath;
    
    switch (model) {
        case SheathModel::STANGEBY:
            sheath = new StangebyModel();
            break;
        case SheathModel::COULETTE_MANFREDI:
            sheath = new CouletteManfrediModel();
            break;
    }

    gitr_precision result = sheath->compute(minDistance, boundaryVector, minIndex);
    delete sheath;  
    return result;
}

//END SHEATH MODEL