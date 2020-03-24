#ifndef zeroPressureGradient_h
#define zeroPressureGradient_H

#include "flowCondition.H"

namespace Foam
{

class zeroPressureGradient: public flowCondition
{

    eulerMesh& emesh_;

    bool addPressureGradientField_;
    
    void initiateDPdx();
    
public:

    //- Runtime type information
        TypeName("zeroPressureGradient");

    //- Constructor
        zeroPressureGradient
        (
            eulerMesh& emesh,
            const dictionary& dict
        );

    //- Destructor
        virtual ~zeroPressureGradient() {} 

    //- Member funcs

        virtual void correctDPdx();

        virtual volVectorField pressureField(); 

}; // End class 

} // End namespace Foam

#endif
