#ifndef zeroPressureGradient_h
#define zeroPressureGradient_H

#include "flowCondition.H"
#include "dynamicFvMesh.H"

namespace Foam
{

class zeroPressureGradient: public flowCondition
{

    dynamicFvMesh& mesh_;

    bool addPressureGradientField_;
    
    void initiateDPdx();
    
public:

    //- Runtime type information
        TypeName("zeroPressureGradient");

    //- Constructor
        zeroPressureGradient
        (
            dynamicFvMesh& mesh,
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
