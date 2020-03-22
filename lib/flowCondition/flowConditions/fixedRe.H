#ifndef fixedRe_H
#define fixedRe_H

#include "flowCondition.H"
#include "dynamicFvMesh.H"

namespace Foam
{

class fixedRe
:
    public flowCondition
{

    dynamicFvMesh& mesh_;

    scalar Reynolds_;
    scalar dPdx_;
    scalar Dh_; 
    scalar nu_;
    scalar tol_;
    scalar prevReynolds_;

    word flowDir_;
    word outletName_;

    bool addPressureGradientField_;
    
    void initiateDPdx();
    
public:

    //- Runtime type information
        TypeName("fixedRe");

    //- Constructor
        fixedRe
        (
            dynamicFvMesh& mesh,
            const dictionary& dict
        );

    //- Destructor
        virtual ~fixedRe() {} 

    //- Member funcs

        virtual void correctDPdx();

        virtual volVectorField pressureField(); 

}; // End class 

} // End namespace Foam

#endif
