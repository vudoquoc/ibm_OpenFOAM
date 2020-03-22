#ifndef fixedDP_H
#define fixedDP_H

#include "flowCondition.H"
#include "dynamicFvMesh.H"

namespace Foam
{

class fixedDP
:
    public flowCondition
{

    dynamicFvMesh& mesh_;

    scalar dPdx_;
    
    word flowDir_;
    
public:

    //- Runtime type information
        TypeName("fixedDP");

    //- Constructor
        fixedDP
        (
            dynamicFvMesh& mesh,
            const dictionary& dict
        );

    //- Destructor
        virtual ~fixedDP() {} 

    //- Member funcs

		virtual void correctDPdx();
		
        virtual volVectorField pressureField(); 

}; // End class 

} // End namespace Foam

#endif
