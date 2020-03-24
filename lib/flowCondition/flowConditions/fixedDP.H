#ifndef fixedDP_H
#define fixedDP_H

#include "flowCondition.H"

namespace Foam
{

class fixedDP
:
    public flowCondition
{

    eulerMesh& emesh_;

    scalar dPdx_;
    
    word flowDir_;
    
public:

    //- Runtime type information
        TypeName("fixedDP");

    //- Constructor
        fixedDP
        (
            eulerMesh& emesh,
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
