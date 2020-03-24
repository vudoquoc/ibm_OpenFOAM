/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       flowCondition class
       +===   \===\\ =//        OpenFOAM 6.0
\*---------------------------------------------------------------------------*/

#ifndef flowCondition_H
#define flowCondition_H

#include "dictionary.H"
#include "eulerMesh.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class flowCondition
{

public:

    TypeName("flowCondition");

    declareRunTimeSelectionTable
    (
        autoPtr,
        flowCondition,
        dictionary,
        (
            eulerMesh& emesh,
            const dictionary& dict
        ),
        (emesh, dict)
    );

    flowCondition
    (
        eulerMesh& emesh,
        const dictionary& dict
    );

    static autoPtr<flowCondition> New  
    (
        eulerMesh& emesh,
        const dictionary& dict
    );


    virtual ~flowCondition(){}

    virtual void correctDPdx() = 0;

    virtual volVectorField pressureField() = 0;

}; // End class

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
