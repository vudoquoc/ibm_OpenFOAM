/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - IBModel class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#ifndef flowCondition_H
#define flowCondition_H

#include "dictionary.H"
#include "dynamicFvMesh.H"

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
            dynamicFvMesh& mesh,
            const dictionary& dict
        ),
        (mesh, dict)
    );

    flowCondition
    (
        dynamicFvMesh& mesh,
        const dictionary& dict
    );

    static autoPtr<flowCondition> New  
    (
        dynamicFvMesh& mesh,
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
