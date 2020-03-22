/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - newton class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#include "newton.H"
#include "wallPolyPatch.H"
#include "IBParticle.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(newton, 0);
    addToRunTimeSelectionTable(IBMotionSolver, newton, dictionary);
}

//---------------------------------Constructors------------------------------//
Foam::newton::newton
(
    word typeName,
	dynamicFvMesh& mesh,
    eulerMesh& emesh,
    PtrList<IBObject>& ibo,
    autoPtr<IBModel>& ibmodelPtr
)
:
    IBMotionSolver(typeName, mesh, emesh, ibo, ibmodelPtr),
	mesh_(mesh),
    emesh_(emesh),
    ibo_(ibo),
    ibmodelPtr_(ibmodelPtr)
{}

// -------------------------------Member Functions----------------------------//

void  Foam::newton::moveObjects()
{
    for(int i=0; i<ibo_.size(); i++)
    {
        if (ibo_[i].movable())
        {
            Info<< "Moving object "<<ibo_[i].name()<<endl;
            
            vector repulsiveForce(vector::zero);
            
            for(int j=0; j<ibo_.size(); j++)
            {
                if (i == j)
                    continue;
                else
                {
                    repulsiveForce += objMutualRepulsive(i, j);
                }
            }
            
            ibo_[i].updateObjectMotionUhlmann
            (
                ibmodelPtr_->FLagrange()[i],
                repulsiveForce,
                emesh_.rhoFluid(),
                emesh_.g()
            );

            ibo_[i].updateVelocity();
            ibo_[i].movePoints();
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //