/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       immersedBoundaryMethod class
       +===   \===\\ =//        OpenFOAM 6.0
\*---------------------------------------------------------------------------*/

#include "immersedBoundaryMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(immersedBoundaryMethod, 0);
}
//----------------------------Private Member Functions-----------------------//
void Foam::immersedBoundaryMethod::createObjects(const dictionary& dict)
{
    label i=0;
    scalar nObjects = dict.size();
    objects_.setSize(nObjects);
    Info<<"Found "<<nObjects<<" IB objects!!!"<<endl;
    forAllConstIter(IDLList<entry>, dict, iter)
    {
        if (iter().isDict())
        {
            objects_.set
            (
                i++,
                IBObject::New 
                (
                    iter().keyword(),
                    emesh_,
                    iter().dict()
                )
            );
        }
    }
    objects_.setSize(i);
}

//---------------------------------Constructors------------------------------//
Foam::immersedBoundaryMethod::immersedBoundaryMethod
(
	dynamicFvMesh& mesh
)
:	
	IOdictionary
	(
		IOobject
		(
			"IBMDict",
			mesh.time().constant(),
			mesh,
        	IOobject::MUST_READ,
        	IOobject::NO_WRITE
		)
	),
	emesh_(mesh)
{
	// Create flow condition for microchannels
	flowConditionPtr_ = flowCondition::New(emesh_,subDict("flowCondition"));

	// Create immersed objects
	createObjects(subDict("IBObjects"));

	// Create immersed boundary model
	modelPtr_  = IBModel::New
				 (
					emesh_, 
					subDict("IBModel"), 
					objects_
				 );
	
	// Create motion solver for ibobjects motions
	motionSolverPtr_ = IBMotionSolver::New
				 (
					subDict("IBModel").lookup("motionSolver"), 
					emesh_, 
					objects_, 
					modelPtr_
				 );
}

// -------------------------------Member Functions----------------------------//

Foam::volVectorField Foam::immersedBoundaryMethod::ibForce
(
	const volVectorField& U
)
{
   return modelPtr_->ibForce(U);
}

void Foam::immersedBoundaryMethod::update()
{
	return motionSolverPtr_->moveObjects();
}
/*
void Foam::immersedBoundaryMethod::write()
{
	return modelPtr_->write();
}
*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
