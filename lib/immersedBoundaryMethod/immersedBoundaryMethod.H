/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - immersedBoundaryMethod class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryMethod_H
#define immersedBoundaryMethod_H

#include "dynamicFvMesh.H"
#include "autoPtr.H"
#include "eulerMesh.H"
#include "flowCondition.H"
#include "IBObject.H"
#include "IBModel.H"
#include "IBMotionSolver.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
							Class immersedBoundaryMethod declaration
\*---------------------------------------------------------------------------*/

class immersedBoundaryMethod
:
	IOdictionary
{
	// Private data

		//- Euler mesh including cell size, cell volume, rho fluid, ...
		eulerMesh emesh_;

		//- Pointer to selected flow condition
		autoPtr<flowCondition> flowConditionPtr_;

		//- IB Objects
		PtrList<IBObject> objects_;

		//- Pointer to selected IB model 
		autoPtr<IBModel> modelPtr_;

		//- Pointer to selected motion solver		
		autoPtr<IBMotionSolver> motionSolverPtr_;

	// Private member functions

		// Create ib objects
		void createObjects(const dictionary& dict);
	
public:
	
	TypeName("immersedBoundaryMethod");
	
	immersedBoundaryMethod(dynamicFvMesh& mesh);

	virtual ~immersedBoundaryMethod()
	{}

	virtual volVectorField ibForce(const volVectorField& U);

	virtual void update();
	
	//virtual void write();

}; //End class

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
// ************************************************************************* //
