/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - newton class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
newton class:
    - Solve newtons equations for 6DoF motion.
    - Cannot handle objects that have density closed to fluids density
    - M.Uhlmann, J. Comput. Phys. 209, 448–476 (2005).
\*---------------------------------------------------------------------------*/

#ifndef newton_H
#define newton_H

// #include "IOmanip.H"
#include "autoPtr.H"
#include "dynamicFvMesh.H"
#include "eulerMesh.H"
#include "IBObject.H"
#include "IBModel.H"
#include "IBMotionSolver.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
							Class newton declaration
\*---------------------------------------------------------------------------*/
class IBModel;

class newton
:
	public IBMotionSolver
{
	dynamicFvMesh& mesh_;

	eulerMesh& emesh_;

	PtrList<IBObject>& ibo_;
	
	autoPtr<IBModel>& ibmodelPtr_;

public:
	
	TypeName("newton");

	//- Constructor
		newton
		(
			word typeName,
			dynamicFvMesh& mesh, 
			eulerMesh& emesh,
			PtrList<IBObject>& ibo,
			autoPtr<IBModel>& ibmodelPtr
		);
	//- Destructor
		virtual ~newton()
		{}
	
	//- Member functions
		//- Move objects
		virtual void moveObjects();

}; //End class

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //