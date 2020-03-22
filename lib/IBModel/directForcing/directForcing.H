/*---------------------------------------------------------------------------*\
       +===   /===\  ==   ===   Hanoi University of Science and Technology
       ||  \\     \\ ||  //     School of Transportation Engineering
       ||   ))     )) | //      Computational Engineering Group
       ||  //    \\/  |//       IBMlibs - directForcing class
       +===   \===\\ =//        OpenFOAM 5.0 - 13/4/2018
\*---------------------------------------------------------------------------*/

#ifndef directForcing_H
#define directForcing_H

#include "IBModel.H"
#include "IBObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
							Class IBModel declaration
\*---------------------------------------------------------------------------*/

class directForcing
: 
	public IBModel
{
	//- Private data
		//- Reference to mesh
		dynamicFvMesh& mesh_;
		
		//- Reference to euler mesh
		eulerMesh& emesh_;

		//- Reference to IB objects
		PtrList<IBObject>& ibo_;

		//- Number of directForcing loop
		label nMDF_;
		
		List<vectorField> uLagrang_;

		List<vectorField> fLagrang_;
		
	//- Disallow default bitwise copy construct
	directForcing(const directForcing&);
	
	//- Disallow default bitwis assignment
	void operator=(const directForcing&);

	
public:

	//- Runtime type information
	TypeName("directForcing");

	//- Constructor
		directForcing
		(
			dynamicFvMesh& mesh,
			eulerMesh& emesh,
			const dictionary& dict,
			PtrList<IBObject>& ibo
		);

	//- Destructor
		virtual ~directForcing()
		{}

	//- Member functions

		virtual const List<vectorField>& ULagrange(){return uLagrang_;}
		virtual const List<vectorField>& FLagrange(){return fLagrang_;}
		virtual volVectorField ibForce(const volVectorField& U);

		// virtual volVectorField ibForceInt();
		
		// //- Perform multi-Firect forcing
		// virtual void multiDirectForcing
		// (
		// 	volVectorField& u,
		// 	volVectorField& ibForce
		// );

}; //End class

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
// ************************************************************************* //
