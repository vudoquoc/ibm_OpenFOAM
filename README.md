# ibm_OpenFOAM

### Immersed Boundary Library in OpenFOAM for fluid-structure interaction simulations

Resolved:
16/3: Add flowCondition class
21/3: Remove IBObjectRegistry class.

Remaining Issues:
- Access mesh_ from emesh_ so as to avoid passing mesh_ variable in all major classes' constructors
- Move updateObjectMotionUhlmann/Tobias from IBParticle to corresponding IBMotionSolver class
