// Solve the Momentum equation

fvVectorMatrix UEqn
(
    fvm::ddt(U) 
  + fvm::div(phi, U)
  - fvm::laplacian(nu, U)
);

UEqn.relax();

//volVectorField gradP = IBM.pressureGradientField();

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p) );//+ gradP);
}

volVectorField ibForce = IBM.ibForce(U);

solve(UEqn == -fvc::grad(p) + ibForce);

//IBM.multiDirectForcing(U, ibForce);
