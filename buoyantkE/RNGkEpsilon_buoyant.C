/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "RNGkEpsilon_buoyant.H"
#include "fvOptions.H"
#include "bound.H"
#include "gravityMeshObject.H"
#include "fvcGrad.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void RNGkEpsilon_buoyant<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = Cmu_*sqr(k_)/epsilon_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> RNGkEpsilon_buoyant<BasicTurbulenceModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()
            /dimTime
        )
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> RNGkEpsilon_buoyant<BasicTurbulenceModel>::epsilonSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            epsilon_,
            dimVolume*this->rho_.dimensions()*epsilon_.dimensions()
            /dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
RNGkEpsilon_buoyant<BasicTurbulenceModel>::RNGkEpsilon_buoyant
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    Cmu_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.0845
        )
    ),
    C1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.42
        )
    ),
    C2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C2",
            this->coeffDict_,
            1.68
        )
    ),
    C3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C3",
            this->coeffDict_,
            0
        )
    ),

    Cg_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cg",
            this->coeffDict_,
            1.0
        )
    ),

    sigmak_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmak",
            this->coeffDict_,
            0.71942
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            0.71942
        )
    ),
    eta0_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "eta0",
            this->coeffDict_,
            4.38
        )
    ),
    beta_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "beta",
            this->coeffDict_,
            0.012
        )
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool RNGkEpsilon_buoyant<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
	Cg_.readIfPresent(this->coeffDict());
        sigmak_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());
        eta0_.readIfPresent(this->coeffDict());
        beta_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void RNGkEpsilon_buoyant<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    const volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));
    
    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    const volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))().v()
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    const volScalarField::Internal GbyNu
    (
        IOobject::scopedName(this->type(), "GbyNu"),
        tgradU().v() && devTwoSymm(tgradU().v())
    );
    tgradU.clear();

    const volScalarField::Internal G(this->GName(), nut()*GbyNu);

    const volScalarField::Internal eta(sqrt(mag(GbyNu))*k_/epsilon_);
    const volScalarField::Internal eta3(eta*sqr(eta));

    const volScalarField::Internal R
    (
        ((eta*(-eta/eta0_ + scalar(1)))/(beta_*eta3 + scalar(1)))
    );

    // --- Declare gravity vector field here ---
    const uniformDimensionedVectorField& g =
        meshObjects::gravity::New(this->mesh_.time());
    const volScalarField& rhoVol =
    this->mesh_.objectRegistry::template lookupObject<volScalarField>("rho");
//    const volScalarField Gcoeff
//    (
//    (Cg_ * this->Cmu_) * this->alpha_ * this->k_ * (g & fvc::grad(rhoVol))
//   / (this->epsilon_ + this->epsilonMin_)
//    );

const dimensionedScalar epsilonFloor(
    "epsilonFloor",
    this->epsilon_.dimensions(),  // same dimensions as epsilon_
    1e-6                          // numerical value
);
// Compute Gcoeff safely
    // Compute Gcoeff safely using this-> members
const volScalarField Gcoeff
(
    ((this->Cg_ * this->Cmu_) * this->alpha_ * this->k_ * (g & fvc::grad(rhoVol))) /
    max(this->epsilon_ + this->epsilonMin_, epsilonFloor) // prevents huge values
);
    const dimensionedScalar dummy("dummy", Gcoeff.dimensions(), 1.0);
    const dimensionedScalar dummyT("dummy", dimensionSet(0, 0, 1, 0, 0, 0, 0), 1.0);
    const volScalarField GcoeffDimless = Gcoeff / dummy / dummyT;
// Compute statistics of GcoeffDimless
scalar G_min = gMin(GcoeffDimless);      // minimum value
scalar G_max = gMax(GcoeffDimless);      // maximum value
scalar G_avg = gSum(GcoeffDimless) / GcoeffDimless.internalField().size();  // average

Info << "GcoeffDimless stats: min = " << G_min
     << ", max = " << G_max
     << ", avg = " << G_avg << endl;

// Compute the dot product g · ∇ρ
const volScalarField gDotGradRho = g & fvc::grad(rhoVol);

// Compute min, max, average
scalar minVal = gMin(gDotGradRho);
scalar maxVal = gMax(gDotGradRho);
scalar avgVal = gSum(gDotGradRho) / gDotGradRho.internalField().size();

// Print to terminal/log
Info << "g · grad(rho) stats: "
     << "min = " << minVal << ", "
     << "max = " << maxVal << ", "
     << "avg = " << avgVal << endl;

    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();
    // Push any changed cell values to coupled neighbours
    epsilon_.boundaryFieldRef().template evaluateCoupled<coupledFvPatch>();


    vector gHat(g.value() / mag(g.value()));

    const volScalarField v = gHat & this->U_;  // tmp<volScalarField>
    const volScalarField u = mag(this->U_ - gHat*v) + dimensionedScalar("SMALL", dimVelocity, SMALL);


    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        (C1_ - R)*alpha()*rho()*GbyNu*Cmu_*k_()
      - fvm::SuSp(((2.0/3.0)*C1_ - C3_)*alpha()*rho()*divU, epsilon_)
      - fvm::Sp(C2_*alpha()*rho()*epsilon_()/k_(), epsilon_)
      -fvm::SuSp(this->C1_*0.5 * GcoeffDimless, this->epsilon_)
      + epsilonSource()
      + fvOptions(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvOptions.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvOptions.correct(epsilon_);
    bound(epsilon_, this->epsilonMin_);


    // Turbulent kinetic energy equation

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha()*rho()*G
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*epsilon_()/k_(), k_)
      -fvm::SuSp(GcoeffDimless, this->k_)
      + kSource()
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
