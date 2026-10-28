/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "Farassat1AFormulation.H"
#include "FfowcsWilliamsHawkings.H"
#include "Time.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "mathematicalConstants.H"
#include "Pstream.H"

#include <cmath>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::Farassat1AFormulation::Farassat1AFormulation
(
    const FfowcsWilliamsHawkings& fwh
)
:
    fwhFormulation(fwh, false),
    acousticDeltaT_(0),
    tOrigin_(0)
{
    this->initialize();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::Farassat1AFormulation::~Farassat1AFormulation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::Farassat1AFormulation::initialize()
{
    const label nObs = fwh_.observers_.size();
    const label nSurf = fwh_.controlSurfaces_.size();
    const bool adjust = fwh_.obr_.time().isAdjustTimeStep();

    Un_.resize(nSurf);
    L_.resize(nSurf);
    M_.resize(nSurf);
    dotUn_.resize(nSurf);
    dotL_.resize(nSurf);
    dotM_.resize(nSurf);
    derivTime_.resize(nSurf, -GREAT);

    signal_.resize(nObs);
    tauPrev_.resize(nObs);
    integrandPrev_.resize(nObs);
    hasPrev_.resize(nObs);

    forAll(signal_, iObs)
    {
        signal_[iObs].resize(nSurf);
        tauPrev_[iObs].resize(nSurf);
        integrandPrev_[iObs].resize(nSurf);
        hasPrev_[iObs].resize(nSurf);
    }

    if (fwh_.acousticDeltaT_ > VSMALL)
    {
        acousticDeltaT_ = fwh_.acousticDeltaT_;
    }

    forAll(Un_, iSurf)
    {
        ensureSurfaceStorage(iSurf, fwh_.controlSurfaces_[iSurf].Cf().size());
        Un_[iSurf].setTimeStepVariable(adjust);
        L_[iSurf].setTimeStepVariable(adjust);
        M_[iSurf].setTimeStepVariable(adjust);
    }
}

void Foam::functionObjects::Farassat1AFormulation::ensureAcousticStep()
{
    if (acousticDeltaT_ > VSMALL)
    {
        return;
    }

    acousticDeltaT_ = fwh_.obr_.time().deltaT().value();

    if (acousticDeltaT_ < VSMALL)
    {
        acousticDeltaT_ = 1e-6;
    }

    if (Pstream::master())
    {
        Info<< "Farassat1A acousticDeltaT = " << acousticDeltaT_ << nl;
    }
}

void Foam::functionObjects::Farassat1AFormulation::ensureSurfaceStorage
(
    const label iSurf,
    const label nFaces
)
{
    if
    (
        iSurf < dotUn_.size()
     && dotUn_[iSurf].size() == nFaces
     && Un_[iSurf].size() == nFaces
    )
    {
        return;
    }

    const bool adjust = fwh_.obr_.time().isAdjustTimeStep();

    Un_[iSurf].setTimeStepVariable(adjust);
    Un_[iSurf].resize(nFaces);
    L_[iSurf].setTimeStepVariable(adjust);
    L_[iSurf].resize(nFaces);
    M_[iSurf].setTimeStepVariable(adjust);
    M_[iSurf].resize(nFaces);

    dotUn_[iSurf].resize(nFaces);
    dotL_[iSurf].resize(nFaces);
    dotM_[iSurf].resize(nFaces);

    forAll(hasPrev_, iObs)
    {
        tauPrev_[iObs][iSurf].resize(nFaces);
        integrandPrev_[iObs][iSurf].resize(nFaces);
        hasPrev_[iObs][iSurf].resize(nFaces);
        hasPrev_[iObs][iSurf] = 0;
        signal_[iObs][iSurf].clear();
    }
}

void Foam::functionObjects::Farassat1AFormulation::deposit
(
    scalarField& signal,
    const scalar tau0,
    const scalar v0,
    const scalar tau1,
    const scalar v1
)
{
    if (tau1 <= tau0 + VSMALL || acousticDeltaT_ < VSMALL)
    {
        return;
    }

    const scalar x0 = (tau0 - tOrigin_)/acousticDeltaT_;
    const scalar x1 = (tau1 - tOrigin_)/acousticDeltaT_;
    const label b0 = label(std::floor(x0));
    const label b1 = label(std::floor(x1));

    if (b1 < 0)
    {
        return;
    }

    const label start = max(b0 + 1, 0);

    if (start > b1)
    {
        return;
    }

    // A too-fine grid would bring the per-face history cost back as bins.
    if (b1 > 20000000)
    {
        FatalErrorInFunction
            << "Farassat 1A observer signal needs bin " << b1
            << " at acousticDeltaT " << acousticDeltaT_
            << ". Set a larger acousticDeltaT."
            << exit(FatalError);
    }

    if (signal.size() <= b1)
    {
        signal.resize(b1 + 1, 0.0);
    }

    for (label b = start; b <= b1; ++b)
    {
        const scalar tb = tOrigin_ + scalar(b)*acousticDeltaT_;
        const scalar alpha = (tb - tau0)/(tau1 - tau0);
        signal[b] += v0 + (v1 - v0)*alpha;
    }
}

Foam::scalar Foam::functionObjects::Farassat1AFormulation::sampleSignal
(
    const scalarField& signal,
    const scalar t
) const
{
    if (acousticDeltaT_ < VSMALL || t < tOrigin_ || signal.empty())
    {
        return 0.0;
    }

    scalar x = (t - tOrigin_)/acousticDeltaT_;
    if (x < 0)
    {
        return 0.0;
    }

    const scalar nearest = std::floor(x + 0.5);
    if (mag(x - nearest) < 1e-6)
    {
        x = nearest;
    }

    const label i0 = label(std::floor(x));
    if (i0 < 0 || i0 >= signal.size())
    {
        return 0.0;
    }

    const label i1 = i0 + 1;
    if (i1 >= signal.size())
    {
        return signal[i0];
    }

    const scalar alpha = x - scalar(i0);
    return (1.0 - alpha)*signal[i0] + alpha*signal[i1];
}

Foam::scalar Foam::functionObjects::Farassat1AFormulation::observerAcousticPressure
(
    const vectorField& Sf,
    const vectorField& uS,
    const scalarField& rhoS,
    const scalarField& pS,
    label iObs,
    label iSurf,
    scalar ct
)
{
    calculateAcousticPressure(Sf, uS, rhoS, pS, iObs, iSurf, ct);

    // Same query point as the previous per-face interpolation.
    const scalar ct1 = ct + fwh_.obr_.time().deltaT().value()*1.0e-6;

    scalar local = 0.0;
    if (iObs < signal_.size() && iSurf < signal_[iObs].size())
    {
        local = sampleSignal(signal_[iObs][iSurf], ct1);
    }

    reduce(local, sumOp<scalar>());

    return local/(4.0*constant::mathematical::pi);
}

void Foam::functionObjects::Farassat1AFormulation::calculateAcousticPressure
(
    const vectorField& Sf,
    const vectorField& uS,
    const scalarField& rhoS,
    const scalarField& pS,
    label iObs,
    label iSurf,
    scalar ct
)
{
    ensureAcousticStep();
    ensureSurfaceStorage(iSurf, Sf.size());

    // The first observer to reach this surface at time ct advances the
    // derivatives. Later observers reuse them.
    const bool updateDeriv = (derivTime_[iSurf] != ct);

    forAll(Sf, iFace)
    {
        const vector r = robs_[iObs][iSurf][iFace];
        const scalar magr = magrobs_[iObs][iSurf][iFace];
        const vector rh = r/magr;
        const scalar dS = mag(Sf[iFace]);
        const vector n = ni_[iSurf].value(iFace);

        vector L;
        vector M;
        scalar Un;
        vector dL;
        vector dM;
        scalar dUn;

        if (updateDeriv)
        {
            M = fwh_.vS_[iSurf][iFace]/fwh_.c0_;
            const vector U =
                (1.0 - rhoS[iFace]/fwh_.rhoRef_)*fwh_.vS_[iSurf][iFace]
              + rhoS[iFace]*uS[iFace]/fwh_.rhoRef_;
            const tensor Pf =
                pS[iFace]*tensor::I
              + rhoS[iFace]*uS[iFace]*(uS[iFace] - fwh_.vS_[iSurf][iFace]);

            L = Pf & n;
            Un = U & n;

            Un_[iSurf].value(iFace) = Un;
            L_[iSurf].value(iFace) = L;
            M_[iSurf].value(iFace) = M;

            dUn = Un_[iSurf].dot(ct, iFace);
            dL = L_[iSurf].dot(ct, iFace);
            dM = M_[iSurf].dot(ct, iFace);

            dotUn_[iSurf][iFace] = dUn;
            dotL_[iSurf][iFace] = dL;
            dotM_[iSurf][iFace] = dM;
        }
        else
        {
            Un = Un_[iSurf].value(iFace);
            L = L_[iSurf].value(iFace);
            M = M_[iSurf].value(iFace);
            dUn = dotUn_[iSurf][iFace];
            dL = dotL_[iSurf][iFace];
            dM = dotM_[iSurf][iFace];
        }

        const scalar Mr = M & rh;
        const scalar magM = mag(M);
        const scalar lM = L & M;
        const scalar lr = L & rh;
        const scalar dotlr = dL & rh;
        const scalar dotMr = dM & rh;
        const scalar OneByOneMr = 1.0/(1.0 - Mr);
        const scalar OneByOneMrSq = OneByOneMr*OneByOneMr;

        const scalar q =
        (
            fwh_.rhoRef_*dUn*OneByOneMrSq/magr
          + fwh_.rhoRef_*Un*(magr*dotMr + fwh_.c0_*(Mr - magM*magM))
            *OneByOneMrSq*OneByOneMr/magr/magr
        )*dS;

        const scalar fpart1 = dotlr*(dS/magr/fwh_.c0_)*OneByOneMrSq;
        const scalar fpart2 = (lr - lM)*(dS/magr/magr)*OneByOneMrSq;
        const scalar fpart3 =
            lr*(dS/magr/magr/fwh_.c0_)*OneByOneMrSq*OneByOneMr
           *(magr*dotMr + fwh_.c0_*Mr - fwh_.c0_*magM*magM);

        const scalar integrand = q + fpart1 + fpart2 + fpart3;
        const scalar tau = tobs_[iObs][iSurf][iFace];

        if (hasPrev_[iObs][iSurf][iFace])
        {
            deposit
            (
                signal_[iObs][iSurf],
                tauPrev_[iObs][iSurf][iFace],
                integrandPrev_[iObs][iSurf][iFace],
                tau,
                integrand
            );
        }

        tauPrev_[iObs][iSurf][iFace] = tau;
        integrandPrev_[iObs][iSurf][iFace] = integrand;
        hasPrev_[iObs][iSurf][iFace] = 1;
    }

    if (updateDeriv)
    {
        derivTime_[iSurf] = ct;
    }
}

void Foam::functionObjects::Farassat1AFormulation::update()
{
    fwhFormulation::update();
}

void Foam::functionObjects::Farassat1AFormulation::clearExpiredData()
{}

// ************************************************************************* //
