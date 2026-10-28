#ifndef Farassat1AFormulation_H
#define Farassat1AFormulation_H

#include "fwhFormulation.H"
#include "vector.H"
#include "label.H"
#include "List.H"
#include "scalarField.H"
#include "vectorField.H"
#include "AcousticAnalogy.H"

namespace Foam
{

namespace functionObjects
{

class Farassat1AFormulation : public fwhFormulation
{

private:

    template<typename T>
    using tdPrimitive = AcousticAnalogy::tdPrimitive<T>;

    //- Observer-signal spacing. 0 until the first probe adopts deltaT.
    scalar acousticDeltaT_;

    //- Time origin of signal bin 0.
    scalar tOrigin_;

    //- Surface time derivatives. One copy: L, M and Un do not depend on
    //  the observer. dot() advances the stored state, so it runs once per step.
    List<tdPrimitive<scalar>> Un_;
    List<tdPrimitive<vector>> L_;
    List<tdPrimitive<vector>> M_;

    List<scalarField> dotUn_;
    List<vectorField> dotL_;
    List<vectorField> dotM_;

    //- Source time at which the derivatives above were computed.
    List<scalar> derivTime_;

    //- Sum of thickness and loading integrands on a uniform observer grid.
    //  Indexed [observer][surface].
    List<List<scalarField>> signal_;

    //- Previous emission knot, one per local face.
    List<List<scalarField>> tauPrev_;
    List<List<scalarField>> integrandPrev_;

    //- 1 once a previous emission knot exists. label, not bool: UList<bool>
    //  does not give an assignable reference.
    List<List<labelList>> hasPrev_;

    void ensureAcousticStep();

    void ensureSurfaceStorage(label iSurf, label nFaces);

    void deposit
    (
        scalarField& signal,
        scalar tau0,
        scalar v0,
        scalar tau1,
        scalar v1
    );

    scalar sampleSignal(const scalarField& signal, scalar t) const;

    //- Disable default constructor
    Farassat1AFormulation();

    //- Disable copy constructor
    Farassat1AFormulation(const Farassat1AFormulation&);

public:

    //-
    Farassat1AFormulation(const FfowcsWilliamsHawkings& fwh);

    //-
    ~Farassat1AFormulation();

    //-
    scalar observerAcousticPressure (const vectorField& Sf, const vectorField& uS, const scalarField& rhoS, const scalarField& pS, label iObs, label iSurf, scalar ct);

    //-
    inline void calculateAcousticPressure (const vectorField& Sf, const vectorField& uS, const scalarField& rhoS, const scalarField& pS, label iObs, label iSurf, scalar ct);

    //- Face histories are not stored, so there is nothing to drop.
    void clearExpiredData();

    //-
    void update();

    //-
    void initialize();
};

}

} // namespace Foam

#endif
