#include "forces.H"
#include "volFieldsFwd.H"
#include "Tuple2.H"
#include "OFstream.H"
#include "SoundObserver.H"

#ifndef AcousticAnalogy_H
#define AcousticAnalogy_H

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;


namespace functionObjects
{

class AcousticAnalogy
:
    public forces
{

protected:

        //-
        autoPtr<OFstream> analogyOutPtr_;

        //-
        label probeFreq_;

        //-
        scalar timeStart_;

        //-
        scalar timeEnd_;

        //-
        bool writeFft_;

        //- Sonic speed, m/s
        scalar c0_;

        //- List of observers
        List<SoundObserver> observers_;

        //-
        label probeI_;

protected:

        //-
        void makeFile();

        //-
        void writeFft();

        //-
        virtual void correct();
        
        //-
        virtual void updateSurfaceCf();

public:

    template<class T>
    class tdPrimitive
    {

    private:

        //-
        bool varTimeStep_;

        //- True after one complete pass over the faces
        bool haveOld_;

        //- True after two complete passes; enables the fixed-step stencil
        bool haveOldOld_;

        //-
        scalar prevt_;

        //-
        scalar deltat_;

        //-
        List<T> new_;

        //-
        List<T> old_;

        //-
        List<T> oldOld_;

    public:

        //-
        tdPrimitive(T initv, scalar tprev)
        :
            varTimeStep_(false),
            haveOld_(false),
            haveOldOld_(false),
            prevt_(tprev),
            deltat_(tprev),
            new_(1, initv),
            old_(0),
            oldOld_(0)
        {
        };

        //-
    
        tdPrimitive()
        :
            varTimeStep_(false),
            haveOld_(false),
            haveOldOld_(false),
            prevt_(0.0),
            deltat_(GREAT),
            new_(1, pTraits<T>::zero),
            old_(0),
            oldOld_(0)
        {
        };
        
        ~tdPrimitive()
        {
        };

        void resize(label newsz)
        {
            if (new_.size() != newsz)
            {
                new_.resize(newsz);
                haveOld_ = false;
                haveOldOld_ = false;
            }

            // One allocation for the whole face list. dot() must not grow
            // these one face at a time: List::resize copies the array.
            if (old_.size() != newsz)
            {
                old_.resize(newsz);
            }
            if (oldOld_.size() != newsz)
            {
                oldOld_.resize(newsz);
            }
        }
        
        void setTimeStepVariable(bool adjustTimeStep)
        {
            varTimeStep_ = adjustTimeStep;
        }

        label size() const
        {
            return new_.size();
        }

        const T& value(label idx) const
        {
            return new_[idx];
        }

        T& value(label idx)
        {
            return new_[idx];
        }

        T dot(scalar newt, label idx)
        {
            if (old_.size() != new_.size() || oldOld_.size() != new_.size())
            {
                old_.resize(new_.size());
                oldOld_.resize(new_.size());
                haveOld_ = false;
                haveOldOld_ = false;
            }

            if (newt > prevt_)
            {
                deltat_ = newt - prevt_;
                prevt_ = newt;
            }

            T dvdt = pTraits<T>::zero;

            if (!haveOld_)
            {
                dvdt = new_[idx] / deltat_;
            }
            else if (varTimeStep_ || !haveOldOld_)
            {
                dvdt = (new_[idx] - old_[idx]) / deltat_;
            }
            else
            {
                dvdt =
                    (3.0*new_[idx] - 4.0*old_[idx] + oldOld_[idx])
                    / 2.0 / deltat_;
            }

            oldOld_[idx] = old_[idx];
            old_[idx] = new_[idx];

            if (idx == new_.size() - 1)
            {
                if (haveOld_)
                {
                    haveOldOld_ = true;
                }
                haveOld_ = true;
            }

            return dvdt;
        }

    };

        //- Runtime type information
        TypeName("AcousticAnalogy");

        // Constructors

        //- Construct from given Time and dictionary.
        AcousticAnalogy
        (
            const word& name,
            const Time&,
            const dictionary&
        );

        //- Construct from given objectRegistry and dictionary.
        AcousticAnalogy
        (
            const word& name,
            const objectRegistry&,
            const dictionary&
        );

        //- Destructor
        virtual ~AcousticAnalogy();

public:
        //- Get observers list
        List<SoundObserver> getSoundObservers();

        //- Read the Curle data
        virtual bool read(const dictionary&);

        //-
        virtual bool write();

        //-
        virtual bool execute();
};

} //namespace functionObjects

} //namespace Foam

#endif

//
//END-OF-FILE
//



