#ifndef __REMOLLGLOBALFIELD_HH
#define __REMOLLGLOBALFIELD_HH

/*!
   \class remollGlobalField
   \brief Global field interface
*/

#include "G4MagneticField.hh"

#include <vector>

class G4GenericMessenger;
class G4EquationOfMotion;
class G4MagIntegratorStepper;
class G4ChordFinder;
class G4FieldManager;
class G4PropagatorInField;

class remollMagneticField;
class remollDetectorConstruction;

class remollGlobalField : public G4MagneticField {
    private:
        static G4ThreadLocal remollGlobalField* fObject;

    public:
        /// GetObject() returns the single remollGlobalField object.
        /// It is constructed, if necessary.
        static remollGlobalField* GetObject(remollDetectorConstruction* const);
        static remollGlobalField* GetObject();

    private:
        remollGlobalField(remollDetectorConstruction* const);
        remollGlobalField(const remollGlobalField&);
        remollGlobalField& operator=(const remollGlobalField&);

    public:
	virtual ~remollGlobalField();

	void AddNewField(G4String& name);
	void SetFieldScaleByString(G4String& name_scale);
	void SetFieldScale(const G4String& name, G4double scale);
        void SetMagnetCurrentByString(G4String& name_scale);
	void SetMagnetCurrent(const G4String& name, G4double scale);

	void GetFieldValue(const G4double[], G4double*) const;

    private:
        /// Set the stepper
        void SetEquation();

        /// Set the stepper
        void SetStepper();

        /// Set the stepper
        void SetChordFinder();

    public:
        /// Set the stepper types
        void SetEquationType(G4int i) { fEquationType = i; SetEquation(); }

        /// Set the stepper types
        void SetStepperType(G4int i) { fStepperType = i; SetStepper(); }

        /// Set the minimum step length
        void SetMinStep(G4double stp) { fMinStep = stp; SetAccuracyParameters(); }

        /// Set the delta chord length
        void SetDeltaChord(G4double dcr) { fDeltaChord = dcr; SetAccuracyParameters(); }

        /// Set the delta one step length
        void SetDeltaOneStep(G4double stp) { fDeltaOneStep = stp; SetAccuracyParameters(); }

        /// Set the delta intersection length
        void SetDeltaIntersection(G4double its) { fDeltaIntersection = its; SetAccuracyParameters(); }

        /// Set the minimum eps length
        void SetEpsMin(G4double eps) { fEpsMin = eps; SetAccuracyParameters(); }

        /// Set the maximum eps length
        void SetEpsMax(G4double eps) { fEpsMax = eps; SetAccuracyParameters(); }

        /// Set the accuracy parameters
        void SetAccuracyParameters();

        /// Print the accuracy parameters
        void PrintAccuracyParameters();

    private:
        G4int fEquationType;
        G4int fStepperType;

        G4double fMinStep;
        G4double fDeltaChord;
        G4double fDeltaOneStep;
        G4double fDeltaIntersection;
        G4double fEpsMin;
        G4double fEpsMax;

        G4EquationOfMotion*     fEquation;
	G4int                   fEquationDoF;

        G4FieldManager*         fFieldManager;
        G4PropagatorInField*    fFieldPropagator;
        G4MagIntegratorStepper* fStepper;
        G4ChordFinder*          fChordFinder;

    private:
	std::vector<remollMagneticField*> fFields;

        G4int fNfp;
        G4bool fFirst;
        const remollMagneticField **fFieldp;

	remollMagneticField* GetFieldByName(const G4String& name);

	G4GenericMessenger* fMessenger;
	G4GenericMessenger* fGlobalFieldMessenger;

        remollDetectorConstruction* fDetectorConstruction;
};


#endif//__REMOLLGLOBALFIELD_HH
