#ifndef PLUMED_kenref_KEnRefBias_h
#define PLUMED_kenref_KEnRefBias_h

#include "core/ActionAtomistic.h"
#include "bias/Bias.h"
#include "config.h"
#include "core/KEnRef.h"
#include "core/Table.h"
#include <chrono>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace PLMD::kenref {
    class KEnRefBias : public bias::Bias, public ActionAtomistic {
        // ------------------------------------------------------------------ //
        //  Energy model                                                        //
        // ------------------------------------------------------------------ //
        KEnRef<KEnRef_Real_t>::energyModel selectedEnergyModel =
                KEnRef<KEnRef_Real_t>::energyModel::UNKNOWN;

        // ------------------------------------------------------------------ //
        //  Common scalar parameters                                            //
        // ------------------------------------------------------------------ //
        KEnRef_Real_t k_ = 1.0;
        KEnRef_Real_t n_ = 0.25;
        KEnRef_Real_t maxForce_ = 999.0;
        KEnRef_Real_t maxForceSquared_ = maxForce_ * maxForce_;

        // ------------------------------------------------------------------ //
        //  SIGMA model parameters                                              //
        // ------------------------------------------------------------------ //
        KEnRef_Real_t proton_mhz_ = 600.0;
        NamedRowVector<KEnRef_Real_t> rates_ = Table( //TODO provide a way to change it
            {{"5.0e+08", "2.5e+08", "1.0e+12", "1.0e+04"}},
            {{"kens", "kc", "kmethyl", "karo"}}
        ).toNamedRowVector<KEnRef_Real_t>();

        std::vector<SpecDenData<KEnRef_Real_t> > spec_den_data_list_;
        std::string exp_data_folder_; // EXP_DATA_FOLDER keyword value

        // ------------------------------------------------------------------ //
        //  PLATEAUS model parameters                                           //
        // ------------------------------------------------------------------ //
        Eigen::MatrixX<KEnRef_Real_t> g0_;
        Table experimental_data_table_;
        std::string experimental_data_file_; // EXP_DATA_FILE keyword value

        // ------------------------------------------------------------------ //
        //  Atom lists (PLUMED AtomNumber)                                      //
        // ------------------------------------------------------------------ //
        std::vector<AtomNumber> atoms_; // full list passed to requestAtoms()
        std::vector<AtomNumber> guideAtoms_; // alignment atoms
        std::vector<AtomNumber> subAtoms_; // restrained H-pair atoms

        // ------------------------------------------------------------------ //
        //  Atom-name / index mappings                                          //
        // ------------------------------------------------------------------ //
        std::string atomname_mapping_file_; // ATOMNAME_MAPPING keyword value

        // atom name (normalised) → 1-based serial from the mapping PDB
        std::map<std::string, int> atomName_to_globalSerial_map_;

        // atom name → compact sub-index (0-based)
        std::map<std::string, int> atomName_to_sub0Id_map_;

        // compact sub-index → 1-based serial
        std::vector<int> sub0Id_to_global1Serial_;

        // 1-based serial → local index in the positions vector from PLUMED
        std::map<int, int> serial_to_localIdx_;

        // atom-pair index list (sub0Id based)
        std::vector<std::tuple<int, int> > atomId_pairs_;

        // atom-name pair list (strings)
        std::vector<std::tuple<std::string, std::string> > atomName_pairs_;

        // ------------------------------------------------------------------ //
        //  Reference structure (for Kabsch alignment)                          //
        // ------------------------------------------------------------------ //
        std::string reference_pdb_;
        CoordsMatrixType<KEnRef_Real_t> guideAtomsReferenceCoords_;
        CoordsMatrixType<KEnRef_Real_t> guideAtomsReferenceCoordsCentered_;

        // ------------------------------------------------------------------ //
        //  Per-step (Current ?) coordinate state                                           //
        // ------------------------------------------------------------------ //
        CoordsMatrixType<KEnRef_Real_t> subAtomsX_;
        CoordsMatrixType<KEnRef_Real_t> lastFrameSubAtomsX_;
        CoordsMatrixType<KEnRef_Real_t> lastFrameGuideAtomsX_;

        // ------------------------------------------------------------------ //
        //  Derivatives buffers (pre-allocated once in initializeParameters)  //
        //  – allDerivatives_buffer_: holds derivatives for ALL replicas      //
        //    (size = numSubAtoms * 3 * numSimulations)                       //
        //  – derivatives_buffer_:   THIS replica's derivatives               //
        //    (size = numSubAtoms * 3); CoordsMapType is placed over this     //
        //    in calculate() so the inverse-Kabsch transform is zero-copy.   //
        // ------------------------------------------------------------------ //
        std::vector<KEnRef_Real_t> allDerivatives_buffer_;
        std::vector<KEnRef_Real_t> derivatives_buffer_;

        // ------------------------------------------------------------------ //
        //  Grouping for PLATEAUS model                                         //
        // ------------------------------------------------------------------ //
        std::vector<std::vector<std::vector<int> > > simulated_grouping_list_;

        // ------------------------------------------------------------------ //
        //  Control flags                                                       //
        // ------------------------------------------------------------------ //
        bool paramsInitialized = false; // use instead of (step==0): handles restarts because the simulation may be continuing after 0
        bool fit_to_reference_ = true;
        bool saturate_forces_ = false;

        // ------------------------------------------------------------------ //
        //  Timing                                                              //
        // ------------------------------------------------------------------ //
        long long calculateForces_time = 0; //TODO fix this

        // ------------------------------------------------------------------ //
        //  Private helpers                                                     //
        // ------------------------------------------------------------------ //
        /** One-time setup (mirrors fillParamsStep0). Called on the first step. */
        void initializeParameters();

        /** Fill `out` with current sub-atom positions in Å (mirrors fillSubAtomsX). */
        void fillSubAtomsX(CoordsMatrixType<KEnRef_Real_t> &out) const;

        /** Return current guide-atom positions in Å (mirrors getGuideAtomsX). */
        CoordsMatrixType<KEnRef_Real_t> getGuideAtomsX() const;

        /**
         * Periodic-image correction so atoms don't jump across box boundaries
         * (mirrors KEnRefForceProvider::restoreNoJump).
         * @param box_nm   PLUMED box Tensor in nm; scaled internally by toAngstrom.
         */
        void restoreNoJump(CoordsMatrixType<KEnRef_Real_t> &atoms, const CoordsMatrixType<KEnRef_Real_t> &reference,
            const Tensor &box_nm, bool toAngstrom, int numOmpThreads, const bool printStatistics) const;

    public:
        //May remove this macro later
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        explicit KEnRefBias(const ActionOptions &);
        static void registerKeywords(Keywords &);
        void calculate() override;

        // Explicit overrides to resolve multiple-inheritance ambiguity
        void lockRequests() override;
        void unlockRequests() override;
        void calculateNumericalDerivatives(ActionWithValue *a) override;
    };
} // namespace PLMD::kenref

#endif // PLUMED_kenref_KEnRefBias_h
