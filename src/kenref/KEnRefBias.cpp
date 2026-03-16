#include "bias/Bias.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/Vector.h"
#include "tools/Tensor.h"
#include "tools/Communicator.h"
#include "tools/PDB.h"
// Include your KEnRef headers
#include "core/KEnRef.h"
#include "core/kabsch.h"   // Kabsch_Umeyama
#include "core/IoUtils.h"  // IoUtils::getAtomMappingFromPdb, load_spec_den_data, readTable, normalizeName …
#include "KEnRefBias.h"

#include <unistd.h>

#include "tools/OpenMP.h"

namespace PLMD::kenref {
    // Register the action
    PLUMED_REGISTER_ACTION(KEnRefBias, "KENREF")

    // ============================================================
    //  registerKeywords
    // ============================================================
    void KEnRefBias::registerKeywords(Keywords &keys) {
        //First call the parent
        Bias::registerKeywords(keys);
        ActionAtomistic::registerKeywords(keys);

        keys.setDisplayName("KENREF");
        // keys.use("ARG");  // Commented out - ARG not needed for this bias
        // 2nd: Add your keywords here
        // Required parameters
        keys.add("compulsory", "MODEL", "SIGMA", "The energy model to use (SIGMA or PLATEAUS)");
        keys.add("compulsory", "K", "1.0", "Force constant");
        keys.add("compulsory", "N", "0.25", "Power scaling factor");

        // SIGMA model
        keys.add("optional", "EXP_DATA_FOLDER", "Folder with experimental data (SIGMA model)");
        keys.add("optional", "PROTON_MHZ", "Proton MHz (SIGMA model)");
        keys.add("optional", "RATES", "Relaxation rates, comma-separated (SIGMA model)");
        // keys.add("optional", "SPEC_DEN_DATA_FILE", "Spectral-density data file (SIGMA model)");

        // PLATEAUS model
        keys.add("optional", "EXP_DATA_FILE", "Experimental data file (PLATEAUS model)");
        keys.add("optional", "G0", "Target g values (PLATEAUS model)"); //TODO is this  right?

        // Atom selection and reference
        keys.add("atoms", "GUIDE_ATOMS", "Atoms used for alignment to reference");
        keys.add("optional", "REF", "Reference structure PDB for alignment");
        keys.add("optional", "ATOMNAME_MAPPING", "PDB file with atom-name → atom-index mapping");
        // keys.add("optional", "INDEX", "Index file (created using 'gmx make_ndx' containing atom groups");
        // keys.add("optional", "GUIDE_ATOMS_GROUP_NAME", "Name of guide group in index file to use for alignment");

        // Advanced
        keys.add("optional", "MAX_FORCE", "Maximum force magnitude (default 999)");
        // keys.add("optional", "STRIDE", "frequency of application of KEnRef"); //NO NEED
        keys.addFlag("FIT_TO_REFERENCE", false, "Fit coordinates to reference before calculating energy");
        keys.addFlag("SATURATE_FORCES", false, "Clamp forces to MAX_FORCE");

        // Output components
        // keys.addOutputComponent("bias", "default", "The instantaneous value of the bias potential"); //already added by Bias
        keys.addOutputComponent("energy", "default", "Total KEnRef restraint energy");
        keys.addOutputComponent("rmsd", "default", "RMSD from reference (after fitting)");

        // more keywords? Claude says NO
        keys.use("RESTART");
        keys.use("UPDATE_FROM");
        keys.use("UPDATE_UNTIL");
    }

    // ============================================================
    //  Constructor
    // ============================================================
    KEnRefBias::KEnRefBias(const ActionOptions &ao) : PLUMED_BIAS_INIT(ao), ActionAtomistic(ao) {

        // volatile bool holdToDebug = true;
        // while (/*simulationIndex > 0 &&*/ holdToDebug) {
        //     sleep(1);
        // }

        // Parse your keywords here

        // --- Energy model ---
        std::string energyModelStr;
        parse("MODEL", energyModelStr);
        const auto it = KEnRef<KEnRef_Real_t>::energyModelMap.find(energyModelStr);
        if (it == KEnRef<KEnRef_Real_t>::energyModelMap.end() ||
            it->second == KEnRef<KEnRef_Real_t>::energyModel::UNKNOWN)
            error("MODEL must be SIGMA or PLATEAUS");
        selectedEnergyModel = it->second;

        // --- Common parameters ---
        parse("K", k_);
        parse("N", n_);
        KEnRef_Real_t maxForce = maxForce_;
        parse("MAX_FORCE", maxForce);
        maxForce_ = maxForce;
        maxForceSquared_ = maxForce * maxForce;

        parseFlag("FIT_TO_REFERENCE", fit_to_reference_);
        parseFlag("SATURATE_FORCES", saturate_forces_);

        // --- Model-specific parameters ---
        if (selectedEnergyModel == KEnRef<KEnRef_Real_t>::energyModel::SIGMA) {
            parse("PROTON_MHZ", proton_mhz_);
            parse("EXP_DATA_FOLDER", exp_data_folder_);
            // TODO: parse custom RATES if provided
        } else {
            parse("EXP_DATA_FILE", experimental_data_file_);
        }

        // parseArgumentList("ARG");// use it if you will take input from other actions

        // --- Reference PDB and atom-name mapping ---
        // Parse ATOMNAME_MAPPING first: it is required to build the atom-name → serial
        // index map, which is needed NOW (in the constructor) to derive the sub-atom list.
        parse("ATOMNAME_MAPPING", atomname_mapping_file_);
        if (atomname_mapping_file_.empty())
            error("ATOMNAME_MAPPING is required");

        // REF is optional. If omitted, the ATOMNAME_MAPPING PDB doubles as the
        // reference structure for Kabsch alignment (it already contains coordinates).
        parse("REF", reference_pdb_);
        if (reference_pdb_.empty()) {
            reference_pdb_ = atomname_mapping_file_;
            log.printf("  REF not specified – using ATOMNAME_MAPPING file as reference: %s\n",
                       reference_pdb_.c_str());
        }

        // ---- Load atom-name → serial map (1-based) -----------------------------------
        // We need this here (not deferred to initializeParameters) because we must call
        // requestAtoms() before checkRead(), and the atom list is derived from the data files.
        atomName_to_globalSerial_map_ = IoUtils::getAtomMappingFromPdb<std::string, int>(
            atomname_mapping_file_, IoUtils::fill_atomId_to_index_Map);
        if (atomName_to_globalSerial_map_.empty())
            error("No atom mapping found in ATOMNAME_MAPPING file: " + atomname_mapping_file_);

        // ---- Load experimental data and extract atom-name pairs ----------------------
        // Doing this here lets us derive the sub-atom list without the user having to
        // specify ATOMS explicitly. initializeParameters() will use the already-loaded
        // data to build the sub-indexing structures.
        atomName_pairs_.clear();
        if (selectedEnergyModel == KEnRef<KEnRef_Real_t>::energyModel::SIGMA) {
            const bool handleNames = IoUtils::should_handleNames(atomName_to_globalSerial_map_);
            // Mirrors: SpecDenData_ = IoUtils::load_spec_den_data(Settings::experimentalDataFolder, handleNames)
            spec_den_data_list_ = IoUtils::load_spec_den_data(exp_data_folder_, handleNames);
            if (spec_den_data_list_.empty())
                error("No spectral-density data loaded from EXP_DATA_FOLDER: " + exp_data_folder_);

            // Collect all atom-name pairs across all SpecDenData objects
            for (const auto &sdd : spec_den_data_list_)
                for (const auto &ap : sdd.get_atom_pairs())
                    atomName_pairs_.emplace_back(ap);
        } else { // PLATEAUS
            // Mirrors: experimentalData_table_ = IoUtils::readTable(Settings::experimentalDataFileName, ...)
            experimental_data_table_ =
                    IoUtils::readTable(experimental_data_file_, true, false, "\\s*,\\s*", -1 /*no limit*/);
            if (experimental_data_table_.rowCount() == 0)
                error("No experimental data found in EXP_DATA_FILE: " + experimental_data_file_);
            const Table &table = experimental_data_table_;
            bool handleUnpreparedAtomNames = false;
            for (int row = 0; row < table.rowCount(); row++) {
                if (IoUtils::isNotPrepared(table(row, "atom1")) ||
                    IoUtils::isNotPrepared(table(row, "atom2"))) {
                    log.printf("  WARNING: It seems that your data is from an unprepared file. We will try to handle it, but we can not guarantee the results.\n");
                    handleUnpreparedAtomNames = true;
                    break;
                }
            }
            for (int row = 0; row < table.rowCount(); row++) {
                atomName_pairs_.emplace_back(
                    IoUtils::normalizeName(table(row, "atom1"), handleUnpreparedAtomNames),
                    IoUtils::normalizeName(table(row, "atom2"), handleUnpreparedAtomNames));
            }
            // Fill g0_ matrix (numPairs × 2)
            g0_.resize(table.rowCount(), 2);
            for (int i = 0; i < table.rowCount(); ++i) {
                std::istringstream t1(table(i, "g1")), t2(table(i, "g2"));
                t1 >> g0_(i, 0);
                t2 >> g0_(i, 1);
            }
        }

        // ---- Derive subAtoms_ from atom-name pairs (no ATOMS keyword needed) ---------
        // Collect the unique 1-based serials that appear in any pair, in sorted order
        // (std::set gives us uniqueness + ordering automatically).
        std::set<int> subSerials;
        for (const auto &[a1, a2] : atomName_pairs_) {
            subSerials.insert(atomName_to_globalSerial_map_.at(a1));
            subSerials.insert(atomName_to_globalSerial_map_.at(a2));
        }
        subAtoms_.clear();
        for (int s : subSerials) {
            // AtomNumber an;
            // an.setSerial(s);
            // subAtoms_.emplace_back(an);
            subAtoms_.emplace_back(AtomNumber().setSerial(s)); //TODO validate this
        }
        if (subAtoms_.empty())
            error("No restrained atoms derived from experimental data – check your data files and ATOMNAME_MAPPING");

        // --- Atom lists ---
        // Guide atoms are used for Kabsch alignment only.
        parseAtomList("GUIDE_ATOMS", guideAtoms_);

        // Merge guide + sub atoms into the single atoms_ list that PLUMED will fetch.
        // initializeParameters() will sort out the sub-indexing within this merged list.
        atoms_.insert(atoms_.end(), guideAtoms_.begin(), guideAtoms_.end());
        atoms_.insert(atoms_.end(), subAtoms_.begin(), subAtoms_.end());
        requestAtoms(atoms_);

        // --- Output components ---
        // addComponent("bias"); //already added by Bias
        // componentIsNotPeriodic("bias");
        addComponent("energy");
        componentIsNotPeriodic("energy");
        addComponent("rmsd");
        componentIsNotPeriodic("rmsd");

        checkRead();

        log.printf("  KEnRef bias  model=%s  k=%f  n=%f  proton_mhz=%f\n", energyModelStr.c_str(), k_, n_, proton_mhz_);
        log.printf("  %zu guide atoms, %zu restrained atoms\n", guideAtoms_.size(), subAtoms_.size());
        // log.printf("  %zu relaxation rates\n", rates_.size());

        cite("Restraining interproton angular and distance dynamics with KEnRef. "
            "Amr Alhossary & Colin Smith, J Phys Chem B. 130, 11 (2026) DOI: 10.1021/acs.jpcb.5c08554");
    }

    // ============================================================
    //  initializeParameters   (≈ fillParamsStep0)
    //
    //  Called once on the first step (paramsInitialized == false).
    //  Mirrors KEnRefForceProvider::fillParamsStep0() but uses
    //  PLUMED atom numbering instead of GROMACS domain-decomp IDs.
    // ============================================================
    void KEnRefBias::initializeParameters() {
        log.printf("KEnRefBias::initializeParameters()\n - at log");
        log.printf("MPI DEBUG: multi_sim_comm size = %zu rank = %zu\n", multi_sim_comm.Get_size(), multi_sim_comm.Get_rank());

        // NOTE: atomName_to_globalSerial_map_, atomName_pairs_, spec_den_data_list_,
        // experimental_data_table_, and g0_ are already populated by the constructor.
        // They were loaded there so that requestAtoms() could be called with the correct
        // sub-atom list before checkRead().  Nothing to reload here.

        // ---- 1. Build sub-atom index structures --------------------------------------
        //   Mirrors the globalAtomIdFlags / sub0Id_to_global1Id block in fillParamsStep0.
        //
        //   In PLUMED there is no domain decomposition to worry about.
        //   We work entirely with PLUMED's local atom list (atoms_ requested above).
        //   We build:
        //     global1Serial_to_sub0Id_  :  serial (1-based) → compact index
        //     sub0Id_to_global1Serial_  :  compact index    → serial (1-based)
        //
        //   "global serial" == AtomNumber::serial() which is 1-based in PLUMED.

        int maxSerialOfInterest = -1;
        // First pass: find the highest serial we care about
        for (const auto &[a1, a2]: atomName_pairs_) {
            int s1 = atomName_to_globalSerial_map_.at(a1);
            int s2 = atomName_to_globalSerial_map_.at(a2);
            maxSerialOfInterest = std::max(maxSerialOfInterest, std::max(s1, s2));
        }

        std::vector serialFlags(maxSerialOfInterest + 1, false);
        for (const auto &[a1, a2]: atomName_pairs_) {
            serialFlags[atomName_to_globalSerial_map_.at(a1)] = true;
            serialFlags[atomName_to_globalSerial_map_.at(a2)] = true;
        }

        std::vector global1Serial_to_sub0Id(maxSerialOfInterest + 1, -1);
        sub0Id_to_global1Serial_.clear();
        int subId = 0;
        for (int s = 0; s <= maxSerialOfInterest; ++s) {
            if (serialFlags[s]) {
                global1Serial_to_sub0Id[s] = subId;
                sub0Id_to_global1Serial_.push_back(s);
                subId++;
            }
        }
        const int numSubAtoms = subId;

        // atomName → sub0Id  (mirrors atomName_to_atomSub0Id_map_)
        atomName_to_sub0Id_map_.clear();
        for (const auto &[name, serial]: atomName_to_globalSerial_map_) {
            if (serial <= maxSerialOfInterest && global1Serial_to_sub0Id[serial] != -1)
                atomName_to_sub0Id_map_[name] = global1Serial_to_sub0Id[serial];
        }

        // Build integer atom-pair index list  (mirrors atomId_pairs_)
        atomId_pairs_ = KEnRef<KEnRef_Real_t>::atomNamePairs_2_atomIdPairs(
            atomName_pairs_, atomName_to_sub0Id_map_);

        // Cache per-SpecDenData atomId pair sub-indices (SIGMA only)
        if (selectedEnergyModel == KEnRef<KEnRef_Real_t>::energyModel::SIGMA) {
            for (auto &sdd: spec_den_data_list_) {
                const auto &ap = sdd.get_atom_pairs();
                sdd.set_atomIdPairs_to_sub0Atom_id_pairs_cache(
                    {KEnRef<KEnRef_Real_t>::atomNamePairs_2_atomIdPairs(ap, atomName_to_sub0Id_map_)});
            }
        }

        // ---- 2. PLUMED serial → local index in atoms_ --------------------------------
        //   PLUMED delivers positions in the order of atoms_ (the list we gave requestAtoms).
        //   We need to map each sub0Id back to the corresponding position in atoms_.
        //   atoms_ = [guideAtoms_ | subAtoms_] in the order we appended them.
        serial_to_localIdx_.clear();
        for (int li = 0; li < static_cast<int>(atoms_.size()); ++li)
            serial_to_localIdx_[atoms_[li].serial()] = li;

        // ---- 3. grouping list for PLATEAUS -------------------------------------------
        //   Mirrors the switch(numSimulations) block.
        //   PLUMED does not have multi-simulation MPI the same way GROMACS does.
        //   We obtain the replica count from the Communicator; fall back to 1.
        int numSimulations = 1;
        if (Communicator::plumedHasMPI() && Communicator::initialized())
            numSimulations = multi_sim_comm.Get_size();

        if (selectedEnergyModel == KEnRef<KEnRef_Real_t>::energyModel::PLATEAUS) {
            switch (numSimulations) {
                case 1:
                    simulated_grouping_list_ = {{{0}}, {{0}}};
                    break;
                case 2:
                    simulated_grouping_list_ = {{{0, 1}}, {{0}, {1}}};
                    break;
                case 3:
                    simulated_grouping_list_ = {{{0, 1, 2}}, {{0}, {1}, {2}}};
                    break;
                default:
                    error("PLATEAUS model does not yet support more than 3 replicas yet.");
            }
        }

        // ---- 4. Allocate coordinate and derivative buffers ---------------------------
        // All three buffers are allocated ONCE here so that calculate() never performs
        // heap allocations on the hot path.
        //
        // Memory layout (mirrors KEnRefForceProvider):
        //   allDerivatives_buffer_  – flat buffer holding derivatives for ALL replicas
        //                             (size = subSize * numSimulations).
        //   derivatives_buffer_     – flat buffer for THIS replica's derivatives
        //                             (size = subSize).
        //   For a single simulation both have the same size (numSimulations == 1).
        //   For multi-sim, allDerivatives_buffer_ is filled on rank 0 and scattered;
        //   derivatives_buffer_ is the scatter destination on every rank.
        //
        // The Eigen CoordsMapType in calculate() creates a zero-copy view directly
        // over derivatives_buffer_.data(), so no extra copy is needed before the
        // inverse-Kabsch transform.
        subAtomsX_.resize(numSubAtoms, 3);
        lastFrameSubAtomsX_.resize(subAtomsX_.rows(), subAtomsX_.cols());
        const int subSize = subAtomsX_.rows() * subAtomsX_.cols();
        allDerivatives_buffer_.resize(static_cast<size_t>(subSize) * numSimulations, 0);
        derivatives_buffer_.resize(subSize, 0); //TODO No need to resize it if ! isMultiSim

        //TODO I want the program to be high performance.
        //  1) allocate a new buffer every cycle then copy the data between memory locations?
        //  2) Why don't we use the same buffer for both the Eigen Matrix and the send or receive buffer?

        // ---- 5. Load reference structure and initialise Kabsch state ------------------
        if (fit_to_reference_ && !reference_pdb_.empty()) {
            // Read reference coords for guide atoms from PDB.
            // readAtomsFromPDB() fills the positions in PLUMED's internal units (nm).
            // KEnRef works in Ångström; we scale by 10 when filling subAtomsX_.
            PDB pdb;
            if (!pdb.read(reference_pdb_, usingNaturalUnits(), 0.1))
                error("Cannot read reference PDB: " + reference_pdb_);

            const int nGuide = static_cast<int>(guideAtoms_.size());
            guideAtomsReferenceCoords_.resize(nGuide, 3);
            for (int i = 0; i < nGuide; ++i) {
                Vector pos = pdb.getPosition(guideAtoms_[i]);
                // pdb positions are in nm (PLUMED internal); convert to Å
                if constexpr (std::is_same_v<KEnRef_Real_t, double>) {
                    guideAtomsReferenceCoords_(i, 0) = pos[0] * 10;
                    guideAtomsReferenceCoords_(i, 1) = pos[1] * 10;
                    guideAtomsReferenceCoords_(i, 2) = pos[2] * 10;
                }else {
                    guideAtomsReferenceCoords_(i, 0) = static_cast<KEnRef_Real_t>(pos[0]) * 10;
                    guideAtomsReferenceCoords_(i, 1) = static_cast<KEnRef_Real_t>(pos[1]) * 10;
                    guideAtomsReferenceCoords_(i, 2) = static_cast<KEnRef_Real_t>(pos[2]) * 10;
                }
            }
            // Pre-center the reference (mirrors setGuideAtomsReferenceCoords)
            guideAtomsReferenceCoordsCentered_ =
                    Kabsch_Umeyama<KEnRef_Real_t>::translateCenterOfMassToOrigin(guideAtomsReferenceCoords_);

            // Initialise lastFrame buffers with the reference (matches fillParamsStep0 end)
            lastFrameGuideAtomsX_.resize(nGuide, 3);
            lastFrameGuideAtomsX_ = guideAtomsReferenceCoords_;
        }

        // Seed lastFrameSubAtomsX_ with the *current* positions.
        // We do this here so restoreNoJump has a valid reference on step 0.
        // (Mirrors fillParamsStep0's final fillSubAtomsX / lastFrameGuideAtomsX_ assignment.)
        fillSubAtomsX(lastFrameSubAtomsX_);

        log.printf("  initializeParameters done:  %d sub-atoms, %zu atom pairs\n", numSubAtoms, atomId_pairs_.size());
        std::cout << "  initializeParameters done:  "<< numSubAtoms <<" sub-atoms, "<< atomId_pairs_.size() <<" atom pairs\n";
    }

    // ============================================================
    //  Helpers
    // ============================================================

    // Fill subAtomsX from the current PLUMED positions.
    // Mirrors KEnRefForceProvider::fillSubAtomsX() (without GROMACS domain-decomp).
    // Positions in PLUMED are in nm; we convert to Å (* 10).
    void KEnRefBias::fillSubAtomsX(CoordsMatrixType<KEnRef_Real_t> &out) const {
        const std::vector<Vector> &pos = getPositions();
        for (int i = 0; i < out.rows(); ++i) {
            int serial = sub0Id_to_global1Serial_[i]; // 1-based serial
            const int li = serial_to_localIdx_.at(serial); // index into pos[]
            out(i, 0) = static_cast<KEnRef_Real_t>(pos[li][0]) * 10; // nm → Å
            out(i, 1) = static_cast<KEnRef_Real_t>(pos[li][1]) * 10;
            out(i, 2) = static_cast<KEnRef_Real_t>(pos[li][2]) * 10;
        }
    }

    // Fill guide-atom positions from the current PLUMED positions.
    // Mirrors KEnRefForceProvider::getGuideAtomsX().
    CoordsMatrixType<KEnRef_Real_t> KEnRefBias::getGuideAtomsX() const {
        const std::vector<Vector> &pos = getPositions();
        const int nGuide = static_cast<int>(guideAtoms_.size());
        CoordsMatrixType<KEnRef_Real_t> gx(nGuide, 3);
        for (int i = 0; i < nGuide; ++i) {
            const int li = serial_to_localIdx_.at(guideAtoms_[i].serial());
            gx(i, 0) = static_cast<KEnRef_Real_t>(pos[li][0]) * 10;
            gx(i, 1) = static_cast<KEnRef_Real_t>(pos[li][1]) * 10;
            gx(i, 2) = static_cast<KEnRef_Real_t>(pos[li][2]) * 10;
        }
        return gx;
    }

    // restoreNoJump – mirrors KEnRefForceProvider::restoreNoJump()
    // but accepts a PLUMED Tensor (row-major, nm) instead of a GROMACS `matrix`.
    // We pass in the box already converted to Å (scale_factor applied by caller).
    void KEnRefBias::restoreNoJump(CoordsMatrixType<KEnRef_Real_t> &atoms, const CoordsMatrixType<KEnRef_Real_t> &reference,
        const Tensor &box_nm, const bool toAngstrom, int numOmpThreads, const bool printStatistics) const {

        // Convert box to the same unit as the coordinates
        const KEnRef_Real_t scale = toAngstrom ? 10.0 : 1.0;

        // Build a local 3×3 box array in the target unit (mirrors `msmul(box_, scale_factor, box)`).
        KEnRef_Real_t box[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                box[i][j] = static_cast<KEnRef_Real_t>(box_nm[i][j]) * scale;

        // Half the diagonal: used as the threshold for detecting a box-crossing.
        // (Matches the original: only diagonal elements are used for the threshold check.)
        const Eigen::RowVector3<KEnRef_Real_t> box_half(
            box[0][0] * 0.5,
            box[1][1] * 0.5,
            box[2][2] * 0.5);

        // Precompute active (non-zero) dimensions to avoid branching in inner loop .
        // Iterate from highest to lowest dimension (matches original DIM-1 → 0 order).
        std::vector<int> active_dims;
        for (int m = 2; m >= 0; --m){
            if (box_half[m] != 0){
                active_dims.push_back(m);
            }
        }

        int updated_count = 0;

        // Parallelize over atoms; each row is independent.
        // updated_count is a shared reduction (summed across threads).
#pragma omp parallel for num_threads(numOmpThreads) default(none) reduction(+:updated_count) \
    shared(atoms, reference, box, box_half, active_dims)
        for (int i = 0; i < atoms.rows(); ++i) {
            // Get row references once per atom
            auto atom = atoms.row(i); // Eigen row-expression: writes back to atoms
            const auto &refAtom = reference.row(i);
            bool atom_updated = false;

            // Triclinic algorithm (mirrors the original exactly):
            // When dimension m is corrected, ALL (ACTIVE) dimensions d ≤ m receive the
            // corresponding box vector component, which handles off-diagonal (shear) terms.
            for (const int m: active_dims) {
            // Check if atom jumped across the box in this dimension
                while (atom[m] - refAtom[m] <= -box_half[m]) {
                    // Jumped to negative image, correct by adding box size
                    for (int d = 0; d <= m; ++d)
                        atom[d] += box[m][d];
                    atom_updated = true;
                }
                while (atom[m] - refAtom[m] > box_half[m]) {
                    // Jumped to positive image, correct by subtracting box size
                    for (int d = 0; d <= m; ++d)
                        atom[d] -= box[m][d];
                    atom_updated = true;
                }
            }
            // Count atoms updated (not correction steps), matching the original's semantics.
            if (atom_updated) updated_count++;
        }

        if (printStatistics && updated_count > 0)
            log.printf("  INFO: restoreNoJump corrected %d atoms (out of %d)\n", updated_count, static_cast<int>(atoms.rows()));
    }

    // ============================================================
    //  calculate()   (≈ KEnRefForceProvider::calculateForces())
    // ============================================================
    void KEnRefBias::calculate() {
        const auto begin = std::chrono::high_resolution_clock::now();
        const long step = getStep();

        // ---- Step 0: one-time initialisation -----------------------------------------
        if (!paramsInitialized) {
            initializeParameters();
            paramsInitialized = true;
        }

        // ---- MPI / replica setup (mirrors the simulationContext_ / numSimulations block)
        const bool isMultiSim = Communicator::plumedHasMPI() && Communicator::initialized() && multi_sim_comm.Get_size() > 1;
        const int numSimulations = isMultiSim ? multi_sim_comm.Get_size() : 1;
        const int simulationIndex = isMultiSim ? multi_sim_comm.Get_rank() : 0;

        if (step % 10 == 0) {
            std::cout << "step:" << getStep() << std::endl;
            log.printf("  --> numSimulations %d  simulationIndex %d  step %ld\n", numSimulations, simulationIndex, step);
            log.flush();
        }

        // ---- Kabsch alignment --------------------------------------------------------
        //   Mirrors the "fit all models to reference" block in calculateForces().

        // 1) Get current guide-atom positions (Å)
        CoordsMatrixType<KEnRef_Real_t> guideAtomsX = getGuideAtomsX();

        // 2) NoJump correction on guide atoms
        const Tensor &box = getPbc().getBox(); // nm
        restoreNoJump(guideAtomsX, lastFrameGuideAtomsX_, box, /*toAngstrom=*/true,
            static_cast<int>(OpenMP::getNumThreads()), step % 10 == 0);

        // 3) Find the 3-D affine (Kabsch) transform aligning current guide atoms to reference
        const auto affine = Kabsch_Umeyama<KEnRef_Real_t>::find3DAffineTransform( guideAtomsX,
            guideAtomsReferenceCoordsCentered_, false, false, false);

        // 4) Fill subAtomsX with current positions (Å)
        fillSubAtomsX(subAtomsX_);

        // 5) NoJump correction on sub atoms
        restoreNoJump(subAtomsX_, lastFrameSubAtomsX_, box, /*toAngstrom=*/true,
            static_cast<int>(OpenMP::getNumThreads()), step % 10 == 0);

        // 6) Update "last frame" buffers
        //TODO Do the next 2 lines copy the data only or allocate **new** buffers as well?
        lastFrameSubAtomsX_ = subAtomsX_;
        lastFrameGuideAtomsX_ = guideAtomsX;

        // 7) Apply Kabsch rotation/translation to sub atoms
        subAtomsX_ = Kabsch_Umeyama<KEnRef_Real_t>::applyTransform(affine, subAtomsX_);

        // ---- Multi-simulation: gather all replicas' coordinates on rank 0 ------------
        //   Mirrors the MPI_Gather block.
        const int subSize = static_cast<int>(subAtomsX_.size()); // rows*3

        // allSimulationsSubAtomsX holds coords for ALL replicas (on rank 0).
        // For single sim it is just subAtomsX_ itself.
        CoordsMatrixType<KEnRef_Real_t> allSimulationsSubAtomsX; //TODO Save a buffer allSimulationsSubAtomsX_ and create a reference to reuse it every step.
        if (isMultiSim) {
            allSimulationsSubAtomsX.resize(numSimulations * subAtomsX_.rows(), 3);
            multi_sim_comm.Gather(subAtomsX_.data(), subSize, allSimulationsSubAtomsX.data(), subSize, /*root=*/0);
        } else {
            allSimulationsSubAtomsX = subAtomsX_;
        }

        // ---- Energy + derivative calculation (rank 0 only for multi-sim) -------------
        //   Mirrors the "if (!isMultiSimulation || simulationIndex == 0)" block.
        KEnRef_Real_t energy = 0.0;
        std::optional<std::vector<CoordsMatrixType<KEnRef_Real_t> > > all_derivatives_opt;

        if (!isMultiSim || simulationIndex == 0) {
            // Decompose the gathered matrix into a vector of per-replica matrices
            std::vector<CoordsMatrixType<KEnRef_Real_t> > allSims_vec;
            allSims_vec.reserve(numSimulations);
            for (int i = 0; i < numSimulations; ++i) {
                allSims_vec.emplace_back(
                    CoordsMapType<KEnRef_Real_t>(
                        &allSimulationsSubAtomsX.data()[i * subSize],
                        subAtomsX_.rows(), 3));
            }

            switch (selectedEnergyModel) {
                case KEnRef<KEnRef_Real_t>::energyModel::PLATEAUS:
                    std::tie(energy, all_derivatives_opt) =
                            KEnRef<KEnRef_Real_t>::coord_array_to_g_energy(
                                allSims_vec, atomId_pairs_,
                                simulated_grouping_list_, g0_,
                                k_, n_, /*gradient=*/true, static_cast<int>(OpenMP::getNumThreads()));
                    break;

                case KEnRef<KEnRef_Real_t>::energyModel::SIGMA:
                    // log.printf("  allSims_vec.size() = %zu\n", allSims_vec.size());
                    // for (size_t i = 0; i < allSims_vec.size(); ++i) {
                    //     log.printf("  allSims_vec[%zu]: rows=%d, cols=%d\n", i, (int)allSims_vec[i].rows(), (int)allSims_vec[i].cols());
                    // }
                    // log.printf("  spec_den_data_list_.size() = %zu\n", spec_den_data_list_.size());
                    // log.printf("  atomName_to_sub0Id_map_.size() = %zu\n", atomName_to_sub0Id_map_.size());
                    // log.flush();

                    std::tie(energy, all_derivatives_opt) =
                            KEnRef<KEnRef_Real_t>::coord_array_to_sigma_energy(
                                allSims_vec, rates_, spec_den_data_list_, proton_mhz_, k_, n_, atomName_to_sub0Id_map_,
                                /*gradient=*/true, static_cast<int>(OpenMP::getNumThreads()));
                    break;

                default:
                    error("Unknown energy model in calculate()");
            }

            if (step % 10 == 0)
                log.printf("  Step %ld  KEnRef energy = %g\n", step, energy);
        }

        // ---- Scatter derivatives back to each replica --------------------------------
        //   allDerivatives_buffer_ is pre-allocated in initializeParameters(); no heap
        //   allocation here.  Zero only the region we are about to write (rank 0 only,
        //   and only for the models actually returned).
        if (!isMultiSim || simulationIndex == 0) {
            // Pack the per-replica derivative matrices into the flat buffer
            const auto &dvec = all_derivatives_opt.value();
            for (int i = 0; i < static_cast<int>(dvec.size()); ++i)
                std::copy_n(dvec[i].data(), subSize, &allDerivatives_buffer_[i * subSize]);
        }

        // Pointer to the buffer to use for derivative mapping
        KEnRef_Real_t* deriv_data_ptr = nullptr;

        if (isMultiSim) {
            // Multi-simulation: scatter derivatives from root to all replicas
            multi_sim_comm.Scatter(allDerivatives_buffer_.data(), subSize, derivatives_buffer_.data(), subSize, /*root=*/0);
            deriv_data_ptr = derivatives_buffer_.data();
        } else {
            // TODO should do nothing: derivatives are already in allDerivatives_buffer_[0..subSize).
            // Single simulation: derivatives are already in allDerivatives_buffer_[0..subSize).
            // No copy needed - point directly to allDerivatives_buffer_
            deriv_data_ptr = allDerivatives_buffer_.data();
        }

        // ---- Transform derivatives back through the inverse Kabsch rotation ----------
        //   Mirrors: derivatives_rectified = Kabsch_Umeyama::applyInverseOfTransform(affine, map)
        const CoordsMapType<KEnRef_Real_t> deriv_map(deriv_data_ptr, subAtomsX_.rows(), 3);
        //TODO I am not sure The above line is the right way to use the pointer. The original line is KEnRefForceProvider.cpp:292

        CoordsMatrixType<KEnRef_Real_t> derivatives_rectified =
                Kabsch_Umeyama<KEnRef_Real_t>::applyInverseOfTransform(affine, deriv_map);

        // ---- Optional force saturation -----------------------------------------------
        if (saturate_forces_)
            KEnRef<KEnRef_Real_t>::saturate(derivatives_rectified, maxForceSquared_, OpenMP::getNumThreads());

        // ---- Apply forces to PLUMED atoms --------------------------------------------
        //   Mirrors: force[*piLocal] -= { dx, dy, dz }
        //   The gradient returned by KEnRef is dE/dx. We need to apply force = -dE/dx.
        //   Since we're using ActionAtomistic::addForce which adds forces directly,
        //   we pass -derivatives_rectified (negated).
        //
        //   NOTE: KEnRef returns derivatives in Å⁻¹ (energy / Å), but PLUMED works in
        //   nm. We scale by 1/10 to convert Å⁻¹ → nm⁻¹.
        constexpr KEnRef_Real_t toNm = 0.1;
        //TODO parallel for loop here
        for (int i = 0; i < derivatives_rectified.rows(); ++i) {
            //TODO is there a faster way than the next 2-step mapping? if so, use it everywhere
            int serial = sub0Id_to_global1Serial_[i];
            int li = serial_to_localIdx_.at(serial);
            Vector d(
                //Keep the static cast (although it looks not useful) for future need.
                static_cast<double>(derivatives_rectified(i, 0)) * toNm,
                static_cast<double>(derivatives_rectified(i, 1)) * toNm,
                static_cast<double>(derivatives_rectified(i, 2)) * toNm);
            // Get the atom index for this local index and add force (negated gradient)
            auto atomIndex = getValueIndices(getAbsoluteIndex(li));
            addForce(atomIndex, -d);
        }

        // Virial is set to zero (no explicit virial calculation – consistent with the
        // original GROMACS code which also ignores virial for this bias).
        // For Bias actions inheriting from ActionAtomistic, virial handling is different
        // from Colvar actions, and can be omitted if not using PBCs.

        // ---- RMSD (after fitting, Å → nm for PLUMED) ---------------------------------
        // TODO This measure is useless. I included it just to show something. Later, replace it with RMSDr or remove it.
        double rmsd = 0.0;
        if (fit_to_reference_ && guideAtomsReferenceCoords_.rows() > 0) {
            CoordsMatrixType<KEnRef_Real_t> diff =
                    guideAtomsX - guideAtomsReferenceCoords_;
            rmsd = std::sqrt(diff.rowwise().squaredNorm().mean());
            rmsd *= 0.1; // Å → nm
        }

        // ---- Set output components ---------------------------------------------------
        const double bias_value = energy;
        getPntrToComponent("bias")->set(bias_value);
        getPntrToComponent("energy")->set(bias_value);
        getPntrToComponent("rmsd")->set(rmsd);
        setBias(bias_value);

        // ---- Timing ------------------------------------------------------------------
        const auto end = std::chrono::high_resolution_clock::now();
        const auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        calculateForces_time += elapsed.count();
        if (step % 10 == 0 && simulationIndex == 0)
            log.printf("  This step (%ld): %.5f s  Total walltime: %.3f s\n",
                       step,
                       static_cast<double>(elapsed.count()) * 1e-9,
                       static_cast<double>(calculateForces_time) * 1e-9);
    }

    // ============================================================
    //  Ambiguity resolution
    // We need to lock both the Atoms (ActionAtomistic) and any Arguments/CVs (Bias)
    // ============================================================
    void KEnRefBias::lockRequests() {
        ActionAtomistic::lockRequests();
        Bias::lockRequests();
    }

    void KEnRefBias::unlockRequests() {
        ActionAtomistic::unlockRequests();
        Bias::unlockRequests();
    }

    // Resolve ambiguity for numerical derivatives (used for 'plumed driver --debug-forces')
    // Since KEnRef drives forces based on atomic coordinates, we use the Atomistic implementation.
    void KEnRefBias::calculateNumericalDerivatives(ActionWithValue *a) {
        ActionAtomistic::calculateNumericalDerivatives(a);
    }
} // namespace PLMD::kenref

//TODO setupConstantValues() may be used to setup the values in the first step?
//TODO link output file to action, and maybe to another file
//TODO use readAtomsFromPDB() from ActionAtomistic
//TODO use const `std::vector<AtomNumber>& PLMD::ActionAtomistic::getAbsoluteIndexes() const` and `AtomNumber PLMD::ActionAtomistic::getAbsoluteIndex(int i) const`
