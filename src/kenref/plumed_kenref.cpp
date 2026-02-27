////#ifdef __PLUMED_HAS_KENREF
//#include "bias/Bias.h"
//#include "tools/Matrix.h" // For coordinate handling
//#include "tools/OpenMP.h"
////#include "tools/Parser.h"   // For reading PLUMED input
//
//#include "plumed_kenref.h"
//
//// Include your KEnRef headers
//#include "core/KEnRef.h"
//#include "core/Table.h"
//// ... other KEnRef headers
//
//namespace PLMD {
//namespace bias {
//
//// Use the appropriate precision (double is standard for PLUMED/MD)
//using Real = double;
////using KEnRef_Real = double;
//
//class KEnRefSigmaBias : public Bias {
//private:
//    // Parameters read from PLUMED.dat
//    KEnRef_Real_t k_force;
//    KEnRef_Real_t n_power;
//    KEnRef_Real_t proton_mhz;
//
//    // Static KEnRef data structures
//    std::vector<SpecDenData<KEnRef_Real_t>> spec_den_data_list;
//    // ... other static data (rates, atomNames_2_atomIds, etc.)
//
//public:
//    // Constructor (Initialization)
//    KEnRefSigmaBias(const ActionOptions& ao) : Bias(ao) {
//        // --- 1. Read parameters from PLUMED.dat ---
//        pld::read("K", k_force);
//        pld::read("N", n_power);
//        pld::read("PROTON_MHZ", proton_mhz);
//
//        // --- 2. Load static data (Requires custom code here) ---
//        // Example: load files to populate spec_den_data_list
//        // This is where you would use pld::file or similar tools to read external files
//        // and populate the 'spec_den_data_list', 'rates', and 'atomNames_2_atomIds'.
//        // ... LoadSpecDenData(spec_den_data_list, ...);
//
//        // This action needs coordinates for all atoms
//        addRequiredAtomicPositions();
//    }
//
//    // Main calculation function, called every step
//    void calculate() override {
//        // --- 3. Get Coordinates from PLUMED ---
//        // PLUMED provides coordinates for a single model (the current MD step)
//        // You may need to adapt this section for ensemble simulations.
//
//        // Get positions as a raw Eigen-compatible matrix (must be done carefully)
//        const auto& p = getPositions(); // PLUMED's position wrapper
//
//        // Create the coordinate structure expected by KEnRef.
//        // Assuming a single structure simulation for simplicity:
//        CoordsMatrixType<KEnRef_Real_t> current_coords(p.rows(), 3);
//        for(int i = 0; i < p.rows(); ++i) {
//            current_coords(i, 0) = p[i][0];
//            current_coords(i, 1) = p[i][1];
//            current_coords(i, 2) = p[i][2];
//        }
//        std::vector<CoordsMatrixType<KEnRef_Real_t>> coord_array_for_kenref = {current_coords};
//
//        // --- 4. Call the core KEnRef function ---
//        bool gradient = true; // Always true for a Bias
//        int numOmpThreads = 0; // Use PLUMED's threading model
//
//        auto [energy, d_energy_d_coord_array_opt] =
//            KEnRef<KEnRef_Real_t>::coord_array_to_sigma_energy(
//                coord_array_for_kenref, rates, spec_den_data_list,
//                proton_mhz, k_force, n_power, atomNames_2_atomIds, gradient, numOmpThreads
//            );
//
//        // --- 5. Write Energy back to PLUMED ---
//        setEnergy(energy);
//
//        // --- 6. Write Forces back to PLUMED ---
//        if (d_energy_d_coord_array_opt.has_value()) {
//            const auto& dE_dx_array = d_energy_d_coord_array_opt.value();
//
//            // Assuming single model/structure in the array (index 0)
//            const auto& dE_dx_matrix = dE_dx_array.at(0);
//
//            for (int i = 0; i < dE_dx_matrix.rows(); ++i) {
//                // Force = -Gradient
//                Real force_x = -dE_dx_matrix(i, 0);
//                Real force_y = -dE_dx_matrix(i, 1);
//                Real force_z = -dE_dx_matrix(i, 2);
//
//                // Add the calculated force to the total force on atom 'i'
//                addForce(i, {force_x, force_y, force_z});
//            }
//            //Adding something to break
//            for (int i = 0; i < dE_dx_matrix.rows(); ++i) {
//            } //TODO This line should NOT break
//        }
//    }
//};
//
//// Required PLUMED registration boilerplate
//PLUMED_REGISTER_ACTION(KEnRefSigmaBias, "KENREF_SIGMA_BIAS")
//
//}
//}
//
////#endif