#include "core/PlumedMain.h"
#include "core/ActionRegister.h"

namespace PLMD {

    void registerKEnRefModule(ActionRegister& reg);

    extern "C" void plumed_module_init(PlumedMain* plumed) {
//        ActionRegister& reg = plumed->getActionRegister();
//        registerKEnRefModule(reg);

        // Registration happens via PLUMED_REGISTER_ACTION macros
        // No need to do anything here
        (void) plumed; // Avoid unused parameter warning
    }
}
