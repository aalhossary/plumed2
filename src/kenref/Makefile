USE=core tools vesselbase multicolvar bias

# generic makefile
#This line MUST be at the beginning. Before (NOT AFTER) the variables.
include ../maketools/make.module

# make KEnRef double
CPPFLAGS += -DKENREF_DOUBLE

# Debug output
$(info ===== Debug Info =====)
# Debug: Print variables
$(info >>>>>>>>>>>>> KENREF_INCLUDE is: $(KENREF_INCLUDE))
$(info >>>>>>>>>>>>> KENREF_LIBS is: $(KENREF_LIBS))
$(info >>>>>>>>>>>>> HAVE_KENREF: $(HAVE_KENREF))
$(info >>>>>>>>>>>>> Before setting: CPPFLAGS=$(CPPFLAGS))
$(info >>>>>>>>>>>>> Before setting: LIBS=$(LIBS))
$(info =====================================)
# Set flags BEFORE including make.module
ifneq ($(KENREF_INCLUDE),)
  CPPFLAGS+=$(KENREF_INCLUDE)
else
  $(warning KENREF_INCLUDE is empty! We will hardcode it)
  # Hardcode for testing
  CPPFLAGS+=-I/smithlab/opt/kenref-dev/2025/2025.3/debug/AVX_512/include/core -I/smithlab/opt/eigen/3.4.0/include/eigen3
endif

ifneq ($(KENREF_LIBS),)
  LIBS+=$(KENREF_LIBS)
else
  $(warning KENREF_LIBS is empty! We will hardcode it)
  # Hardcode for testing
  LIBS+=-L/smithlab/opt/kenref-dev/2025/2025.3/debug/AVX_512/lib -lkenref_core -L/smithlab/home/opt/rhel08/openmpi/5.0.7_clang20.1.1_cuda12.4.1/lib/libmpi.so
endif

$(info >>>>>>>>>>>>> After setting: CPPFLAGS=$(CPPFLAGS))
$(info >>>>>>>>>>>>> After setting: LIBS=$(LIBS))
$(info =====================================)
