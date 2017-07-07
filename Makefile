
# 1. Adapt PATHS, COMPILER, and FLAGS if necessary
# 2a. Use 'make' to compile everything (the libShc, the extension libShx for openCV image loading/saving, and examples)
# 2b. Use 'make lib_shc' to compile the lib only (no dependencies on OpenCV)

# PATHS
pathKissFFT = ./kissfft

# COMPILER
Warnings = -Wall -Wextra -pedantic
Compile  = -march=native -std=c++11 -DEIGEN_NO_DEBUG -O3 

# FLAGS
ShcFlags     = -I./include
EigenFlags   = `pkg-config eigen3 --cflags`
OpenCVFlags  = `pkg-config opencv --cflags`
KissFFTFlags = -I$(pathKissFFT)
Flags        = $(Compile) $(Warnings) $(ShcFlags) $(EigenFlags) $(OpenCVFlags) $(KissFFTFlags)

# LIBS
OpenCVLibs   = `pkg-config opencv --libs`
Libs         = $(OpenCVLibs)



.SUFFIXES:

all: lib_shc lib_shx examples

lib_shc: build/libShc.a
lib_shx: build/libShx.a
examples: examples/compass examples/tangentDistance examples/compass examples/compass_masked examples/hemispherical_continuation examples/insect_eye examples/ocamcalib examples/pointclouds examples/simple examples/seqSlam 

build/libShc.a: build/Shc_cg.o build/Shc_compass.o build/Shc_conversions.o build/Shc_feature.o build/Shc_fft.o build/Shc_file.o build/Shc_filter.o build/Shc_indices.o build/Shc_init.o build/Shc_misc.o build/Shc_print.o build/Shc_rotate.o build/Shc_surface.o build/Shc_tangentDistance.o build/Shc_transform.o build/Shc_interface_conversions.o build/Shc_interface_func.o build/Shc_interface_init.o build/Shc_interface_misc.o build/Shc_interface_options.o build/Shc_interface_print.o build/Shc_interface_transform.o build/kiss_fft.o build/kiss_fftr.o

build/libShx.a: build/Shx_wrapper.o build/Shx_intern.o

examples/%: examples/%.o
	g++ -L./build $@.o -lShx -lShc -lrt $(Libs) -o $@

%.a:
	ar rvs $@ $+

build/%.o: src/%.cpp
	g++ $(Warnings) $(Flags) -c $< -o $@

build/%.o: $(pathKissFFT)/%.c
	g++ $(Warnings) $(Flags) -c $< -o $@

examples/%.o: src_examples/%.cpp
	g++ $(Warnings) $(Flags) -c $< -o $@

clean:
	find build/ -maxdepth 1 -type f -delete
	find examples/ -maxdepth 1 -type f -delete