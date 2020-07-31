#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux
CND_DLIB_EXT=so
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/d_mcarlo.o \
	${OBJECTDIR}/fft_convolve_3d.o \
	${OBJECTDIR}/input_module.o \
	${OBJECTDIR}/link_gaze_grid.o \
	${OBJECTDIR}/lm_exact.o \
	${OBJECTDIR}/main.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=-fdefault-integer-8 -m64 -I/opt/intel/mkl/include -I/opt/intel/mkl/include/fftw -g -fno-range-check

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/nss_core

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/nss_core: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	gfortran -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/nss_core ${OBJECTFILES} ${LDLIBSOPTIONS} -Wl,--start-group /opt/intel/mkl/lib/intel64/libmkl_gf_ilp64.a /opt/intel/mkl/lib/intel64/libmkl_sequential.a /opt/intel/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

${OBJECTDIR}/d_mcarlo.o: d_mcarlo.f08
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/d_mcarlo.o d_mcarlo.f08

${OBJECTDIR}/fft_convolve_3d.o: fft_convolve_3d.f08
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/fft_convolve_3d.o fft_convolve_3d.f08

${OBJECTDIR}/input_module.o: input_module.f08
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/input_module.o input_module.f08

${OBJECTDIR}/link_gaze_grid.o: link_gaze_grid.f08
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/link_gaze_grid.o link_gaze_grid.f08

${OBJECTDIR}/lm_exact.o: lm_exact.f08
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/lm_exact.o lm_exact.f08

${OBJECTDIR}/main.o: main.f08
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/main.o main.f08

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} *.mod

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
