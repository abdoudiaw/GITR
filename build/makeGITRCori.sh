source ../env.cori.sh

cmake -DTHRUST_INCLUDE_DIR=/global/homes/t/tyounkin/code/thrust/ \
-DLIBCONFIGPP_INCLUDE_DIR=/global/homes/t/tyounkin/code/libconfigBuild/gnu/include \
-DLIBCONFIGPP_LIBRARY=/global/homes/t/tyounkin/code/libconfigBuild/gnu/lib/libconfig++.so \
-DCMAKE_C_COMPILER=gcc \
-DCMAKE_CXX_COMPILER=g++ \
-DNETCDF_LIBRARY=$NETCDF_DIR/lib/libnetcdf.so \
-DNETCDF_INCLUDE_DIR=$NETCDF_DIR/include \
-DNETCDF_CXX_INCLUDE_DIR=/opt/cray/pe/netcdf/4.6.1.3/GNU/8.2/include \
-DNETCDF_CXX_LIBRARY=$NETCDF_DIR/lib/libnetcdf_c++4_gnu_82.so \
-DUSE_CUDA=0 \
    -DUSEMPI=0 \
    -DUSE_MPI=1 \
    -DUSE_OPENMP=0 \
    -DUSE_BOOST=1 \
    -DUSEIONIZATION=1 \
    -DUSERECOMBINATION=1 \
    -DUSEPERPDIFFUSION=1 \
    -DUSEPARDIFFUSION=1 \
    -DUSECOULOMBCOLLISIONS=1 \
    -DUSEFRICTION=1 \
    -DUSEANGLESCATTERING=1 \
    -DUSEHEATING=1 \
    -DUSETHERMALFORCE=0 \
    -DUSESURFACEMODEL=1 \
    -DUSESHEATHEFIELD=1 \
    -DBIASED_SURFACE=1 \
    -DUSEPRESHEATHEFIELD=1 \
    -DBFIELD_INTERP=0 \
    -DLC_INTERP=0 \
    -DGENERATE_LC=0 \
    -DEFIELD_INTERP=0 \
    -DPRESHEATH_INTERP=0 \
    -DDENSITY_INTERP=2 \
    -DTEMP_INTERP=2 \
    -DFLOWV_INTERP=0 \
    -DGRADT_INTERP=0 \
    -DODEINT=0 \
    -DFIXEDSEEDS=0 \
    -DPARTICLESEEDS=1 \
    -DGEOM_TRACE=0 \
    -DGEOM_HASH=3 \
    -DGEOM_HASH_SHEATH=3 \
    -DPARTICLE_TRACKS=0 \
    -DPARTICLE_SOURCE_SPACE=0 \
    -DPARTICLE_SOURCE_ENERGY=0 \
    -DPARTICLE_SOURCE_ANGLE=0 \
    -DPARTICLE_SOURCE_FILE=1 \
    -DSPECTROSCOPY=3 \
    -DUSE3DTETGEOM=1 \
    -DUSECYLSYMM=1 \
    -DFLUX_EA=1 \
    -DUSEFIELDALIGNEDVALUES=0 \
    -DUSE_SORT=0 \
    -DFORCE_EVAL=0 \
    -DCHECK_COMPATIBILITY=1 \
..
#-DMPI_C_LIBRARIES=/opt/cray/pe/mpt/7.7.6/gni/mpich-gnu/8.2/lib/libmpich.so \
#-DMPI_C_LIB_NAMES=/opt/cray/pe/mpt/7.7.6/gni/mpich-gnu/8.2/lib/libmpich.so \
#-DMPI_C_INCLUDE_PATH=/opt/cray/pe/mpt/7.7.6/gni/mpich-gnu/8.2/include \
#-DMPI_CXX_LIBRARIES=/opt/cray/pe/mpt/7.7.6/gni/mpich-gnu/8.2/lib/libmpichcxx_gnu_82.so \
#-DMPI_CXX_LIB_NAMES=/opt/cray/pe/mpt/7.7.6/gni/mpich-gnu/8.2/lib/libmpichcxx_gnu_82.so \
#-DMPI_CXX_INCLUDE_PATH=/opt/cray/pe/mpt/7.7.6/gni/mpich-gnu/8.2/include \
#-DMPI_C_COMPILER=mpich \
#-DMPI_CXX_COMPILER=mpichcxx \
