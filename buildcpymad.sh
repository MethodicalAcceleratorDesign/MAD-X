mkdir build
cd build
cmake \
    -DMADX_X11=OFF \
    -DMADX_STATIC=ON \
    -DMADX_ONLINE=OFF \
    -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_INSTALL_PREFIX=$MADXDIR \
    -DCMAKE_C_FLAGS="-fvisibility=hidden" \
    -DCMAKE_BUILD_TYPE=Release \
    -DMADX_INSTALL_DOC=OFF \
    ..
make install
echo $MADXDIR