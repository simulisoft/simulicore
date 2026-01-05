#!/bin/bash

ROOT_DIR=`dirname $0`

sed -i 's/#culite_tutorial_datatypes\>//g' ${ROOT_DIR}/html/*.js
sed -i 's/#culite_tutorial_vectors\>//g' ${ROOT_DIR}/html/*.js
sed -i 's/#culite_tutorial_vectors_dense\>//g' ${ROOT_DIR}/html/*.js
sed -i 's/#culite_tutorial_matrices\>//g' ${ROOT_DIR}/html/*.js
sed -i 's/#culite_tutorial_matrices_dense\>//g' ${ROOT_DIR}/html/*.js
sed -i 's/#culite_tutorial_matrices_sparse\>//g' ${ROOT_DIR}/html/*.js
sed -i 's/#culite_tutorial_linsol\>//g' ${ROOT_DIR}/html/*.js
sed -i 's/#culite_tutorial_linsol_lapack\>//g' ${ROOT_DIR}/html/*.js
sed -i 's/#culite_tutorial_linsol_cudss\>//g' ${ROOT_DIR}/html/*.js

exit 0
