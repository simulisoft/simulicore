#!/bin/bash

ROOT_DIR=`dirname $0`

sed -i 's/#cmpi_tutorial_comm\>//g' ${ROOT_DIR}/html/*.js

exit 0
