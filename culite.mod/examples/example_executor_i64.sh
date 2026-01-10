#!/bin/bash
set -e
ROOT_DIR=`dirname $0`
if [ ! -d "${ROOT_DIR}/i64/out" ]; then
  mkdir -p ${ROOT_DIR}/i64/out
fi

for file in ${ROOT_DIR}/i64/bin/*.sh
do
	echo "Running ${file%.sh}"
  ofile="$(basename  "${file%.sh}").out"
  "$file" > "${ROOT_DIR}/i64/out/$ofile"
done

exit 0
