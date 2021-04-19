#!/bin/bash

if [[ ! -e ${PWD}/config.yml ]]
then
  cp /usr/local/bin/config.yml ${PWD}
fi

exec "${@}"
