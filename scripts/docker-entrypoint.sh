#!/bin/bash

if [[ ! -e ${PWD}/config.yml ]]
then
  cp /config.yml ${PWD}
fi

exec "${@}"
