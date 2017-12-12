#!/usr/bin/env bash

cat goana.R \
    | sed 's/Hs/Mm/g' \
    > goana.mouse.R
