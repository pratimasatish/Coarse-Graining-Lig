#!/bin/bash --login

awk '(/^3/) && $3 >= 5.7 {print NR-2, $0}' first-frame.xyz > surf-Cd.txt
# awk '(/^2/) && $3 <= 7.7 && $3>0 {print NR-2, $0}' first-frame.xyz > surf-lig.txt
