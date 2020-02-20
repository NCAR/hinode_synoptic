#!/bin/bash

# Gather all the HOP 336 pointings from the observation plan comments
find $HOME/ops/cp/cmdpln | egrep '/20(1[7-9]|2[0-9])[0-9]{4}[_revise]*/re-point_[0-9]{12}\.txt$'| xargs -n1 egrep "HOP[ -]*(336|393)" | uniq | tee pointings.txt
