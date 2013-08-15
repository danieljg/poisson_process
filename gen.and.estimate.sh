#!/bin/bash
echo generating data
time source dg.sh > dg.log
echo estimating error
time source estimate.sh > estimate.log
echo by daniel jim\'enez
