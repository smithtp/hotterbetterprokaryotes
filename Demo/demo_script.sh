#!/bin/bash
# Run the demo script to fit Schoolfield model to TPC data

echo "Fitting from database..."

python3.4 ../Code/demo.py

echo "Finished fitting! Figures and parameter outputs written to /Results/"
