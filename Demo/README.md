# Fitting Code Demo

###################################

## Running the demo

This folder contains a demo_data.csv file which is a subset of our full dataset, containing data from just 5 TPCs. This also contains a /Results/ folder where demo outputs will be written to. The contents of the Results folder can be deleted before running the demo, to see them written out again.

To fit the 4-parameter Sharpe-Schoolfield model to these data, simply run the demo_script.sh wrapper script from the command line.

    > bash demo_script.sh
    
## Expected outputs

The command line will feed back parameter estimates from each fit and the folder locations of the outputs. Figures from each of the 5 fits should be found in the /Results/ directory, along with summary.csv detailing the fitting results.
