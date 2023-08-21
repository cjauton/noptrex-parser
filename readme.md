# NOPTREX Parser

The `noptrex.cpp` program is a physics data analysis tool that performs various analysis tasks on raw binary data files and produces results for further investigation.

## Description

The `noptrex.cpp` program includes and uses several libraries such as `TCanvas.h`, `TLeaf.h`, `TH2F.h`, `TSystem.h`, `TTreeFormula.h`, `TStopwatch.h`, `TFile.h`, and `TTree.h`. It also defines several constants and macros to handle file paths, configuration settings, and data processing parameters.

The program reads data from input files, applies decimation, performs windowing, and calculates average signals for specified ranges. It then matches the average signals to predefined thresholds to determine specific states and sequences. The results are stored in a ROOT TTree object.

## Usage

To use the `noptrex.cpp` program, follow these steps:

1. Compile the program using a ROOT-compatible compiler.
2. Ensure that the required data files are present in the appropriate directories (`CONF`, `RUN`, and `SPIN`).
3. Modify the `CONF` file to set the necessary environment variables (configurations) for your analysis. Refer to the provided example values in the `noptrex.cpp` program and adjust them according to your specific requirements.
4. Run the compiled executable with the desired parameters. For example:

    ```bash
    ./noptrex -run 123 -verbose 1 -opt "some_option"
    ```

    Replace 123 with the desired run number. Adjust the verbose parameter (0 or 1) to control the verbosity of the program's output. Modify "some_option" with any additional options you want to pass to the program.

5. The program will process the data files based on the provided configurations and parameters.
6. The results will be stored in a ROOT TTree object.

Or run from the ROOT interactive terminal by running `root` then load the `noptrex.cpp` macro:

```c++
root .L noptrex.cpp
noptrex(<run_number>)
```

Use the generated results for further analysis or visualization using ROOT libraries or tools.

## Requirements

- ROOT (CERN) library

## Contributing

Contributions to this project are welcome. If you have any suggestions, improvements, or bug fixes, please create an issue or submit a pull request.

## License

This program is provided under an open-source license. Please see the LICENSE file for more information.
