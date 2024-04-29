# Dumpy: A Compact and Adaptive Index for Large Data Series Collections

Dumpy is an innovative data series index that utilizes an adaptive multi-ary data structure along with a node packing algorithm, specifically engineered to consolidate small leaf nodes. This approach significantly enhances the efficiency and search precision of SAX-based indexes. Additionally, there are two distinct variants of Dumpy: Fuzzy and Dumpy-Memory. These variants are tailored to deliver more precise approximate search outcomes, with Fuzzy optimized for on-disk datasets and Dumpy-Memory optimized for in-memory datasets.


## Build

We executed the program using Windows Subsystem for Linux (WSL).

To build the project,

1. create a "build" directory under the project

2. cd build

3. cmake ..

4. make

5. cd .. && ./bin/Dumpy

## Suggestions

1. All configuration details, including dataset information, are stored in the config.ini file. To complete your task, please carefully review and update the config.ini file with the necessary changes.

2. Before proceeding, it is essential to create a basic graph structure. This initial graph structure is denoted by index=0.

3. It is highly recommended to construct the SAX table before creating the entire index structure. This step brings benefits to all types of indexes.

4. Utilizing SIMD (Single Instruction, Multiple Data) techniques is crucial for optimizing the Dynamic Time Warping (DTW) distance calculations. However, if your machine does not support SIMD or you prefer not to use it, you can comment out the compile command labeled 'haswell' in the CMakeList.txt file. Additionally, you should comment out all code sections that utilize AVX2 commands. It's worth noting that all SIMD functions in our repository have an equivalent SISD (Single Instruction, Single Data) version available.

## Reproducibility
After compiling this project, you can produce the results by following procedures.
All the configurations are integrated into `config.ini`, including all the parameters and all the functions of this repo.

0. Build a skeleton graph. Set in `GraphConstruction.cpp` segment number (e.g., 16), and number of bits for fetching the neighbors of a SAX word, (e.g., 3, at most 4).
Set `index = 0` in `config.ini`, and run `./bin/Dumpy`.
After that, update the [other] section in `config.ini`.

1. Build the index. Set `index = 1` if Dumpy, and `index = 2` if Dumpy-Fuzzy. Set `ops = 0` for building index.
`th` is for the leaf node capacity, usually 10,000.
`segmentNum` is the number of segments, also set in `Const.h`, usually 16.
`bitsCardinality` is for the cardinality of bits of SAX words, usually 8.
`fbl_size` is the buffer size, in MB, depending on your machine.
`max_diff`, `fuzzy`, `delta` are for Dumpy-Fuzzy.
`small_perc` is to define a small leaf node for leaf node packing, 1.0 for an exhausted packing.
`max_mask_bit_percentage` is for ensuring the quality of the leaf pack, usually 0.8.
`f_low`, `f_high` are for adaptive split, usually 0.5 and 1.5.
`alpha` is the weighting factor in adaptive split

- [Dataset], the dataset you use should have a section in `config.ini`. Update the paths and the basic info inside.

After confirming these parameters, run `./bin/Dumpy`.

2. Run the queries.
- `index = 2 or 8` exact queries, ED and DTW, respectively.

After confirming these parameters, run `./bin/Dumpy`, then the results will be printed to your screen.

## Datasets

The dataset we used (Deep) are now in OneDrive.(Link to be added here)


# Reference
Zeyu Wang, Qitong Wang, Peng Wang, Themis Palpanas, and Wei Wang. 2023. Dumpy: A Compact and Adaptive Index for Large Data Series Collections. Proc. ACM Manag. Data 1, 1, Article 111 (May 2023), 27 pages. https://doi.org/10.1145/3588965


