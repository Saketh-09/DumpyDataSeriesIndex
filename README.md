# Dumpy: A Compact and Adaptive Index for Large Data Series Collections

Dumpy is an innovative data series index that utilizes an adaptive multi-ary data structure along with a node packing algorithm, specifically engineered to consolidate small leaf nodes. This approach significantly enhances the efficiency and search precision of SAX-based indexes. Additionally, there are two distinct variants of Dumpy: Fuzzy and Dumpy-Memory. These variants are tailored to deliver more precise approximate search outcomes, with Fuzzy optimized for on-disk datasets and Dumpy-Memory optimized for in-memory datasets.


## Suggestions

1. All configuration details, including dataset information, are stored in the Const.py file.

2. Before proceeding, it is essential to create a basic graph structure. This initial graph structure is denoted by index=0.

3. It is highly recommended to construct the SAX table before creating the entire index structure. This step brings benefits to all types of indexes.

4. Utilizing SIMD (Single Instruction, Multiple Data) techniques is crucial for optimizing the Dynamic Time Warping (DTW) distance calculations. 

## Reproducibility
We can produce the results by following procedure.
All the configurations are integrated into `Const.py`, including all the parameters and all the functions of this repo.

Create a Virtual Environment by using:
python3 -m venv myenv 
source myenv/bin/activate

Install the dependencies:
pip install -r requirements.txt

Build the index by running the Python script with appropriate parameters specified in `Const.py`.

Run queries using the Python script to perform exact or approximate searches based on the specified index type (Dumpy or Dumpy-Fuzzy).


## Datasets

The dataset we used (Deep) are now in OneDrive.(https://cometmail-my.sharepoint.com/:f:/g/personal/sxa230044_utdallas_edu/EhlOXDuieqhPt1aEgq0ppK8BgS2NfBBY2ZuuE0JYYBrRrA?e=H7KCrb)
Download the data set & update the paths in Const.py

complete dataset link: Deep: http://sites.skoltech.ru/compvision/noimi

# Reference
Zeyu Wang, Qitong Wang, Peng Wang, Themis Palpanas, and Wei Wang. 2023. Dumpy: A Compact and Adaptive Index for Large Data Series Collections. Proc. ACM Manag. Data 1, 1, Article 111 (May 2023), 27 pages. https://doi.org/10.1145/3588965


