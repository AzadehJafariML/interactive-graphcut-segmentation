# Graph-Cut Segmentation: Usage Instructions

## Requirements

Before running the code, make sure the **NIfTI toolbox** is installed and added to the MATLAB path. If you do not already have it, download the toolbox and add it to your MATLAB library/path before starting the segmentation pipeline.

## Running the Pipeline

1. Open and run the main script:

```matlab
MainFlow
```

All major commands in the script include inline comments explaining their purpose. Please review these comments and the corresponding references for additional methodological details.

2. In `MainFlow`, manually enter the input file name in **line 13**.

3. Enter the desired slice number in **prompt line 24**.

## Preprocessing

The script `skull_stripping_preprocess` removes non-brain tissues, including the skull and surrounding structures.

In this script, the radius parameter of the `imopen` command can be adjusted to improve skull-stripping performance for each slice. The default radius value is **55**.

```matlab
imopen(..., strel('disk', 55))
```

Depending on the image quality and slice characteristics, this value may need to be optimized manually for each dataset.

## Seed Selection

After preprocessing, an image will be displayed on the screen for manual seed selection.

Seeds must be selected for:

1. Object/foreground region
2. Background region

To draw seeds:

1. Hold the **right mouse button**.
2. Move the mouse over the desired region.
3. First select the object/foreground seeds.
4. Press **Enter** after finishing the object seed selection.
5. Then select the background seeds in the same way.
6. Press **Enter** again when finished.

The selected seed points are saved into matrices and used to construct the probability distribution functions (PDFs) for the regional term in the graph-cut segmentation model (cost function).

## Parameter Tuning

In the script `AlphaGraph`, the values of **lambda** and **gamma** must be determined experimentally to obtain the best segmentation result.

These parameters remain constant for the whole data.

## Runtime

The segmentation process may take some time depending on the image size and parameter settings. In some cases, processing may take approximately **15–20 minutes**.

## Visualizing the Results

After the segmentation process is complete, run:

```matlab
Segment_Results
```

This script displays the final segmentation output.

## Leakage Correction

If the final segmentation contains leakage or unwanted segmented regions, run:

```matlab
Remove_leakage
```

This script can be used to remove leakage artifacts from the final result. **Two examples** of leakage are porivded for better understanding the leakage issue

## Evaluation

To perform quantitative evaluation and statistical analysis, use these parameters to compare with Gold_standard

* Hausdorff distance
* Absolute volume difference (AVD)
* Statistical testing using t-test

```matlab
Lvalidationmode
```
