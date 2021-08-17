# Picocystis nanoSIMS analysis

## Directory
- `nanosims-processor.py`: Primary script for generating ROI-specific isotope ratios.
- `loop-nanosims-processor.sh`: Shell script that allows `nanosims-processor.py` to operate over many nanoSIMS datasets at once.
- `June_2021/`: Directory containing nanoSIMS data and ROIs.
- `shell_scripts`: Shell scripts to help with processing nanoSIMS ROIs.
- `roi_tbl.tsv`: Master dataframe for 15N analysis.
- `ROI_data.Rmd`: Processing and plotting nanoSIMS 15N data.
- `nanoSIMSiso.Rmd`: Calculating generation times.
- `checking_nh4_dark.Rmd`: Showing that activity in NH4_dark condition is not from erroneous light exposure.
- `asv_processing.Rmd`: 16S amplicon sequencing analysis.
- `asv_analyzed/`: 16S ASV dataframes.
- `mono_8m_96hr/`: Contains micrographs.

# Using the nanosims-processor script

### 0. Starting File Structure and naming scheme
For ease of use, it is best to nest the nanoSIMS data in the following way.
1. Main directory containing all nanoSIMS for a given project.
2. Sub-directories for different nanoSIMS runs
3. Sub-sub-directories for individual nanoSIMS images

### 1. Start with a `.im` file

```
~/
├── project_directory/ # nanoSIMS data directory
    ├── run1/
    |   ├── image1/
    |   |   └── image1.im
    |   ├── ...
    |   └── image2/
    ├── ...
    └── run2/
```

### 2. Draw your ROIs
This process requires that you have manually drawn ROIs for each nanoSIMS image you intend to analyze.
Most image manipulation software will work for this, but some of the more design-focused ones will use odd colors, add backgrounds, or have other issues. Adobe, GIMP, MS Paint, will all work just fine. The most important thing is consistency in ROI colors, and that the ROI colors must correspond to strict RGB codes.
- Green: `RGB: 0, 255, 0` `HEX: #00ff00`
- Red: `RGB: 255, 0, 0` `HEX: #ff0000`
- Blue: `RGB: 0, 0, 255` `HEX: #0000ff`

Append `_ROI.png` to your `image` name for ease of processing!

```
~/
├── project_directory/ # nanoSIMS data directory
    ├── run1/
    |   ├── image1/
    |   |   ├── image1_ROI.png
    |   |   └── image1.im
    |   ├── ...
    |   └── image2/
    ├── ...
    └── run2/
```
Great, now `nanosims-processor.py` has everything it needs to generate channel-specific images.

### 3. Generate channel-specific images

Now we can run `nanosims-processor.py` either individually for each image, or we can loop it over a series of images using `loop-nanosims-processor.sh`.

```
nanosims-processor.py --roi image1_ROI.png -i image1.im -o image1;
```
Returns:
```
~/
├── project_directory/ # nanoSIMS data directory
    ├── run1/
    |   ├── image1/
    |   |   ├── image1_12C.png
    |   |   ├── image1_13C.png
    |   |   ├── image1_14N12C.png
    |   |   ├── image1_15N12C.png
    |   |   ├── image1_31P.png
    |   |   ├── image1_32S.png
    |   |   ├── image1_34S.png
    |   |   ├── image1_SE.png
    |   |   └── image1_ROI.png
    |   ├── ...
    |   └── image2/
    ├── ...
    └── run2/
```
In the case of the Picocystis data, expanding out the `NH4_light_chain10A_1` image, we see the final result:
```
picocystis/
└── June_2021/
    ├── NH4_light/
    |   ├── GB21_L2_NH4_light_A_1_f0
    |   ├── GB21_L2_NH4_light_A_1_f1
    |   ├── GB21_L2_NH4_light_B_1_f0
    |   ├── GB21_L2_NH4_light_B_2_f0
    |   ├── GB21_L2_NH4_light_B_2_f0
    |   ├── ...
    |   └── GB21_L2_NH4_light_chain10A_1_f1
    |       ├── GB21_L2_NH4_light_chain10A_1_f1_12C.png
    |       ├── GB21_L2_NH4_light_chain10A_1_f1_13C.png
    |       ├── GB21_L2_NH4_light_chain10A_1_f1_14N 12C.png
    |       ├── GB21_L2_NH4_light_chain10A_1_f1_15N 12C.png
    |       ├── GB21_L2_NH4_light_chain10A_1_f1_31P.png
    |       ├── GB21_L2_NH4_light_chain10A_1_f1_32S.png
    |       ├── GB21_L2_NH4_light_chain10A_1_f1_34S.png
    |       ├── GB21_L2_NH4_light_chain10A_1_f1_SE.png
    |       ├── GB21_L2_NH4_light_chain10A_2_f1_ROI.png # The ROI you drew! Note that it shares its filename with the image directory.
    |       ├── GB21_L2_NH4_light_chain10A_2_f1_ratio15N_12C-x-14N_12C.tsv # OUTPUT ROI data summary
    |       └── GB21_L2_NH4_light_chain10A_2_f1_ratio15N_12C-x-14N_12C.png # OUTPUT ROI visualization
    ├── NH4_dark/
    ├── NO3_dark/
    ├── NO3_light/
    ├── ...
    └── gly_dark
```
### 4. Collecting .tsv output
Now we have lots of `.tsv` files, each corresponding to a single image, containing ROI-specific data. We can work with these individually or, if we prefer, iteratively collect all of the .tsv files into a folder called `roi_data` using a shell one-liner like the one found in `extract_tsv.sh`. In short, this code searches two directories down `**/**/` for files ending with `.tsv`.
```
n | cp -i **/**/*.tsv roi_data/
```

# Handoff to R
Once all of our ROI data is neatly contained within `roi_data/`, importing into Rstudio (or other) is straightforward.

In this repo, `ROI_data.Rmd` walks us through this. In brief:

**Import ROI Data**
```
temp_list = list.files(path = "June_2021/roi_data", pattern="*.tsv")

# Create a list of columns corresponding to those in the .tsv files
tbl_colnames <- c("N_source", "light", "inc_time_hr", "ROI", "12C", "13C", "14N_12C", "15N_12C", "31P", "32S", "34S", "SE", "Ratio_15N_12Cx14N_12C" )
# Create empty tibble to be filled
roi_tbl <- tibble()
roi_tbl <- roi_tbl[-1,] # Remove the garbage row of the tibble

# Read in each .tsv, iteratively append each to our `roi_tbl`
for (tsv in temp_list) {
  tmp_df <- read_tsv(paste0("June_2021/roi_data/", tsv))
  tmp_df <- tmp_df %>% mutate(filename = tsv)
  roi_tbl <- bind_rows(roi_tbl, tmp_df)
}
```
**Inspect**
Inspecting the data, we see it is indeed quite tidy and ready for analysis!
```
> roi_tbl
# A tibble: 438 × 16
   Group   ROI   `12C` `13C` `14N_12C` `15N_12C` `31P`  `32S` `34S`       SE Ratio_15N_12Cx14… filename  
   <chr> <dbl>   <dbl> <dbl>     <dbl>     <dbl> <dbl>  <dbl> <dbl>    <dbl>             <dbl> <chr>     
 1 red       0 7692870 81630    769335      2901  1309 133310  5906 11880832           0.00376 1_CN_ligh…
 2 red       1 5421421 57308    589241      2198  1380 159254  6941  9780419           0.00372 1_CN_ligh…
 3 blue      0  553515  5935     88941       359   150  10528   443  1048372           0.00402 1_CN_ligh…
 4 blue      1  985215 10570     94778       355   263  17104   777  1516020           0.00373 1_CN_ligh…
 5 blue      2 1851684 19650   1753653     12111 20884 256457 11370  2832207           0.00686 1_CN_ligh…
 6 red       0 4774159 49723   1798080     13140  3022 176625  7583 12558414           0.00725 1_CN_ligh…
 7 red       0 7393259 78017   2535672     25564 21557 584673 26096 11187109           0.00998 1_CN_ligh…
 8 blue      0 1380422 14699    150733       571   416  16807   759  2333582           0.00377 1_CN_ligh…
 9 blue      1 1167079 12612     44057       206   188  13096   579  1613073           0.00465 1_CN_ligh…
10 blue      2  701304  7526     41563       184   161  16677   756  1137139           0.00441 1_CN_ligh…
# … with 428 more rows, and 4 more variables: sample <chr>, layer <chr>, atom % N15 <dbl>, notes <chr>
```
