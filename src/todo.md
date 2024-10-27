- [DONE] loop on all mRNA channels
- [DONE] use list generators
- [DONE] fix pyCharm python highlighting
- [DONE] check that it still works from the command line
- [DONE] make notebooks?
- [DONE] oneliner to open vsi picture + DIC layer (layer order, scale, color, blending)
- [DONE] fix Git situation
- [DONE] deduplicate images 240927
- [DONE] automatically set slice to in-focus slice (with value from spot detection)
- [DONE] add timings to output config
- [WONTDO] move (DIC, DAPI, all?) output files to their own channel folder (to accomodate Omnipose)
- [WONTDO] histograms in napari or matplotlib; automatic contrast adapted to mRNA channels and DIC
- [LATER] fix Git situation for real (git-lfs + Zenodo)
- [LATER] pytest tests
- [LATER] parallelize (dask)
- [LATER] fix envs mess
- [LATER] when using symlinks, check that they work after a round trip to OneDrive
- [DONE] analysis: per condition/growth rate, plot:
  - number of cells
  - cell area
  - cell eccentricity
  - number of spots/Txs
  - number of spots/Txs per cell
- [DONE] analysis: normalize all counts by number of cells
- [DONE] process just one picture
- [DONE] change size and color of spots 
- [NEXT] split all pictures in 4
- [NEXT] tune pipeline:
  1. spot detection
     1. [DONE] Remove slices that are further than +- n slices from focus
     2. [DONE] Intensity histogram for selected spots: ideally, we expect a bimodal distribution (individual spots and Tx's)
     3. [WONTDO] Remove large blotches: replace noisy rectangle with mean value of picture excluding all noisy rectangles / exclude them from spot assignment
     4. [DO] Check how many out-of-cells spots we detect
     5. [DO] Estimate background in the different channels, using patches that show background fluoresence (inside cells)
     6. Automatic threshold detection: 
        1. [DONE] Check if elbow is visible in number of selected spots versus threshold plot(Bigfish example notebook)
        2. [DO] Check whether threshold from different pictures from the same channel are consistent
        3. [DONE] Difference between LoG and remove_background_gausian?
        4. Automatic threshold looks too low for rpoD
           1. [DONE] Why does BigFish detect low-intensity spots?
        5. Automatic threshold looks too high for hipBA
  2. spot decomposition
     3. [DO] Tune distance between spots for Tx's
     4. [DO] Intensity histogram for selected, individual spots: they should be similar (within the same channel)
  3. spot counting/assignment 
     1. [DO] exclude cells on the border
     2. [DO] exclude spots from noisy regions
  4. DIC segmentation
  5. DAPI segmentation
  6. DAPI/DIC registration 
     1. [DO] check out skimage.registration
     2. [DO] split images in quadrants and align separately
  7. other
- [BIO] Are cells with 1 or 2 TA expressed different? Aspect, area, eccentricity...?
- [NEXT] optimise Omnipose segmentations
  - [DO] select "segmentable" images
  - [DO] if enough cells, don't waste time on optimizing segmentation. Just select enough adequate cells based on area/eccentricity. YESSSS!
    - "adequate" cells mean:
      - area and eccentricity are within "normal" range
      - there a nuclear mask within the cell mask (=discard cells that have no DAPI)
  - check ModelZoo for bacteria DIC model
  - check out Cellpose
- [LATER] per condition/growth rate, plot:
  - number of RNA per Tx
  - distance between spots/Txs and nucleoid
  - spots/Txs vs. cell size
  - spots/Txs vs. eccentricity
  - spots/Txs vs. growth rate
- [LATER] nucleoid segmentation
  - was this trip really necessary? Segmenting a 2D projection adds a lot of arbitrariness
  - can we not compute distance between putative transcription site and nucleoid in 3D?
  - NB: there is an issue with the DIC / DAPI alignment: the alignment is not good for the whole picture
- [MAYBE] deskew MG1655_LB images
- [MAYBE] use duplicate images to test analysis pipeline (do I get the same results with two duplicates?)
- [MAYBE] manually add a "valid" mask to exclude parts of the picture where DIC is screwed up (e.g. MG16555_MAN_OD_0.3_left_10)
- [NOTES] Omnnipose notes https://omnipose.readthedocs.io/training.html
  - cyto2_omni uses nuclei as chan2 -> do as well
    - make dapidic.tif composite to give to Omnipose
    - split images in 4 quadrants
    - do DIC/DAPI alignment per quadrant, or
      - correct chromatic aberration with digital filter
  - try bact_fluor model on fluorescence channels
  - points to remember if doing transfer learning:
    - diameter?
    - save_every
    - dataloader
    - DIC / DAPI alignment: "Tip :If using a transmissive modality like phase contrast or brightfield or DIC, use the same filter cube as your fluorescence channel. 
    - This usually removes any offset between the channels. Otherwise, be sure to do multimodal registration between the channels."
  - cell masks
    - dimensions: "You should aim to make training images of roughly size (512,512)"
    - object density: "As a general rule, you want to train on images with densely packed objects"
  - Split input image in 4:
    - 00-read-vsi.py: .vsi + DIC.tif -> {DIC, DAPI, rpoD, rnlAB, hipBA.tif}
      - take quadrants and store as .pt{1,2,3,4} (separate pictures in separate directories)
    - 01-configure.py: 
      - remove .vsi and DIC.tif reading in part
      - process as usual
