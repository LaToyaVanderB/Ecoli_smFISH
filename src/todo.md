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
- [NEXT] figure out if/where pipeline fails: check:
  1. spot detection
     1. Filter out points that are totally out of focus (e.g. GLU01 hipBA blotches)
     2. Blotches screw up processing? (GLU01 hipBA)
  2. spot decomposition
  3. spot counting/assignment 
  4. DIC segmentation
  5. DAPI segmentation
  6. DAPI/DIC registration 
  7. other
- [NEXT] optimise Omnipose segmentations
  - (Philipp) If you have enough cells, don't waste time on optimizing segmentation. Just select enough adequate cells besed on area/eccentricity. YESSSS!
  - check ModelZoo for bacteria DIC model
  - Omnnipose notes https://omnipose.readthedocs.io/training.html
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