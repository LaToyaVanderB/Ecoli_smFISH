- decompose_spots
- assign_spots
- jsonpickle handler for pandas
- increase IDE memory
- 
- view
- loop on all images in an experiment
- plot results: make results Excel sheet
- make test cases
      - without crop
      - with crop
- napari plugin 
- got into BigFish detect_spots
  - filtering?
  - get spots with threshold and actual BigFish intensities

- env:
1. napari
   1. check that console works
   2. open vsi files: conda install aicsimageio
2. napari plugin development 
   1. conda install napari-flofish
3. project
    1. omnipose: install from git clone on MacOS
        1. cd /Users/adele/Projects/omnipose && pip install -e .
        2. cd /Users/adele/Projects/cellpose-omni && pip install -e .
        3. GUI: omnipose
    2. pip install jsonpickle bioio bioio-ome-tiff bioio-ome-zarr bioio-bioformats
    3. pip install big-fish
   4. conda install scyjava
   5. export JAVAHOME=/Users/adele/miniconda3/envs/omnipose/lib/jvm
   5. set $JAVA_HOME for bioio-bioformats
    4. edit peakdetect.py to: from scipy.fft import ifft
    5. check that tests pass
4. 


Omnipose install: try
- install from Git clone with pip install -e 
  - from Git clone repo: pip install -e .
  - from project dir: pip install -e <path_to_git_clone_repo>
- install from Git repo
  - pip install git+https://github.com/kevinjohncutler/omnipose.git
  - pip install git+https://github.com/kevinjohncutler/cellpose-omni.git
- install pip modules did not work so far
  - you only need to install omnipose, but you need to uninstall omnipose and cellpose_omni: pip uninstall omnipose cellpose_omni 