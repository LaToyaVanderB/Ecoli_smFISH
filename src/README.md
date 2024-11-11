Workflow:

1. Config step
   `python src/01-configure.py tests/input/config.json`
2. Segmentation 
For now, we do this from a terminal, and we need to switch to an Omnipose-capable environment.
   `python src/02-segment.py tests/output/config.json`
3. Segmentation postprocessing
   Run "Cell masks preprocesing"  cells from `segment.ipynb` 
4. Crop
Open segmented image in Napari, record top-left and bottom-right corners in config.json.
`"crop": [ymin, xmin, ymax, xmax]`
If no cropping, there should be no `"crop"` entry (i.e. `"crop": false` will fail)
4. Spot detection
   `python src/03-detect-spots.py tests/output/config.json`
5. Spot decomposition
   `python src/04-decompose-spots.py tests/output/config.json`
6. Spot assignment
   `python src/05-assign-spots.py tests/output/config.json`
7. Result analysis
   `06-plot-results.ipynb`

Other useful commands:
1. View json file:
`python -m json.tool tests/output/MG1655_LB_01/img.json`
2. View image: in Napari viewer console, do:
```
import sys; sys.path.insert(0, ''); import tools.view; from importlib import reload
tools.view.import_layers('/Volumes/Flodrive/Florence/smFISH/analysis/20240927-LB/MG1655_LB_fixed2806_hybed1906_left_01', mode='cells', viewer=viewer)
tools.view.import_layers('tests/exp16/output/MG1655_GLU_OD_0.3_left_02', mode='cells', viewer=viewer)
tools.view.import_layers('tests/output-seg_ok/MG1655_GLU_OD_0.3_left_02', mode='all', viewer=viewer)
tools.view.import_layers('/Volumes/Flodrive/Florence/smFISH/Ecoli_smFISH/tests/output-seg_ok/MG1655_GLU_OD_0.3_left_02', mode='spots', viewer=viewer)
```
3. Most recently modified files on macos:
`find . -type f -exec stat -f "%m %N" "{}" \; | sort -nr | head`
4. 'max project' spots:
`np.delete(spots, 0, 1)`
or use 'out_of_slice_display': True in napari viewer add_points command
