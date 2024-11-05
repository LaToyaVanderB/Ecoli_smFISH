1. Config step
   `python src/01-configure.py tests/input/config.json`
2. Segmentation (for now, need to switch to omnipose environment)
   `python src/02-segment.py tests/output/config.json`
3. Spot detection
   `python src/03-detect-spots.py tests/output/config.json`
4. Spot decomposition
   `python src/04-decompose-spots.py tests/output/config.json`
5. Spot assignment
   `python src/05-assign-spots.py tests/output/config.json`
6. Result analysis
   `06-plot-results.ipynb`
7. Bonus
`python -m json.tool tests/output/MG1655_LB_01/img.json`
```
import sys; sys.path.insert(0, ''); import tools.view
tools.view.import_layers('tests/exp16/output/MG1655_GLU_OD_0.3_left_02', mode='cells', viewer=viewer)
tools.view.import_layers('tests/output-seg_ok/MG1655_GLU_OD_0.3_left_02', mode='all', viewer=viewer)
tools.view.import_layers('/Volumes/Flodrive/Florence/smFISH/Ecoli_smFISH/tests/output-seg_ok/MG1655_GLU_OD_0.3_left_02', mode='spots', viewer=viewer)
```
alors qu'on voudrait:
`python tools/view.py tests/output-seg_ok/MG1655_GLU_OD_0.3_left_02 all`
but:
`10/07/2024 03:38:28PM [INFO] No OpenGL_accelerate module loaded: No module named 'OpenGL_accelerate'`
DAPI/DIC alignment and combined image for Omnipose segmentation:
```
dy, dx = dic.translate; dy, dx = int(dy), int(dx); io.imsave('dicdapi.tif', np.stack([dic.data[dy:, :dx], dapi.data[:-dy, -dx:], np.zeros[], np.zeros[]]))
```
Most recently modified files on macos:
`find . -type f -exec stat -f "%m %N" "{}" \; | sort -nr | head`

Max project spots:
`np.delete(spots, 0, 1)`
or use 'out_of_slice_display': True in napari viewer add_points command
