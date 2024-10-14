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
tools.view.import_layers('tests/output-seg_ok/MG1655_LB_fixed2806_hybed1906_left_01', mode='all', viewer=viewer)
```
alors qu'on voudrait:
`python tools/view.py tests/output-seg_ok/MG1655_GLU_OD_0.3_left_02 all`
but:
`10/07/2024 03:38:28PM [INFO] No OpenGL_accelerate module loaded: No module named 'OpenGL_accelerate'`
DAPI/DIC alignment:
```
dx, dy = dic.translate; dx, dy = int(dx), int(dy); io.imsave('dicdapi.tif', np.stack([dic.data[dy:, :dx], dapi.data[:-dy, -dx:]]))
```
Most recently modified files on macos:
`find . -type f -exec stat -f "%m %N" "{}" \; | sort -nr | head`