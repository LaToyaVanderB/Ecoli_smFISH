1. Config step
   `~/PycharmProjects/Ecoli_smFISH % python src/01-configure.py tests/input/config.json`
2. Segmentation (for now, need to switch to omnipose environment)
   `~/PycharmProjects/Ecoli_smFISH % python src/02-segment.py tests/output/config.json`
3. Spot detection
   `~/PycharmProjects/Ecoli_smFISH % python src/03-detect-spots.py tests/output/config.json`
4. Spot decomposition
   `~/PycharmProjects/Ecoli_smFISH % python src/04-decompose-spots.py tests/output/config.json`
5. Spot assignment
   `~/PycharmProjects/Ecoli_smFISH % python src/05-assign-spots.py tests/output/config.json`
6. Result analysis
   `06-plot-results.ipynb`

