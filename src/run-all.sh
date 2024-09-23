python src/01-configure.py tests/input/config.json
python src/02-segment.py tests/output/config.json
python src/03-detect-spots.py tests/output/config.json
python src/04-decompose-spots.py tests/output/config.json
python src/05-assign-spots.py tests/output/config.json
# display results: 06-plot-results.ipynb
