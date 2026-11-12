# Converting `tunedxgb.rds` to `tunedxgb.json`

`tunedxgb.rds` (Goretzko & Buhner, 2020) stores an xgboost booster serialized by
a 2018 xgboost release in the legacy binary format, which xgboost >= 2.0 can no
longer read. `tunedxgb.json` is the same booster (134 rounds x 8 classes = 1072
trees) re-saved in the modern JSON format. It was produced once, as follows:

1. In R (with `mlr` attached), write the raw booster bytes:
   ```r
   m <- readRDS("tunedxgb.rds")
   writeBin(m$learner.model$next.model$learner.model$raw, "legacy.bin")
   ```
2. In an isolated Python environment with the old release (`pip install xgboost==1.7.6`):
   ```python
   import xgboost as xgb
   b = xgb.Booster(); b.load_model("legacy.bin"); b.save_model("tunedxgb.json")
   ```

Checked: class probabilities from `tunedxgb.json` in R (xgboost 3.2.1.1) match
xgboost 1.7.6 on the legacy file exactly (max abs diff 0) when both use the
booster's best iteration (129 of 134 rounds, selected by early stopping).
