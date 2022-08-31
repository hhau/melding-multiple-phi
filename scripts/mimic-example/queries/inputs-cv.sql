-- for cv, we can just pull the amount measurements and the charttime corresponds to the time the fluid was delivered by.
SELECT
  cv.subject_id,
  cv.hadm_id,
  cv.icustay_id,
  cv.charttime AS charttime_equiv,
  cv.itemid,
  ROUND(cv.amount) AS amount,
  cv.amountuom,
  d.label,
  d.dbsource AS orig_table
FROM mimiciii.inputevents_cv cv
INNER JOIN mimiciii.d_items d
  ON cv.itemid = d.itemid
WHERE cv.amountuom = 'ml';
