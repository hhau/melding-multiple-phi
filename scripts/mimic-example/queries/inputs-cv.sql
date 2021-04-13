-- for cv, we can just pull the amount measurements and the charttime corresponds to the time the fluid was delivered by.
SELECT
  subject_id,
  hadm_id, 
  icustay_id, 
  charttime AS charttime_equiv, 
  itemid,
  ROUND(amount) AS amount,
  'cv' AS orig_table
FROM mimiciii.inputevents_cv;
