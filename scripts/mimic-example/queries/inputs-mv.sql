-- for mv, we take the endime as the chartime equivalent, but we can't figure out how much volume is associated with all the drug administrations.
SELECT 
  subject_id, 
  hadm_id, 
  icustay_id, 
  endtime AS charttime_equiv, 
  itemid, 
  ROUND(
    CASE
      WHEN mv.amountuom = 'L'
        THEN mv.amount * 1000.0
      WHEN mv.amountuom = 'ml'
        THEN mv.amount
    ELSE NULL END
  ) AS amount,
  'mv' AS orig_table
FROM mimiciii.inputevents_mv mv
WHERE mv.amountuom IN ('ml', 'L');
