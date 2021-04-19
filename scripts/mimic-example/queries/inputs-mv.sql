-- for mv, we take the endime as the chartime equivalent, but we can't figure out how much volume is associated with all the drug administrations.
SELECT
  mv.subject_id,
  mv.hadm_id,
  mv.icustay_id,
  mv.endtime AS charttime_equiv,
  mv.itemid,
  ROUND(
    CASE
      WHEN mv.amountuom = 'L'
        THEN mv.amount * 1000.0
      WHEN mv.amountuom = 'ml'
        THEN mv.amount
    ELSE NULL END
  ) AS amount,
  'ml' AS amountuom,
  d.label,
  d.dbsource AS orig_table
FROM mimiciii.inputevents_mv mv
INNER JOIN mimiciii.d_items d
  ON mv.itemid = d.itemid
WHERE mv.amountuom IN ('ml', 'L');
