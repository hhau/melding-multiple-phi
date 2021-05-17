SELECT
  oe.icustay_id,
  oe.charttime,
  CASE
    WHEN oe.itemid = 227488 THEN -1 * value
    ELSE value
  END AS value,
  d.label,
  d.dbsource AS orig_table
FROM mimiciii.outputevents oe
INNER JOIN mimiciii.d_items d
  ON oe.itemid = d.itemid
WHERE oe.value < 5000 -- sanity check on urine value
  AND oe.icustay_id IS NOT NULL;
