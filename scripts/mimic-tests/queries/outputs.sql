-- First we drop the table if it exists
SELECT 
  oe.icustay_id, 
  oe.charttime, 
  SUM(
    CASE WHEN oe.itemid = 227488 THEN -1*value
    ELSE value END
  ) AS value
FROM mimiciii.outputevents oe
WHERE oe.value < 5000 -- sanity check on urine value
  AND oe.icustay_id IS NOT NULL
GROUP BY icustay_id, charttime;
