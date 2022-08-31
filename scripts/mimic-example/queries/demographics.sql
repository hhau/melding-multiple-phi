SELECT
  ie.icustay_id,
  p.gender,
  ROUND(DATETIME_DIFF(ie.intime, p.dob, 'YEAR')) AS age_at_icu_adm
FROM mimiciii.patients p
LEFT JOIN mimiciii.icustays ie
  ON p.subject_id = ie.subject_id;