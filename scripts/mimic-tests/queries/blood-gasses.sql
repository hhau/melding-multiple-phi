/*CREATE OR REPLACE FUNCTION DATETIME_ADD(datetime_val TIMESTAMP(3), intvl INTERVAL) RETURNS TIMESTAMP(3) AS $$
BEGIN
RETURN datetime_val + intvl;
END; $$
LANGUAGE PLPGSQL;

--  DATETIME_SUB(datetime, INTERVAL 'n' DATEPART) -> datetime - INTERVAL 'n' DATEPART
CREATE OR REPLACE FUNCTION DATETIME_SUB(datetime_val TIMESTAMP(3), intvl INTERVAL) RETURNS TIMESTAMP(3) AS $$
BEGIN
RETURN datetime_val - intvl;
END; $$
LANGUAGE PLPGSQL;

CREATE OR REPLACE FUNCTION DATETIME_DIFF(endtime TIMESTAMP(3), starttime TIMESTAMP(3), datepart TEXT) RETURNS NUMERIC AS $$
BEGIN
RETURN 
    EXTRACT(EPOCH FROM endtime - starttime) /
    CASE
        WHEN datepart = 'SECOND' THEN 1.0
        WHEN datepart = 'MINUTE' THEN 60.0
        WHEN datepart = 'HOUR' THEN 3600.0
        WHEN datepart = 'DAY' THEN 24*3600.0
        WHEN datepart = 'YEAR' THEN 365.242*24*3600.0
    ELSE NULL END;
END; $$
LANGUAGE PLPGSQL;*/

with pvt as
( -- begin query that extracts the data
  select ie.subject_id, ie.hadm_id, ie.icustay_id
  -- here we assign labels to ITEMIDs
  -- this also fuses together multiple ITEMIDs containing the same data
      , case
        when itemid = 50800 then 'SPECIMEN'
        when itemid = 50816 then 'FIO2'
        when itemid = 50821 then 'PO2'
        else null
        end as label
        , charttime
        , value
        -- add in some sanity checks on the values
        , case
          -- ensure FiO2 is a valid number between 21-100
          -- mistakes are rare (<100 obs out of ~100,000)
          -- there are 862 obs of valuenum == 20 - some people round down!
          -- rather than risk imputing garbage data for FiO2, we simply NULL invalid values
          when itemid = 50816 and valuenum < 20 then null
          when itemid = 50816 and valuenum > 100 then null
          when itemid = 50821 and valuenum > 800 then null -- PO2
           -- conservative upper limit
        else valuenum
        end as valuenum

    FROM  mimiciii.icustays ie
    left join  mimiciii.labevents le
      on le.subject_id = ie.subject_id and le.hadm_id = ie.hadm_id
      and le.charttime between (DATETIME_SUB(ie.intime, INTERVAL '6' HOUR)) and (DATETIME_ADD(ie.intime, INTERVAL '14' DAY))
      and le.ITEMID in
      -- blood gases
      (
        50800, 50816, 50821
      )
)
select pvt.SUBJECT_ID, pvt.HADM_ID, pvt.ICUSTAY_ID, pvt.CHARTTIME
, max(case when label = 'SPECIMEN' then value else null end) as specimen
, max(case when label = 'FIO2' then valuenum else null end) as fio2
, max(case when label = 'PO2' then valuenum else null end) as po2
from pvt
group by pvt.subject_id, pvt.hadm_id, pvt.icustay_id, pvt.CHARTTIME
order by pvt.subject_id, pvt.hadm_id, pvt.icustay_id, pvt.CHARTTIME;