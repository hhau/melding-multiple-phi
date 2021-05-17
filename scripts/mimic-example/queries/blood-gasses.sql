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
    select
        mimiciii.icustays.subject_id,
        mimiciii.icustays.hadm_id,
        mimiciii.icustays.icustay_id,
        -- here we assign labels to ITEMIDs
        -- this also fuses together multiple ITEMIDs containing the same data
        charttime,
        value,
        case
            when itemid = 50800 then 'SPECIMEN'
            when itemid = 50801 then 'AADO2'
            when itemid = 50802 then 'BASEEXCESS'
            when itemid = 50803 then 'BICARBONATE'
            when itemid = 50804 then 'TOTALCO2'
            when itemid = 50805 then 'CARBOXYHEMOGLOBIN'
            when itemid = 50806 then 'CHLORIDE'
            when itemid = 50808 then 'CALCIUM'
            when itemid = 50809 then 'GLUCOSE'
            when itemid = 50810 then 'HEMATOCRIT'
            when itemid = 50811 then 'HEMOGLOBIN'
            when itemid = 50812 then 'INTUBATED'
            when itemid = 50813 then 'LACTATE'
            when itemid = 50814 then 'METHEMOGLOBIN'
            when itemid = 50815 then 'O2FLOW'
            when itemid = 50816 then 'FIO2'
            when itemid = 50817 then 'SO2' -- OXYGENSATURATION
            when itemid = 50818 then 'PCO2'
            when itemid = 50819 then 'PEEP'
            when itemid = 50820 then 'PH'
            when itemid = 50821 then 'PO2'
            when itemid = 50822 then 'POTASSIUM'
            when itemid = 50823 then 'REQUIREDO2'
            when itemid = 50824 then 'SODIUM'
            when itemid = 50825 then 'TEMPERATURE'
            when itemid = 50826 then 'TIDALVOLUME'
            when itemid = 50827 then 'VENTILATIONRATE'
            when itemid = 50828 then 'VENTILATOR'
        end as label,
        -- add in some sanity checks on the values
        case
            -- allow negative baseexcess
            when valuenum <= 0 and itemid != 50802 then null
            when itemid = 50810 and valuenum > 100 then null -- hematocrit
            -- ensure FiO2 is a valid number between 21-100
            -- mistakes are rare (<100 obs out of ~100,000)
            -- there are 862 obs of valuenum == 20 - some people round down!
            -- rather than risk imputing garbage data for FiO2, we simply NULL invalid values
            when itemid = 50816 and valuenum < 20 then null
            when itemid = 50816 and valuenum > 100 then null
            when itemid = 50817 and valuenum > 100 then null -- O2 sat
            when itemid = 50815 and valuenum > 70 then null -- O2 flow
            when itemid = 50821 and valuenum > 800 then null -- PO2
            -- conservative upper limit
            else valuenum
        end as valuenum

    from mimiciii.icustays
    left join mimiciii.labevents
        on
            mimiciii.labevents.subject_id = mimiciii.icustays.subject_id and
            mimiciii.labevents.hadm_id = mimiciii.icustays.hadm_id and  
            mimiciii.labevents.charttime between 
                (DATETIME_SUB(mimiciii.icustays.intime, interval '6' hour)) and
                (DATETIME_ADD(mimiciii.icustays.intime, interval '30' day)) -- ***** HERE IS WHERE WE SPECIFY HOW MANY DAYS WE ARE INTERESTED IN ***** --
            and mimiciii.labevents.itemid in
            -- blood gases
            (
                50800,
                50801,
                50802,
                50803,
                50804,
                50805,
                50806,
                50807,
                50808,
                50809,
                50810,
                50811,
                50812,
                50813,
                50814,
                50815,
                50816,
                50817,
                50818,
                50819,
                50820, 
                50821,
                50822,
                50823,
                50824,
                50825,
                50826,
                50827,
                50828,
                51545
            )
),

bg as
(
    select pvt.subject_id, pvt.hadm_id, pvt.icustay_id, pvt.charttime,
        MAX(case when label = 'SPECIMEN' then value end) as specimen,
        MAX(case when label = 'AADO2' then valuenum end) as aado2,
        MAX(case when label = 'BASEEXCESS' then valuenum end) as baseexcess,
        MAX(case when label = 'BICARBONATE' then valuenum end) as bicarbonate,
        MAX(case when label = 'TOTALCO2' then valuenum end) as totalco2,
        MAX(case when label = 'CARBOXYHEMOGLOBIN' then valuenum end) as carboxyhemoglobin,
        MAX(case when label = 'CHLORIDE' then valuenum end) as chloride,
        MAX(case when label = 'CALCIUM' then valuenum end) as calcium,
        MAX(case when label = 'GLUCOSE' then valuenum end) as glucose,
        MAX(case when label = 'HEMATOCRIT' then valuenum end) as hematocrit,
        MAX(case when label = 'HEMOGLOBIN' then valuenum end) as hemoglobin,
        MAX(case when label = 'INTUBATED' then valuenum end) as intubated,
        MAX(case when label = 'LACTATE' then valuenum end) as lactate,
        MAX(case when label = 'METHEMOGLOBIN' then valuenum end) as methemoglobin,
        MAX(case when label = 'O2FLOW' then valuenum end) as o2flow,
        MAX(case when label = 'FIO2' then valuenum end) as fio2,
        -- OXYGENSATURATION
        MAX(case when label = 'SO2' then valuenum end) as so2,
        MAX(case when label = 'PCO2' then valuenum end) as pco2,
        MAX(case when label = 'PEEP' then valuenum end) as peep,
        MAX(case when label = 'PH' then valuenum end) as ph,
        MAX(case when label = 'PO2' then valuenum end) as po2,
        MAX(case when label = 'POTASSIUM' then valuenum end) as potassium,
        MAX(case when label = 'REQUIREDO2' then valuenum end) as requiredo2,
        MAX(case when label = 'SODIUM' then valuenum end) as sodium,
        MAX(case when label = 'TEMPERATURE' then valuenum end) as temperature,
        MAX(case when label = 'TIDALVOLUME' then valuenum end) as tidalvolume,
        MAX(case when label = 'VENTILATIONRATE' then valuenum end) as ventilationrate,
        MAX(case when label = 'VENTILATOR' then valuenum end) as ventilator
    from pvt
    group by pvt.subject_id, pvt.hadm_id, pvt.icustay_id, pvt.charttime
    order by pvt.subject_id, pvt.hadm_id, pvt.icustay_id, pvt.charttime
),

stg_spo2 as
(
    select subject_id, hadm_id, icustay_id, charttime,
        -- max here is just used to group SpO2 by charttime
        MAX(
            case
                when valuenum <= 0 or valuenum > 100 then null else valuenum
            end
        ) as spo2
    from mimiciii.chartevents
    -- o2 sat
    where itemid in
        (
            646, -- SpO2,
            220277 -- O2 saturation pulseoxymetry
        )
    group by subject_id, hadm_id, icustay_id, charttime
),

stg_fio2 as
(
    select subject_id, hadm_id, icustay_id, charttime,
        -- pre-process the FiO2s to ensure they are between 21-100%
        MAX(
            case
                when itemid = 223835
                     then case
                when valuenum > 0 and valuenum <= 1
                      then valuenum * 100
                -- improperly input data - looks like O2 flow in litres
                when valuenum > 1 and valuenum < 21
                      then null
                when valuenum >= 21 and valuenum <= 100
                      then valuenum end -- unphysiological
                when itemid in (3420, 3422)
                -- all these values are well formatted
                     then valuenum
                when itemid = 190 and valuenum > 0.20 and valuenum < 1
                -- well formatted but not in %
                     then valuenum * 100 end
        ) as fio2_chartevents
    from mimiciii.chartevents
    where itemid in
        (
            3420, -- FiO2,
            190, -- FiO2 set,
            223835, -- Inspired O2 Fraction (FiO2),
            3422 -- FiO2 [measured]
        )
    -- exclude rows marked as error
    and (error is null or error = 0)
    group by subject_id, hadm_id, icustay_id, charttime
),

stg2 as
(
    select bg.*,
        stg_spo2.spo2,
        ROW_NUMBER() over (
             partition by
            bg.icustay_id, bg.charttime
             order by stg_spo2.charttime desc
        ) as lastrowspo2
    from bg
    left join stg_spo2
        -- same patient
        on bg.icustay_id = stg_spo2.icustay_id
    -- spo2 occurred at most 2 hours before this blood gas
             and stg_spo2.charttime >= DATETIME_SUB(
                bg.charttime, interval '2' hour
            )
             and stg_spo2.charttime <= bg.charttime
    where bg.po2 is not null
),

stg3 as
(
    select stg2.*,
        stg_fio2.fio2_chartevents,
        ROW_NUMBER() over (
             partition by
            stg2.icustay_id, stg2.charttime
             order by stg_fio2.charttime desc
        ) as lastrowfio2,

        -- create our specimen prediction
        1 / (1 + EXP(-(-0.02544
                        + 0.04598 * po2
                     + COALESCE(-0.15356 * spo2, -0.15356 * 97.49420 + 0.13429)
                         + COALESCE(
                0.00621 * fio2_chartevents, 0.00621 * 51.49550 + -0.24958
            )
                 + COALESCE( 0.10559 * hemoglobin, 0.10559 * 10.32307 + 0.05954)
                 + COALESCE( 0.13251 * so2, 0.13251 * 93.66539 + -0.23172)
                 + COALESCE(-0.01511 * pco2, -0.01511 * 42.08866 + -0.01630)
                 + COALESCE( 0.01480 * fio2, 0.01480 * 63.97836 + -0.31142)
         + COALESCE(-0.00200 * aado2, -0.00200 * 442.21186 + -0.01328)
            + COALESCE(-0.03220 * bicarbonate, -0.03220 * 22.96894 + -0.06535)
            + COALESCE( 0.05384 * totalco2, 0.05384 * 24.72632 + -0.01405)
            + COALESCE( 0.08202 * lactate, 0.08202 * 3.06436 + 0.06038)
            + COALESCE( 0.10956 * ph, 0.10956 * 7.36233 + -0.00617)
            + COALESCE( 0.00848 * o2flow, 0.00848 * 7.59362 + -0.35803)
                ))) as specimen_prob
    from stg2
    left join stg_fio2
        -- same patient
        on stg2.icustay_id = stg_fio2.icustay_id
    -- fio2 occurred at most 4 hours before this blood gas
             and stg_fio2.charttime between DATETIME_SUB(
            stg2.charttime, interval '4' hour
        ) and stg2.charttime
    -- only the row with the most recent SpO2 (if no SpO2 found lastRowSpO2 = 1)
    where stg2.lastrowspo2 = 1
)

select subject_id, hadm_id,
    icustay_id, charttime,
    specimen -- raw data indicating sample type, only present 80% of the time,

    -- prediction of specimen for missing data
    specimen_prob,
    so2,

    -- oxygen related parameters
    spo2, po2 -- note spo2 is FROM mimiciii.chartevents,
    pco2, fio2_chartevents,
    fio2, aado2,
    ph,
    -- also calculate AADO2
    baseexcess,
    bicarbonate,
    -- acid-base parameters
    totalco2, hematocrit,
    hemoglobin, carboxyhemoglobin,

    -- blood count parameters
    methemoglobin,
    chloride,
    calcium,
    temperature,

    -- chemistry
    potassium, sodium,
    lactate,
    glucose, intubated,
    tidalvolume,
    ventilationrate,

    -- ventilation stuff that's sometimes input
    ventilator, peep, o2flow, requiredo2,
    case
         when specimen is not null then specimen
         when specimen_prob > 0.75 then 'ART' 
    end as specimen_pred,
    case
         when po2 is not null
            and pco2 is not null
            and COALESCE(fio2, fio2_chartevents) is not null
            -- multiple by 100 because FiO2 is in a % but should be a fraction
            then (
                 COALESCE(fio2, fio2_chartevents) / 100
            ) * (760 - 47) - (pco2 / 0.8) - po2
    end as aado2_calc,
    case
         when po2 is not null and COALESCE(fio2, fio2_chartevents) is not null
        -- multiply by 100 because FiO2 is in a % but should be a fraction
        then 100 * po2 / (COALESCE(fio2, fio2_chartevents))
    end as pao2fio2

from stg3
where lastrowfio2 = 1 -- only the most recent FiO2
    -- restrict it to *only* arterial samples
    and (specimen = 'ART' or specimen_prob > 0.75)
order by icustay_id, charttime;
