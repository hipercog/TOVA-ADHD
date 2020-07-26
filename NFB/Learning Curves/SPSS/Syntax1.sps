* I start out with the SESSION.sav data
* Select only those that are not excluded, not transfer and not inverse

USE ALL.
COMPUTE filter_$=(exclude = 0 & normal_not_inverse = 0 & adj_transfer_score = 0).
VARIABLE LABELS filter_$ 'exclude = 0 & normal_not_inverse = 0 & adj_transfer_score = 0 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

USE ALL.
COMPUTE filter_$=(adj_transfer_score ~= 0).
VARIABLE LABELS filter_$ 'adj_transfer_score ~= 0 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

* Aggregate the adj_score

DATASET DECLARE Aggregated_session.
SORT CASES BY patient session_num.
AGGREGATE
  /OUTFILE='Aggregated_session'
  /PRESORTED
  /BREAK=patient session_num
  /adj_score_mean=MEAN(adj_score) 
  /excitement_mean=MEAN(excitement) 
  /hrs_since_sleep_mean=MEAN(hrs_since_sleep) 
  /hours_slept_mean=MEAN(hours_slept) 
  /mood_mean=MEAN(mood) 
  /motivation_mean=MEAN(motivation) 
  /effort_mean=MEAN(effort) 
  /frustration_mean=MEAN(frustration).


*look at the linear trend

MIXED adj_score_mean WITH session_num
  /CRITERIA=CIN(95) MXITER(100) MXSTEP(10) SCORING(1) SINGULAR(0.000000000001) HCONVERGE(0, 
    ABSOLUTE) LCONVERGE(0, ABSOLUTE) PCONVERGE(0.000001, ABSOLUTE)
  /FIXED=session_num | SSTYPE(3)
  /METHOD=REML
  /PRINT=SOLUTION TESTCOV
  /RANDOM=INTERCEPT session_num | SUBJECT(patient) COVTYPE(UN).

*controls mood

MIXED adj_score_mean WITH session_num mood_mean
  /CRITERIA=CIN(95) MXITER(100) MXSTEP(10) SCORING(1) SINGULAR(0.000000000001) HCONVERGE(0, 
    ABSOLUTE) LCONVERGE(0, ABSOLUTE) PCONVERGE(0.000001, ABSOLUTE)
  /FIXED=session_num | SSTYPE(3)
  /METHOD=REML
  /PRINT=SOLUTION TESTCOV
  /RANDOM=INTERCEPT session_num mood_mean | SUBJECT(patient) COVTYPE(UN).

*control motivation

MIXED adj_score_mean WITH session_num motivation_mean
  /CRITERIA=CIN(95) MXITER(100) MXSTEP(10) SCORING(1) SINGULAR(0.000000000001) HCONVERGE(0, 
    ABSOLUTE) LCONVERGE(0, ABSOLUTE) PCONVERGE(0.000001, ABSOLUTE)
  /FIXED=session_num | SSTYPE(3)
  /METHOD=REML
  /PRINT=SOLUTION TESTCOV
  /RANDOM=INTERCEPT session_num motivation_mean | SUBJECT(patient) COVTYPE(UN).

*d prime

MIXED adj_score_mean WITH session_num DPRSST
  /CRITERIA=CIN(95) MXITER(100) MXSTEP(10) SCORING(1) SINGULAR(0.000000000001) HCONVERGE(0, 
    ABSOLUTE) LCONVERGE(0, ABSOLUTE) PCONVERGE(0.000001, ABSOLUTE)
  /FIXED=session_num | SSTYPE(3)
  /METHOD=REML
  /PRINT=SOLUTION TESTCOV
  /RANDOM=INTERCEPT session_num DPRSST | SUBJECT(patient) COVTYPE(UN).

*age

MIXED adj_score_mean WITH session_num Age
  /CRITERIA=CIN(95) MXITER(100) MXSTEP(10) SCORING(1) SINGULAR(0.000000000001) HCONVERGE(0, 
    ABSOLUTE) LCONVERGE(0, ABSOLUTE) PCONVERGE(0.000001, ABSOLUTE)
  /FIXED=session_num | SSTYPE(3)
  /METHOD=REML
  /PRINT=SOLUTION TESTCOV
  /RANDOM=INTERCEPT session_num Age | SUBJECT(patient) COVTYPE(UN).
*gender

MIXED adj_score_mean BY Gender WITH session_num
  /CRITERIA=CIN(95) MXITER(100) MXSTEP(10) SCORING(1) SINGULAR(0.000000000001) HCONVERGE(0, 
    ABSOLUTE) LCONVERGE(0, ABSOLUTE) PCONVERGE(0.000001, ABSOLUTE)
  /FIXED=session_num Gender | SSTYPE(3)
  /METHOD=REML
  /PRINT=SOLUTION TESTCOV
  /RANDOM=INTERCEPT session_num | SUBJECT(patient) COVTYPE(UN).

*gender with interaction

MIXED adj_score_mean BY Gender WITH session_num
  /CRITERIA=CIN(95) MXITER(100) MXSTEP(10) SCORING(1) SINGULAR(0.000000000001) HCONVERGE(0, 
    ABSOLUTE) LCONVERGE(0, ABSOLUTE) PCONVERGE(0.000001, ABSOLUTE)
  /FIXED=session_num Gender Gender*session_num | SSTYPE(3)
  /METHOD=REML
  /PRINT=SOLUTION TESTCOV
  /RANDOM=INTERCEPT session_num | SUBJECT(patient) COVTYPE(UN).

*education

MIXED adj_score_mean BY Education WITH session_num
  /CRITERIA=CIN(95) MXITER(100) MXSTEP(10) SCORING(1) SINGULAR(0.000000000001) HCONVERGE(0, 
    ABSOLUTE) LCONVERGE(0, ABSOLUTE) PCONVERGE(0.000001, ABSOLUTE)
  /FIXED=session_num Education Education*session_num | SSTYPE(3)
  /METHOD=REML
  /PRINT=SOLUTION TESTCOV
  /RANDOM=INTERCEPT session_num | SUBJECT(patient) COVTYPE(UN).

*ADHD/ADD

MIXED adj_score_mean BY ADHD_ADD WITH session_num
  /CRITERIA=CIN(95) MXITER(100) MXSTEP(10) SCORING(1) SINGULAR(0.000000000001) HCONVERGE(0, 
    ABSOLUTE) LCONVERGE(0, ABSOLUTE) PCONVERGE(0.000001, ABSOLUTE)
  /FIXED=session_num ADHD_ADD ADHD_ADD*session_num | SSTYPE(3)
  /METHOD=REML
  /PRINT=SOLUTION TESTCOV
  /RANDOM=INTERCEPT session_num | SUBJECT(patient) COVTYPE(UN).

*comorb diag

MIXED adj_score_mean BY Comorb_diag WITH session_num
  /CRITERIA=CIN(95) MXITER(100) MXSTEP(10) SCORING(1) SINGULAR(0.000000000001) HCONVERGE(0, 
    ABSOLUTE) LCONVERGE(0, ABSOLUTE) PCONVERGE(0.000001, ABSOLUTE)
  /FIXED=session_num Comorb_diag Comorb_diag*session_num | SSTYPE(3)
  /METHOD=REML
  /PRINT=SOLUTION TESTCOV
  /RANDOM=INTERCEPT session_num | SUBJECT(patient) COVTYPE(UN).

*comorb scale

MIXED adj_score_mean BY Comorb_scales WITH session_num
  /CRITERIA=CIN(95) MXITER(100) MXSTEP(10) SCORING(1) SINGULAR(0.000000000001) HCONVERGE(0, 
    ABSOLUTE) LCONVERGE(0, ABSOLUTE) PCONVERGE(0.000001, ABSOLUTE)
  /FIXED=session_num Comorb_scales Comorb_scales*session_num | SSTYPE(3)
  /METHOD=REML
  /PRINT=SOLUTION TESTCOV
  /RANDOM=INTERCEPT session_num | SUBJECT(patient) COVTYPE(UN).









