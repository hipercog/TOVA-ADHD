*Otetaan vain enobio 8/ We take enobio 8


USE ALL.
COMPUTE filter_$=(enobio1_first = 8).
VARIABLE LABELS filter_$ 'enobio1_first = 8 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

*Ajetaan mixed models/ run the mixed model

MIXED adj_score_mean WITH session_num
  /CRITERIA=CIN(95) MXITER(100) MXSTEP(10) SCORING(1) SINGULAR(0.000000000001) HCONVERGE(0, 
    ABSOLUTE) LCONVERGE(0, ABSOLUTE) PCONVERGE(0.000001, ABSOLUTE)
  /FIXED=session_num | SSTYPE(3)
  /METHOD=REML
  /PRINT=SOLUTION TESTCOV
  /RANDOM=INTERCEPT session_num | SUBJECT(patient) COVTYPE(UN).

* Ei konvergenssia: Otetaan vain kun enobio 4/ no convergence achieved, now take enobio 4

USE ALL.
COMPUTE filter_$=(enobio1_first = 4).
VARIABLE LABELS filter_$ 'enobio1_first = 4 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

MIXED adj_score_mean WITH session_num
  /CRITERIA=CIN(95) MXITER(100) MXSTEP(10) SCORING(1) SINGULAR(0.000000000001) HCONVERGE(0, 
    ABSOLUTE) LCONVERGE(0, ABSOLUTE) PCONVERGE(0.000001, ABSOLUTE)
  /FIXED=session_num | SSTYPE(3)
  /METHOD=REML
  /PRINT=SOLUTION TESTCOV
  /RANDOM=INTERCEPT session_num | SUBJECT(patient) COVTYPE(UN).

* Sama juttu, ei konvergenssia/ no convergence achieved



