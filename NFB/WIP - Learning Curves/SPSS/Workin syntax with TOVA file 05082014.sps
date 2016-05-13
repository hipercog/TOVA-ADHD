*Select those where outcome is avalable

USE ALL.
COMPUTE filter_$=(group.2 = 1 | group.2 = 2).
VARIABLE LABELS filter_$ 'group.2 = 1 | group.2 = 2 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

*Run independent t-test between NF and control group for differences in TOVA variables (except Omission errors)

T-TEST GROUPS=group.1(1 2)
  /MISSING=ANALYSIS
  /VARIABLES=changeVAR changeRTM changeCOM changeDP
  /CRITERIA=CI(.95).

*Select those cases where outcome is availbale AND outcome Omission errors is bigger than -1000

DATASET ACTIVATE DataSet3.
USE ALL.
COMPUTE filter_$=((group.2 = 1 | group.2 = 2) & OMSST.2 >  - 1000).
VARIABLE LABELS filter_$ '(group.2 = 1 | group.2 = 2) & OMSST.2 >  - 1000 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

*Run independent t-test on Omission errors

T-TEST GROUPS=group.1(1 2)
  /MISSING=ANALYSIS
  /VARIABLES=changeOM
  /CRITERIA=CI(.95).
