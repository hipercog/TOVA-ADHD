BLOCK
-----
CONTENTS		For every trial (game or transfer), data on performance and coding
NOTES	-	Trials marked bad by the trainer have a code (rejecttype) which was devised in an adhoc way. This might be worth changing to a simpler scheme.
		Note that some trial scores = -1, which is a sign that no game-summary file was found to read the score. This (we believe) happens because the trial is started in CENT, but the EEG acquisition server does not make a connection and the trial is aborted.
		Some trial scores are 0, indicating a game summary file was found but the trial was non-functional.
		For within_session numbers, sometimes the values start at, or increment by more than 1: this is because some sessions contain no trials, and the associated data is all NULL, so no row is written (such NULL sessions are written to the SESSION data file).

patient	-	experiment code for the patient
session_num	-	the number of the DAILY session that these trials come from
within_session	-	number of the CENT session within this day's session
score	-	raw score recorded for a trial
theta_coef	-	trainer-entered coefficient for theta
beta_coef	-	trainer-entered coefficient for beta
gametype	-	name of the game which was played, Empty/Simple Ball, Media (some media trials were Mazes, see trialtype) or AstroComet
gametype_num	-	numeric code for gametype
has_gdf	-	binary code, 0=no EEG data, 1=EEG data is present
gdf_dur_sr250	-	duration of the trial gdf in points; all EEG recorded at sample rate sr=250Hz, gdf_dur/sr = seconds
gdf_good	-	binary code, 0=EEG data was absent or bad, 1=EEG data is good, i.e. long enough to FFT filter and extract band powers, with data recorded in NFB channel 1.
norm_notInv	-	binary code, 0=inverse, 1=normal trials
trialtype	-	1=transfer trial, 2=maze, 3=transfer trial video, 4=transfer trial reading
trainer_says_no	-	binary code 0=trial ok, 1=trial rejected by trainer
rejecttype	-	1=bad signal, 2=crashed mid-block, 3=unusual/untruthful score (e.g. due to failed baseline measurement),
4=subject falling asleep, 5=subject not doing the task, 6=tutorial, 7=score not mentioned in the session diary, 8=info about the specific reason for rejection missing or ambiguous, 9=problem with feedback (e.g. the ball in simple ball game has gotten stuck or does not repond to the changing levels of theta and beta), 11='technical problem', 12=several technical problems mentioned
reject_filter	-	binary code 0=trial ok, 1=trial rejected (all trials with score<1 or marked bad by the trainers)
date_time	-	trial start date and time to 1 second precision

SESSION
-------

NOTES	-	THETA & BETA coefficients were usually set as 1.0, except for three cases 1) early tutorial sessions (normal or inverse), 2) attempts to improve the training at the trainer's discretion, 3) most often, when a session crashed and the trainer wanted to skip baseline recording for the next session, she could enter a coefficient that, when multiplied by the default baseline, would give the value of the day's first recorded baseline.
	-	Since Normal or Inverse protocols are selected when a CENT session is started, all trials in a CENT session will be either Normal or Inverse.

trainer	-	Trainer short name - who was primarily responsible for this patient
patient	-	experiment code for the patient
session_num	-	daily session number
within_session	-	number of the CENT session within this day's session
date_time	-	session start date and time to 1 second precision
raw_score	-	mean of ALL trial scores
basic_score	-	mean of trial scores > 1
TB_ratio	-	beta_base*beta_coef / theta_base*theta_coef
theta_base	-	recorded theta value
beta_base	-	recorded beta value
theta_coef	-	theta coefficient entered by trainer
beta_coef	-	beta coefficient entered by trainer
score	-	mean of trial scores > 1 AND not marked bad by trainer
adj_score	-	adjusted score is the score * TB_ratio; except for inverse sessions, which were score * (theta*coef) / (beta*coef)
adj_transfer_score	-	adjusted score of the transfer trials
excitement	-	self-reported arousal/excitement before the session
hrs_since_sleep	-	self-reported hours since sleeping
hours_slept	-	self-reported duration of last sleep in hours
mood	-	self-reported valence (positive mood) before the session
motivation	-	self-reported motivation before the session
effort	-	self-reported effort after the session
frustration	-	self-reported frustration after the session
trials	-	total number of trials in session
rejected_trials	-	total number of trials in session marked bad by trainer
normal_not_inverse	-	binary code, 0=inverse session, 1=normal session
exclude	-	binary code, 0=1 or more good trials, 1=if all trials are marked bad by trainer!
enobio	-	either 4 or 8 depending on which Enobio was used.

DAILY SESSION
-------------
NOTES	-	Almost all Daily Session fields are the same as Session fields aggregating in the following ways:
		Score, theta/beta, patient condition fields are the mean of matching field from all CENT sessions/day
		Rejected, transfer, normal and inverse trials are the sum of matching field from all CENT sessions/day
		- Fields which don't follow this pattern are explained below.
trainer
substitute_trainer	-	the person who really did the training for that session
patient
session_num
num_CENT_sessions
date_time
raw_score
basic_score
TB_ratio
theta_base
beta_base
theta_coef
beta_coef
score
adj_score
adj_tran_score
mean_excitement
hrs_since_sleep
hours_slept
mean_mood
mean_motivation
mean_effort
mean_frustration
trials
rejected_trials
transfer_trials
inverse_trials
exclude
enobio
Obs_motivaatio	-	observer's estimate from session diary
Obs_asenne	-	observer's estimate from session diary
Obs_keskittymis	-	observer's estimate from session diary
Obs_levottomuus	-	observer's estimate from session diary
Obs_impulsiivisuus	-	observer's estimate from session diary
Obs_turhautumisenkesto	-	observer's estimate from session diary
Obs_itsesaately	-	observer's estimate from session diary
room	-	room number in Malmi

META DATA
---------
NOTES	-	Meta data is aggregated from session data for each patient in the test group.
trainer
patient
num_sessions
num_trials
num_inv_trials
num_trans_trials
adj_normal_score
adj_inverse_score
adj_transfer_score
first_session
last_session
Obs_data_entered	-	Not all observer data from session diaries is recorded. This records how many rows have been entered for each patient.