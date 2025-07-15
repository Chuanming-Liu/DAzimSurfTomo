Input list:
	•	para.in – Parameter file
	•	MODVs.true – Input true Vsv model (same as MOD)
	•	MODGs.true – Input true Gs model (does not contain boundary points); Gs/L (dimensionless)
	•	MODGc.true – Input true Gc model
	•	Surfphase_RV3_5_40s_1s.dat – Blank traveltime data file (only the ray paths are used)

Output list:
	•	surfphase_forward_RV3th.dat – Predicted raypath data file (main output; used in follow-up inversion)
	•	period_Azm_tomo.real – Predicted Rayleigh wave phase velocity map with anisotropy
