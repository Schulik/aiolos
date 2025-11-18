def change_param_entry(PARAM, VALUE, FILE):
	#do sed -i "s/.*REACTION_FILE.*/REACTION_FILE           heliumreac.reac/" runs_2023_co/8me-quicktest_G2-10days-r20-noshad-H2-long-He${a}-kzz-difflimit.par; done
	string = "do sed -i "s/.*REACTION_FILE.*/REACTION_FILE           heliumreac.reac/" runs_2023_co/8me-quicktest_G2-10days-r20-noshad-H2-long-He${a}-kzz-difflimit.par"
