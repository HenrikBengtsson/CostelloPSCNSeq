check: check_shell check_r_parse

check_shell:
	shellcheck -x *.pbs

check_r_parse:
	Rscript -e 'd <- lapply(dir(pattern="[.]R$$"), parse)'
