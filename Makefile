include .make/Makefile

check_shell:
	shellcheck -x *.pbs

check_r_parse:
	Rscript -e 'd <- lapply(dir(pattern="[.]R$$"), parse)'

check_spelling:
	Rscript -e 'spelling::spell_check_files(c("NEWS", "README.md"), ignore = readLines("inst/WORDLIST"))'
