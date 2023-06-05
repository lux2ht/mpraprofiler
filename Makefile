##
##  Makefile tasks for building the vignette and README as Markdown
##
##  Author:   Kevin Ernst <kevin.ernst -at- cchmc.org>
##  Date:     19 October 2022, updated 3 June 2023
##  

NAME = mpraprofiler
# set to 1 to suppress rmarkdown::render's output
QUIET =
# the name of this Makefile
THIS = $(word $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST))


help:  # prints this help
	@perl -e "$$AUTOGENERATE_HELP_PL" $(THIS)

render: vignette readme  # render Markdown versions of README and vignette(s)

readme: README.md  # render Markdown version of README

vignettes: vignette
vignette: vignettes/sample_analysis.md R/*  # render Markdown version of vignette(s)

%.md: %.Rmd
	Rscript -e '.libPaths("$(PWD)/R"); library(rmarkdown); rmarkdown::render("$<", output_file="$@", output_format="github_document", quiet=$(if $(QUIET),TRUE,FALSE))'

clean:  # clean intermediate files
	-rm README.html
	-rm vignettes/*.html
	-rm -r vignettes/*_files/

reallyclean: clean  # clean all generated files
	-rm README.md vignettes/*.md


##
##  Perl script to generate the 'make help' output
##
define AUTOGENERATE_HELP_PL
	use Term::ANSIColor qw(:constants);
	$$max = 0;
	@targets = ();
	print "\n  ", UNDERLINE, "Makefile targets - $(NAME)", RESET, "\n\n";
	while (<>) {
		push @targets, [$$1, $$2] if /^(\w.+):.*#\s*(.*)/;
		$$max = length($$1) if length($$1) > $$max;
	}
	foreach (@targets) {
		printf "    %s%smake %-$${max}s%s    %s\n", BOLD, BLUE, @$$_[0], RESET, @$$_[1];
	}
	print "\n";
endef
export AUTOGENERATE_HELP_PL
