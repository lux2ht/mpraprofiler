##
##  Makefile tasks for building a Markdown version of the vignette (for GitHub)
##
##  It's named 'GNUmakefile' so that devtools::build_vignettes will ignore it,
##  but 'make' will still work for the maintainers.
##
##  Author:   Kevin Ernst <kevin.ernst -at- cchmc.org>
##  Date:     19 October 2022
##  

# set to 1 to suppress rmarkdown::render's output
QUIET =
INPUTFILE = sample_analysis.Rmd
OUTPUTBASE = sample_analysis
RSCRIPT = $(R_HOME)/bin/Rscript
# the name of this Makefile
THIS = $(word $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST))


help:  # prints this help
	@perl -e "$$PRINT_HELP_PERL_SCRIPT" $(THIS)

render: $(OUTPUTBASE).md  # render Markdown output

$(OUTPUTBASE).md: $(INPUTFILE) ../R/*
	$(RSCRIPT) -e "library(rmarkdown); rmarkdown::render('$<', output_file='$@', output_format='github_document', quiet=$(if $(QUIET),TRUE,FALSE))"

# Note: if this file is called 'Makefile', devtools::build_vignettes will run
# this target, which is probably not what you want.
clean:  # clean all generated files
	-rm *.md *.pdf *.html
	-rm -r *_files/


##
##  Perl script to generate the 'make help' output
##
define PRINT_HELP_PERL_SCRIPT
	use Term::ANSIColor qw(:constants);
	$$max = 0;
	@targets = ();
	print "\n  ", UNDERLINE, "Makefile targets", RESET, "\n\n";
	while (<>) {
		push @targets, [$$1, $$2] if /^(\w.+):.*#\s*(.*)/;
		$$max = length($$1) if length($$1) > $$max;
	}
	foreach (@targets) {
		printf "    %s%smake %-$${max}s%s    %s\n", BOLD, BLUE, @$$_[0], RESET, @$$_[1];
	}
	print "\n";
endef
export PRINT_HELP_PERL_SCRIPT
