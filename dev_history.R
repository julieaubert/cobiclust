library(devtools)
library(testthat)
#library(formatR)
#### Only one time ####

# use_git()

# use_r("calculations")
# use_r("date")

# use_build_ignore("dev_history.R")

# use_gpl3_license("Julie Aubert")
# Define License with use_*_license()
#usethis::use_mit_license("Julie Aubert")
# usethis::use_gpl_license(version = 3, include_future = TRUE)


# use_github()
# use_github_links()

# use_testthat()
# use_test("cobiclust")
# use_test("date")

# use_spell_check()

# use_readme_rmd()


# use_lifecycle_badge("maturing")
# badger::badge_last_commit()

# use_news_md()

# use_github_action_check_release()

# use_package_doc()

#### Used regularly ####

load_all()

document()

use_tidy_description()
attachment::att_amend_desc()


test()
spell_check()
spelling::update_wordlist()

# Pour mettre son code au bon format

 #tidy_dir("~/git/cobiclust/R/", width.cutoff = I(80))
# Ise à jour des dependances pour nous
# devtools::install_deps()

check()

goodpractice::goodpractice()
covr::package_coverage()
covr::report()

# install()
# build() install and restart à la place

# use_github_release()
# usethis::use_version()

#### pkgdown ####

pkgdown::build_site()
# pkgdown::template_navbar()
# pkgdown::template_reference()
pkgdown::clean_site()
