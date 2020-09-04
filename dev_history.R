library(devtools)
library(testthat)

#### Only one time ####

# use_git()

# use_r("calculations")
# use_r("date")

# use_build_ignore("dev_history.R")

# use_gpl3_license("Julie Aubert")

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
attachment::att_to_description()


test()
spell_check()
# spelling::update_wordlist()


check()

goodpractice::goodpractice()
covr::package_coverage()
covr::report()

install()
# build()

# use_github_release()
# usethis::use_version()

#### pkgdown ####

pkgdown::build_site()
# pkgdown::template_navbar()
# pkgdown::template_reference()
pkgdown::clean_site()
