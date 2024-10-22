# Git setup with new computer

## Generate PAT
usethis::create_github_token()

## Stash the token
gitcreds::gitcreds_set()

## Create Rproj from Git Repo
usethis::create_from_github("https://github.com/skeyser/Sierra_Biodiv.git",
                            destdir = "C:/Users/srk252/Documents/Rprojs/")
