
# install the package
devtools::install_github("jaymwin/CAbioacoustics")
library(CAbioacoustics)


# set up database credentials; only do once! ------------------------------

# this will cause 4 separate windows to pop up to enter 1) host, 2) user, 3) password, and 4) db; just copy and paste these from LastPass
# these are then stored in your computer's keyring
cb_set_db_credentials()


# test connection ---------------------------------------------------------

# create database connection from credentials stored in keyring
cb_connect_db()

# see that connection worked and list database tables;
# if this works you are all set!
DBI::dbListTables(conn)

# always disconnect from database when finished
cb_disconnect_db()

# now you can connect next time without setting credentials...