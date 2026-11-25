# Steps

When all tests done and successful and PR is approved by codeowners, follow these steps:

1. Ensure that you update the version in `microSALT/__init__.py`
2. Select "Squash and merge" to merge branch into default branch (master).
3. Deploy master to stage:
   1. Log in to appropriate server `ssh <server.scilifelab.se>`
   2. `us`
   3. Install microSALT:
      ```shell
      bash /home/proj/production/servers/resources/hasta.scilifelab.se/update-microsalt.sh -d stage -b master
      ```
      Add `-a` if the release includes an Alembic migration for the SQL database.
   4. Make sure that the installation was successful
4. Deploy master to production:
   1. Log in to appropriate server `ssh <server.scilifelab.se>`
   2. `up`
   3. Install microSALT:
      ```shell
      bash /home/proj/production/servers/resources/hasta.scilifelab.se/update-microsalt.sh -d production -b master
      ```
      Add `-a` if the release includes an Alembic migration for the SQL database.
   4. Make sure that the installation was successful
