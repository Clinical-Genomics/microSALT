# Description
_The features of this PR primarily concerns end-users/bioinformaticians/internals_

_This is a summary of the changes made._

_If not self-evident, also mention what prompted the change._

## Primary function of PR
- [ ] Hotfix
- [ ] Patch
- [ ] Minor functionality improvement
- [ ] New type of analysis
- [ ] Backward-breaking functionality improvement
- [ ] This change requires internal documents to be updated
- [ ] This change requires another repository to be updated

# Testing
_This is a description of the tests necessary to verify the stability of the PR._

_Basic test routine:_
- _`bash /home/proj/production/servers/resources/hasta.scilifelab.se/update-microsalt-stage.sh BRANCHNAME`_
- _`us`_
- _`source activate S_microSALT`_
- _(SITUATIONAL) `export MICROSALT_CONFIG=/home/proj/dropbox/microSALT.json`_
- _Select a relevant subset of the following:_
- _`microSALT analyse project MIC3109`_
- _`microSALT analyse project MIC4107`_
- _`microSALT analyse project MIC4109`_
- _`microSALT analyse project ACC5551`_

_Verify that the results for projects MIC3109, MIC4107, MIC4109 & ACC5551 are consistent with the results attached to AMSystem doc 1490, Microbial_WGS.xlsx_

## Test results
_These are the results of the tests, and necessary conclusions, that prove the stability of the PR._

# Sign-offs
- [ ] Code tested by @sylvinite
- [ ] Approved to run at Clinical-Genomics by @octocat
