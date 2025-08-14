## Description

<!-- Summary of the changes made. If not self-evident, mention what prompted the change. -->

### Primary function of PR

- [ ] Hot-fix
- [ ] Patch
- [ ] Minor functionality improvement
- [ ] New type of analysis
- [ ] Backward-breaking functionality improvement
- [ ] This change requires internal documents to be updated
- [ ] This change requires another repository to be updated

## Testing

<!-- If the update is a hot-fix, it is sufficient to rely on the development testing along with the Travis self-test automatically applied to the PR. -->

- `bash /home/proj/production/servers/resources/hasta.scilifelab.se/install-microsalt-stage.sh BRANCHNAME`
- `us`
- `conda activate S_microSALT`
- `microSALT analyse --input /path/to/fastq/ SAMPLEINFO_FILE`

### Test results

_These are the results of the tests, and necessary conclusions, that prove the stability of the PR._

## Sign-offs

- [ ] Approved to run at Clinical-Genomics by @karlnyr or @Clinical-Genomics/micro
