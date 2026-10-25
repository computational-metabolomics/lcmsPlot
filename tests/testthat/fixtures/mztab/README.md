`MTBLS263-trimmed.mztab` is `examples/2_0-Metabolomics_Release/MTBLS263.mztab`
from the HUPO-PSI mzTab repository (https://github.com/HUPO-PSI/mzTab), trimmed
to the SMF features 6, 7, 13 and 937, the SML rows that reference them, and their
SME evidence. One SML row was added: `529` (isoleucine) references the same
feature as leucine (`528`), so the fixture covers a feature with ambiguous
identifications.
