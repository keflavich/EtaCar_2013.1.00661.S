
import analysisUtils as aU
es = aU.stuffForScienceDataReduction()

importasdm(asdm='raw/uid___A002_X9baf64_X5dc', vis='uid___A002_X9baf64_X5d_raw.ms')

# could not correct antenna positions
es.generateReducScript('uid___A002_X9baf64_X5d_raw.ms', skipSyscalChecks=True,
                       lazy=True, corrAntPos=False)
