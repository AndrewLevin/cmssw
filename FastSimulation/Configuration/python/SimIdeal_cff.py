import FWCore.ParameterSet.Config as cms

from FastSimulation.Configuration.CommonInputs_cff import *

# Famos SimHits producer
from FastSimulation.EventProducer.FamosSimHits_cff import *

# Gaussian Smearing RecHit producer
from FastSimulation.TrackingRecHitProducer.SiTrackerGaussianSmearingRecHitConverter_cfi import *

# Muon simHit sequence
from FastSimulation.MuonSimHitProducer.MuonSimHitProducer_cfi import *

# propagors from muon reconstruction are used in muon simulation
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *

psim = cms.Sequence(
    famosSimHits+
    MuonSimHits
    )
