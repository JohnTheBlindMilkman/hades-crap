#include "TFile.h"
#include "TCutG.h"
#include "TH2F.h"
#include <iostream>
#include <vector>
using namespace std;

//#warning here are latest cuts with gCalor and new TOF digi

TFile *cutfile_betamom_pionCmom = new TFile("/lustre/hades/user/tscheib/apr12/ID_Cuts/BetaMomIDCuts_PionsProtons_gen8_DATA_RK400_PionConstMom.root");

TCutG* betamom_3sig_pim_tof_pionCmom=(TCutG*)cutfile_betamom_pionCmom->Get("BetaCutPiM_TOF_3.0");
TCutG* betamom_3sig_pim_rpc_pionCmom=(TCutG*)cutfile_betamom_pionCmom->Get("BetaCutPiM_RPC_3.0");
TCutG* betamom_3sig_pip_tof_pionCmom=(TCutG*)cutfile_betamom_pionCmom->Get("BetaCutPiP_TOF_3.0");
TCutG* betamom_3sig_pip_rpc_pionCmom=(TCutG*)cutfile_betamom_pionCmom->Get("BetaCutPiP_RPC_3.0");
TCutG* betamom_3sig_p_tof_pionCmom=(TCutG*)cutfile_betamom_pionCmom->Get("BetaCutProton_TOF_3.0");
TCutG* betamom_3sig_p_rpc_pionCmom=(TCutG*)cutfile_betamom_pionCmom->Get("BetaCutProton_RPC_3.0");
TCutG* betamom_25sig_pim_tof_pionCmom=(TCutG*)cutfile_betamom_pionCmom->Get("BetaCutPiM_TOF_2.5");
TCutG* betamom_25sig_pim_rpc_pionCmom=(TCutG*)cutfile_betamom_pionCmom->Get("BetaCutPiM_RPC_2.5");
TCutG* betamom_25sig_pip_tof_pionCmom=(TCutG*)cutfile_betamom_pionCmom->Get("BetaCutPiP_TOF_2.5");
TCutG* betamom_25sig_pip_rpc_pionCmom=(TCutG*)cutfile_betamom_pionCmom->Get("BetaCutPiP_RPC_2.5");
TCutG* betamom_25sig_p_tof_pionCmom=(TCutG*)cutfile_betamom_pionCmom->Get("BetaCutProton_TOF_2.5");
TCutG* betamom_25sig_p_rpc_pionCmom=(TCutG*)cutfile_betamom_pionCmom->Get("BetaCutProton_RPC_2.5");
TCutG* betamom_2sig_pim_tof_pionCmom=(TCutG*)cutfile_betamom_pionCmom->Get("BetaCutPiM_TOF_2.0");
TCutG* betamom_2sig_pim_rpc_pionCmom=(TCutG*)cutfile_betamom_pionCmom->Get("BetaCutPiM_RPC_2.0");
TCutG* betamom_2sig_pip_tof_pionCmom=(TCutG*)cutfile_betamom_pionCmom->Get("BetaCutPiP_TOF_2.0");
TCutG* betamom_2sig_pip_rpc_pionCmom=(TCutG*)cutfile_betamom_pionCmom->Get("BetaCutPiP_RPC_2.0");
TCutG* betamom_2sig_p_tof_pionCmom=(TCutG*)cutfile_betamom_pionCmom->Get("BetaCutProton_TOF_2.0");
TCutG* betamom_2sig_p_rpc_pionCmom=(TCutG*)cutfile_betamom_pionCmom->Get("BetaCutProton_RPC_2.0");
TCutG* betamom_15sig_pim_tof_pionCmom=(TCutG*)cutfile_betamom_pionCmom->Get("BetaCutPiM_TOF_1.5");
TCutG* betamom_15sig_pim_rpc_pionCmom=(TCutG*)cutfile_betamom_pionCmom->Get("BetaCutPiM_RPC_1.5");
TCutG* betamom_15sig_pip_tof_pionCmom=(TCutG*)cutfile_betamom_pionCmom->Get("BetaCutPiP_TOF_1.5");
TCutG* betamom_15sig_pip_rpc_pionCmom=(TCutG*)cutfile_betamom_pionCmom->Get("BetaCutPiP_RPC_1.5");
TCutG* betamom_15sig_p_tof_pionCmom=(TCutG*)cutfile_betamom_pionCmom->Get("BetaCutProton_TOF_1.5");
TCutG* betamom_15sig_p_rpc_pionCmom=(TCutG*)cutfile_betamom_pionCmom->Get("BetaCutProton_RPC_1.5");
TCutG* betamom_1sig_pim_tof_pionCmom=(TCutG*)cutfile_betamom_pionCmom->Get("BetaCutPiM_TOF_1.0");
TCutG* betamom_1sig_pim_rpc_pionCmom=(TCutG*)cutfile_betamom_pionCmom->Get("BetaCutPiM_RPC_1.0");
TCutG* betamom_1sig_pip_tof_pionCmom=(TCutG*)cutfile_betamom_pionCmom->Get("BetaCutPiP_TOF_1.0");
TCutG* betamom_1sig_pip_rpc_pionCmom=(TCutG*)cutfile_betamom_pionCmom->Get("BetaCutPiP_RPC_1.0");
TCutG* betamom_1sig_p_tof_pionCmom=(TCutG*)cutfile_betamom_pionCmom->Get("BetaCutProton_TOF_1.0");
TCutG* betamom_1sig_p_rpc_pionCmom=(TCutG*)cutfile_betamom_pionCmom->Get("BetaCutProton_RPC_1.0");

TFile * cutfile_betamom_simDeltaEle = new TFile("/lustre/hades/user/tscheib/apr12/ID_Cuts/BetaMomIDCuts_PionsProtons_gen8_SIM_newMDCdelta.root");
TCutG* betamom_1sig_pip_tof_simDeltaEle=(TCutG*)cutfile_betamom_simDeltaEle->Get("BetaCutPiP_TOF_1.0");
TCutG* betamom_1sig_pip_rpc_simDeltaEle=(TCutG*)cutfile_betamom_simDeltaEle->Get("BetaCutPiP_RPC_1.0");  //PiM
TCutG* betamom_1sig_pim_tof_simDeltaEle=(TCutG*)cutfile_betamom_simDeltaEle->Get("BetaCutPiM_TOF_1.0");
TCutG* betamom_1sig_pim_rpc_simDeltaEle=(TCutG*)cutfile_betamom_simDeltaEle->Get("BetaCutPiM_RPC_1.0");  //PiM
TCutG* betamom_1sig_p_tof_simDeltaEle=(TCutG*)cutfile_betamom_simDeltaEle->Get("BetaCutProton_TOF_1.0");
TCutG* betamom_1sig_p_rpc_simDeltaEle=(TCutG*)cutfile_betamom_simDeltaEle->Get("BetaCutProton_RPC_1.0");  //Proton
TCutG* betamom_15sig_pip_tof_simDeltaEle=(TCutG*)cutfile_betamom_simDeltaEle->Get("BetaCutPiP_TOF_1.5");
TCutG* betamom_15sig_pip_rpc_simDeltaEle=(TCutG*)cutfile_betamom_simDeltaEle->Get("BetaCutPiP_RPC_1.5");  //PiM
TCutG* betamom_15sig_pim_tof_simDeltaEle=(TCutG*)cutfile_betamom_simDeltaEle->Get("BetaCutPiM_TOF_1.5");
TCutG* betamom_15sig_pim_rpc_simDeltaEle=(TCutG*)cutfile_betamom_simDeltaEle->Get("BetaCutPiM_RPC_1.5");  //PiM
TCutG* betamom_15sig_p_tof_simDeltaEle=(TCutG*)cutfile_betamom_simDeltaEle->Get("BetaCutProton_TOF_1.5");
TCutG* betamom_15sig_p_rpc_simDeltaEle=(TCutG*)cutfile_betamom_simDeltaEle->Get("BetaCutProton_RPC_1.5");  //Proton
TCutG* betamom_2sig_pip_tof_simDeltaEle=(TCutG*)cutfile_betamom_simDeltaEle->Get("BetaCutPiP_TOF_2.0");
TCutG* betamom_2sig_pip_rpc_simDeltaEle=(TCutG*)cutfile_betamom_simDeltaEle->Get("BetaCutPiP_RPC_2.0");  //PiM
TCutG* betamom_2sig_pim_tof_simDeltaEle=(TCutG*)cutfile_betamom_simDeltaEle->Get("BetaCutPiM_TOF_2.0");
TCutG* betamom_2sig_pim_rpc_simDeltaEle=(TCutG*)cutfile_betamom_simDeltaEle->Get("BetaCutPiM_RPC_2.0");  //PiM
TCutG* betamom_2sig_p_tof_simDeltaEle=(TCutG*)cutfile_betamom_simDeltaEle->Get("BetaCutProton_TOF_2.0");
TCutG* betamom_2sig_p_rpc_simDeltaEle=(TCutG*)cutfile_betamom_simDeltaEle->Get("BetaCutProton_RPC_2.0");  //Proton
TCutG* betamom_25sig_pip_tof_simDeltaEle=(TCutG*)cutfile_betamom_simDeltaEle->Get("BetaCutPiP_TOF_2.5");
TCutG* betamom_25sig_pip_rpc_simDeltaEle=(TCutG*)cutfile_betamom_simDeltaEle->Get("BetaCutPiP_RPC_2.5");  //PiM
TCutG* betamom_25sig_pim_tof_simDeltaEle=(TCutG*)cutfile_betamom_simDeltaEle->Get("BetaCutPiM_TOF_2.5");
TCutG* betamom_25sig_pim_rpc_simDeltaEle=(TCutG*)cutfile_betamom_simDeltaEle->Get("BetaCutPiM_RPC_2.5");  //PiM
TCutG* betamom_25sig_p_tof_simDeltaEle=(TCutG*)cutfile_betamom_simDeltaEle->Get("BetaCutProton_TOF_2.5");
TCutG* betamom_25sig_p_rpc_simDeltaEle=(TCutG*)cutfile_betamom_simDeltaEle->Get("BetaCutProton_RPC_2.5");  //Proton
TCutG* betamom_3sig_pip_tof_simDeltaEle=(TCutG*)cutfile_betamom_simDeltaEle->Get("BetaCutPiP_TOF_3.0");
TCutG* betamom_3sig_pip_rpc_simDeltaEle=(TCutG*)cutfile_betamom_simDeltaEle->Get("BetaCutPiP_RPC_3.0");  //PiM
TCutG* betamom_3sig_pim_tof_simDeltaEle=(TCutG*)cutfile_betamom_simDeltaEle->Get("BetaCutPiM_TOF_3.0");
TCutG* betamom_3sig_pim_rpc_simDeltaEle=(TCutG*)cutfile_betamom_simDeltaEle->Get("BetaCutPiM_RPC_3.0");  //PiM
TCutG* betamom_3sig_p_tof_simDeltaEle=(TCutG*)cutfile_betamom_simDeltaEle->Get("BetaCutProton_TOF_3.0");
TCutG* betamom_3sig_p_rpc_simDeltaEle=(TCutG*)cutfile_betamom_simDeltaEle->Get("BetaCutProton_RPC_3.0");  //Proton
