  vector<string> nametriggers;
  vector<string> triggers;


if (isdata){

  long run=fTR->Run;

  if (run>=160404 && run <=161176){
    nametriggers.push_back("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v");
    nametriggers.push_back("HLT_Photon26_IsoVL_Photon18_IsoVL_v");
    nametriggers.push_back("HLT_Photon26_CaloIdL_IsoVL_Photon18_v");
    nametriggers.push_back("HLT_Photon26_IsoVL_Photon18_v");
  }
  else if (run>=161216 && run<=165208){
    nametriggers.push_back("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v");
    nametriggers.push_back("HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v");
    nametriggers.push_back("HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v");
    nametriggers.push_back("HLT_Photon20_R9Id_Photon18_R9Id_v");
    nametriggers.push_back("HLT_Photon26_CaloIdL_IsoVL_Photon18_v");
  }
  else if (run>=165970 && run<=173198){
    nametriggers.push_back("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v");
    nametriggers.push_back("HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v");
    nametriggers.push_back("HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v");
    nametriggers.push_back("HLT_Photon26_R9Id_Photon18_R9Id_v");
    nametriggers.push_back("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v");
    nametriggers.push_back("HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id_v");
    nametriggers.push_back("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v");
    nametriggers.push_back("HLT_Photon36_R9Id_Photon22_R9Id_v");
  }

 }


if (!isdata){

  // MC triggers from DiPhotonJets madgraph Summer11
nametriggers.push_back("HLT_Photon20_R9Id_Photon18_R9Id_v");
nametriggers.push_back("HLT_Photon26_Photon18_v");
nametriggers.push_back("HLT_Photon26_IsoVL_Photon18_v");
nametriggers.push_back("HLT_Photon26_IsoVL_Photon18_IsoVL_v");
nametriggers.push_back("HLT_Photon26_CaloIdL_IsoVL_Photon18_v");
nametriggers.push_back("HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v");
nametriggers.push_back("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v");
nametriggers.push_back("HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v");
nametriggers.push_back("HLT_Photon32_CaloIdL_Photon26_CaloIdL_v");
nametriggers.push_back("HLT_Photon36_CaloIdL_Photon22_CaloIdL_v");
nametriggers.push_back("HLT_DoublePhoton33_v");

 }

for (vector<string>::const_iterator it=nametriggers.begin(); it!=nametriggers.end(); it++){
  for (int j=1;j<10;j++) {
    TString a((*it));
    a+=j;
    triggers.push_back(a.Data());
  }
 }

  for (vector<string>::const_iterator it=triggers.begin(); it!=triggers.end(); it++){
    if ( GetHLTPrescale(*it)>1) {
      cout << "warning: using prescaled trigger!!! " << *it << " " << GetHLTPrescale(*it) << " in run " << fTR->Run << endl;
    }
    if ( GetHLTResult(*it) )        return true;
  }

 return false;



