#include "helper/TTGammaScaleFactor.h"
#include <utility>
#include <cmath>
#include <iostream>
#include "TString.h"

// Due to low stats in the TTGamma sample, we use the yield on the baseline region and 
// extrapolate to each Signal region using a scale factor derived from signal lepton
// ttbar yields.
//     N_ttgamma_SRn = N_ttgamma_sr0 * (N_ttbar_SRn / N_ttbar_SR0)
// where
//     N_ttbar_SRn   = # events with 1 lepton + 2 jets + SRn requirement
//     N_ttgamma_SR0 = baseline prediction for ttgamma from the analysis
// inputs:
//      sr_num        = SignalRegion Number: 0 - 28 (9, 19, 29 excluded) 
//      analysis_type = 0 --> high_pt anlysis, 1 --> low_pt analysis 
// returns pair<float, float> where 
//      first  = scale factor 
//      second = uncertainty 
// MC used: /TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM
std::pair<float, float> TTGammaScaleFactor(const unsigned int sr_num, const unsigned int analysis_type)
{
    // high pT results
    if (analysis_type == 0)
    {
        // number passing ttbar events for Signal Region 0
        const float N_ttbar_SR0 = 751250.62;

        // number passing ttbar events for Signal Region sr_num 
        float N_ttbar_SRn = -999999.0;
        switch (sr_num)
        {
            case 0  : N_ttbar_SRn = 751250.62 ; break;
            case 1  : N_ttbar_SRn = 135849.73 ; break;
            case 2  : N_ttbar_SRn = 11143.52  ; break;
            case 3  : N_ttbar_SRn = 72756.26  ; break;
            case 4  : N_ttbar_SRn = 38384.17  ; break;
            case 5  : N_ttbar_SRn = 36759.14  ; break;
            case 6  : N_ttbar_SRn = 7998.78   ; break;
            case 7  : N_ttbar_SRn = 12780.68  ; break;
            case 8  : N_ttbar_SRn = 16657.09  ; break;
            case 10 : N_ttbar_SRn = 372659.00 ; break;
            case 11 : N_ttbar_SRn = 67766.82  ; break;
            case 12 : N_ttbar_SRn = 5626.79   ; break;
            case 13 : N_ttbar_SRn = 33315.69  ; break;
            case 14 : N_ttbar_SRn = 17304.34  ; break;
            case 15 : N_ttbar_SRn = 18976.98  ; break;
            case 16 : N_ttbar_SRn = 4015.82   ; break;
            case 17 : N_ttbar_SRn = 6182.02   ; break;
            case 18 : N_ttbar_SRn = 7625.88   ; break;
            case 20 : N_ttbar_SRn = 207860.59 ; break;
            case 21 : N_ttbar_SRn = 39718.86  ; break;
            case 22 : N_ttbar_SRn = 3088.90   ; break;
            case 23 : N_ttbar_SRn = 27836.15  ; break;
            case 24 : N_ttbar_SRn = 15404.18  ; break;
            case 25 : N_ttbar_SRn = 7249.26   ; break;
            case 26 : N_ttbar_SRn = 1895.01   ; break;
            case 27 : N_ttbar_SRn = 4222.36   ; break;
            case 28 : N_ttbar_SRn = 6200.34   ; break;
            default : 
                std::cout << Form("[TTGammaScaleFactor] Signal Region %u not valid.  Return bogus value.", sr_num) << std::endl; 
                return std::make_pair(-999999.0, -999999.0); 
        }

        // return the scale factor 
        const float sf_value = N_ttbar_SRn / N_ttbar_SR0;
        const float sf_error = sqrt(sf_value * (1.0 - sf_value) / N_ttbar_SR0);
        return std::make_pair(sf_value, sf_error);
    }
    // low pT results
    else if (analysis_type == 1)
    {
        // number passing ttbar events for Signal Region 0
        const float N_ttbar_SR0 = 312383.12;

        // number passing ttbar events for Signal Region sr_num 
        float N_ttbar_SRn = -999999.0;
        switch (sr_num)
        {
            case 0  : N_ttbar_SRn = 312383.12 ; break; 
            case 1  : N_ttbar_SRn = 63861.72  ; break; 
            case 2  : N_ttbar_SRn = 11143.52  ; break; 
            case 3  : N_ttbar_SRn = 58841.65  ; break; 
            case 4  : N_ttbar_SRn = 38384.17  ; break; 
            case 5  : N_ttbar_SRn = 21909.01  ; break; 
            case 6  : N_ttbar_SRn = 7998.78   ; break; 
            case 7  : N_ttbar_SRn = 11577.39  ; break; 
            case 8  : N_ttbar_SRn = 16657.09  ; break; 
            case 10 : N_ttbar_SRn = 147874.94 ; break; 
            case 11 : N_ttbar_SRn = 31529.23  ; break; 
            case 12 : N_ttbar_SRn = 5626.79   ; break; 
            case 13 : N_ttbar_SRn = 26851.46  ; break; 
            case 14 : N_ttbar_SRn = 17304.34  ; break; 
            case 15 : N_ttbar_SRn = 11225.76  ; break; 
            case 16 : N_ttbar_SRn = 4015.82   ; break; 
            case 17 : N_ttbar_SRn = 5583.87   ; break; 
            case 18 : N_ttbar_SRn = 7625.88   ; break; 
            case 20 : N_ttbar_SRn = 108305.35 ; break; 
            case 21 : N_ttbar_SRn = 19509.91  ; break; 
            case 22 : N_ttbar_SRn = 3088.90   ; break; 
            case 23 : N_ttbar_SRn = 22842.32  ; break; 
            case 24 : N_ttbar_SRn = 15404.18  ; break; 
            case 25 : N_ttbar_SRn = 4641.87   ; break; 
            case 26 : N_ttbar_SRn = 1895.01   ; break; 
            case 27 : N_ttbar_SRn = 3844.01   ; break; 
            case 28 : N_ttbar_SRn = 6200.34   ; break; 
            default : 
                std::cout << Form("[TTGammaScaleFactor] Signal Region %u not valid.  Return bogus value.", sr_num) << std::endl; 
                return std::make_pair(-999999.0, -999999.0); 
        }

        // return the scale factor 
        const float sf_value = N_ttbar_SRn / N_ttbar_SR0;
        const float sf_error = sqrt(sf_value * (1.0 - sf_value) / N_ttbar_SR0);
        return std::make_pair(sf_value, sf_error);
    }
    else
    {
        std::cout << Form("[TTGammaScaleFactor] AnalysisType %u not valid.  Return bogus value.", analysis_type) << std::endl; 
        return std::make_pair(-999999.0, -999999.0); 
    }
} 

// apply the scale factor to the TTGamma prediciton for SR0
// inputs:
//      sr_num        = SignalRegion Number: 0 - 28 (9, 19, 29 excluded) 
//      analysis_type = 0 --> high_pt anlysis, 1 --> low_pt analysis 
//      N_ttg_SR0     = predictions for TTGamma sample in baseline 
// returns pair<float, float> where 
//      first  = predictions for TTGamma in SR sr_num 
//      second = uncertainty 
std::pair<float, float> ApplyTTGammaScaleFactor
(
    const unsigned int sr_num,
    const unsigned int analysis_type,
    const float N_ttg_SR0_value,
    const float N_ttg_SR0_error
)
{
    const std::pair<float, float> sf = TTGammaScaleFactor(sr_num, analysis_type);
    const float N_ttg_SRn_value = N_ttg_SR0_value * sf.first;
    const float N_ttg_SRn_error = N_ttg_SRn_value * sqrt(pow(sf.second/sf.first, 2) + pow(N_ttg_SR0_error/N_ttg_SR0_value, 2));

    if (N_ttg_SR0_value > 0)   return std::make_pair(N_ttg_SRn_value, N_ttg_SRn_error);
    else                       return std::make_pair(N_ttg_SR0_value, N_ttg_SR0_error);
      
}
