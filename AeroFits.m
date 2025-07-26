classdef AeroFits
    methods(Static)
        function [pCL, pCD, pCDi, pCD0, pLD, pEff, ...
                C_L_fit, C_D_fit, C_di_fit, C_d0_fit, L_D_fit, eff_fit] = ...
                AeroFitsFromAlpha(alphaArr, C_LArr, C_DArr, C_diArr, C_d0Arr, L_DArr, effArr)
            [pCL, C_L_fit] = AeroFits.FitfromAlpha(alphaArr, C_LArr);
            [pCD, C_D_fit] = AeroFits.FitfromAlpha(alphaArr, C_DArr);
            [pCDi, C_di_fit] = AeroFits.FitfromAlpha(alphaArr, C_diArr);
            [pCD0, C_d0_fit] = AeroFits.FitfromAlpha(alphaArr, C_d0Arr);
            [pLD, L_D_fit] = AeroFits.FitfromAlpha(alphaArr, L_DArr);
            [pEff, eff_fit] = AeroFits.FitfromAlpha(alphaArr, effArr);
        end
        
        function [pY, fit] = FitfromAlpha(alphaArr, yArr)
            pY = polyfit(alphaArr, yArr, 3);
            fit = polyval(pY, alphaArr);
        end
    end
end